"""Probe-level CNV GWAS of a quantitative trait -- Omar Shanta's PGC CNV GWAS (Freeze_2024) adapted to a linear model.

Kept from Omar's pipeline (scripts/ALL/<DIS>/noFam_noSex/cnv_gwas_scripts + META_noSex/meta_analysis_scripts):
  * probes = CNV breakpoints (BP1 and BP2 of every call in the freeze; master files chr<i>_{del,dup}.bp), restricted to
    probes overlapped by >= 1 call of that type in the analysed cohort (05_filter_breakpoint_probe_file.sh);
  * a person carries a probe if any of their calls of that type has BP1 <= probe <= BP2 (GenomicRanges overlap);
  * DEL (TYPE 1) and DUP (TYPE 3) tested separately; chr 1-23;
  * input calls already pass the freeze filters (10 kb, 10 probes, 10% frequency);
  * probe filters of 09_01_color.R: cohort carrier frequency > 2% (known loci exempt), chr11:134824089-134923878 DUP,
    NCNVs < 12, > 20 QC-filtered calls at the probe (optional input), segmental-duplication overlap;
  * genome-wide threshold 0.05 / n_independent_tests, independent tests by Omar's greedy Jaccard merge (>= 0.04) of the
    call sets of adjacent kept probes (calc_FWER/calculate_FWER_<chr>.sh).
Changed for a quantitative trait: logistic LRT -> OLS t-test of trait ~ carrier + covariates (exact, see common.py).
With several cohorts/platforms, run per cohort and combine with stouffer_n (sample-size weights, as Omar's meta_sa).
"""
import argparse
import heapq
import json
import os
import numpy as np
import pandas as pd
from common import read_plink_cnv, design_matrix, Residualiser

TYPES = {"DEL": 1, "DUP": 3}


def load_probes(bp_dir, cnv, tcode, chroms):
    """Omar's master breakpoint files if given, else the cohort's own call breakpoints (same rule, cohort-only)."""
    tname = "del" if tcode == 1 else "dup"
    out = []
    for c in chroms:
        if bp_dir:
            f = os.path.join(bp_dir, f"chr{c}_{tname}.bp")
            if not os.path.exists(f):
                continue
            p = pd.read_csv(f, sep="\t")["BP"].to_numpy(np.int64)
        else:
            s = cnv[(cnv.CHR == c) & (cnv.TYPE == tcode)]
            p = np.concatenate([s.BP1.to_numpy(), s.BP2.to_numpy()])
        out.append(pd.DataFrame({"CHR": c, "BP": np.unique(p)}))
    return pd.concat(out, ignore_index=True)


def sweep(cnv_t, probes_chr, row_of_iid, R, yr, yryr, keep_sets=None, set_cap=None):
    """Sweep one chromosome: incremental carrier sums -> exact OLS per probe. Optionally keep call-ID sets."""
    starts = cnv_t.sort_values("BP1")
    s_bp1 = starts.BP1.to_numpy(); s_row = starts["row"].to_numpy(); s_bp2 = starts.BP2.to_numpy(); s_id = starts.CNVID.to_numpy()
    k = R.Q.shape[1]
    qsum = np.zeros(k); ysum = 0.0; ncar = 0
    cnt = {}                        # person row -> number of active calls
    active = []                     # heap of (BP2, row, cnvid)
    active_ids = set()
    i = 0
    res = []
    df = R.N - R.k - 1
    n_since_refresh = 0
    for bp in probes_chr:
        while i < len(s_bp1) and s_bp1[i] <= bp:
            r = s_row[i]
            heapq.heappush(active, (s_bp2[i], r, s_id[i])); active_ids.add(s_id[i])
            c = cnt.get(r, 0)
            if c == 0:
                qsum += R.Q[r]; ysum += yr[r]; ncar += 1
            cnt[r] = c + 1
            i += 1
        while active and active[0][0] < bp:
            _, r, cid = heapq.heappop(active); active_ids.discard(cid)
            c = cnt[r] - 1
            if c == 0:
                del cnt[r]; qsum -= R.Q[r]; ysum -= yr[r]; ncar -= 1
            else:
                cnt[r] = c
        if keep_sets is not None and 0 < ncar and (set_cap is None or ncar <= set_cap):
            keep_sets[bp] = frozenset(active_ids)
        n_since_refresh += 1
        if n_since_refresh >= 20000 and ncar > 0:   # guard against floating drift
            rows = np.fromiter(cnt.keys(), int); qsum = R.Q[rows].sum(0); ysum = yr[rows].sum(); n_since_refresh = 0
        if ncar == 0:
            res.append((bp, 0, np.nan, np.nan, np.nan, 1.0)); continue
        xrxr = ncar - qsum @ qsum
        if xrxr <= 1e-10 or ncar == R.N:
            res.append((bp, ncar, np.nan, np.nan, np.nan, 1.0)); continue
        b = ysum / xrxr
        se = np.sqrt((yryr - b * ysum) / df / xrxr)
        t = b / se
        res.append((bp, ncar, b, se, t, np.nan))
    return res


def interval_hits(chrom, pos, beds):
    """Count bed intervals (chr,start,end) covering each position (start <= pos <= end)."""
    out = np.zeros(len(pos), int)
    for c in np.unique(chrom):
        m = chrom == c
        b = beds[beds[0] == c]
        if len(b) == 0:
            continue
        st = np.sort(b[1].to_numpy()); en = np.sort(b[2].to_numpy())
        p = pos[m]
        out[m] = np.searchsorted(st, p, side="right") - np.searchsorted(en, p, side="left")
    return out


def read_bed_any(path):
    b = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], dtype={0: str}, comment="#")
    b = b[b[0].str.match(r"^(chr)?[0-9XY]+$")]
    b[0] = b[0].str.replace("chr", "").replace({"X": "23", "Y": "24"}).astype(int)
    return b


def independent_tests(kept, sets):
    """Omar's calc_FWER greedy merge on one chromosome. kept: DataFrame (BP, P) sorted by BP."""
    bps = kept.BP.to_numpy(); ps = kept.P.to_numpy()
    n = len(bps)
    indep = np.ones(n, bool); group = np.zeros(n, int)
    if n == 0:
        return indep, group
    cur = 0; merge = False; merged = None; g = 0
    for j in range(n - 1):
        t1 = merged if merge else sets[bps[cur]]
        t2 = sets[bps[j + 1]]
        inter = len(t1 & t2); uni = len(t1 | t2)
        jac = np.floor(100 * inter / uni) / 100 if uni else 0.0    # bc scale=2 truncation, as in Omar's script
        if jac >= 0.04:
            merged = t1 | t2; merge = True; group[j + 1] = g
            if ps[cur] <= ps[j + 1]:
                indep[j + 1] = False
            else:
                indep[cur] = False; cur = j + 1
        else:
            cur = j + 1; merge = False; g += 1; group[j + 1] = g
    return indep, group


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno", required=True); ap.add_argument("--trait", required=True)
    ap.add_argument("--cnv", required=True); ap.add_argument("--out", required=True); ap.add_argument("--label", default="")
    ap.add_argument("--bp_dir", default=None, help="Omar's master breakpoint dir (chr<i>_{del,dup}.bp); omit to use cohort breakpoints")
    ap.add_argument("--numeric", default="age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10")
    ap.add_argument("--factors", default="sex,batch")
    ap.add_argument("--chroms", default="1-23")
    ap.add_argument("--known_loci", default=None); ap.add_argument("--segdups", default=None)
    ap.add_argument("--filtered_cnv_beds", default="", help="comma list of QC-filtered call beds (Omar's filteredCNVs)")
    ap.add_argument("--freq_max", type=float, default=0.02); ap.add_argument("--min_ncnv", type=int, default=12)
    ap.add_argument("--max_filtered", type=int, default=20)
    ap.add_argument("--dup_blacklist", default="11:134824089-134923878")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    lo, hi = map(int, a.chroms.split("-")); chroms = list(range(lo, hi + 1))

    ph = pd.read_csv(a.pheno, sep="\t", dtype={"IID": str})
    ph = ph[np.isfinite(ph[a.trait])].reset_index(drop=True)
    num = [c for c in a.numeric.split(",") if c]; fac = [c for c in a.factors.split(",") if c]
    ph = ph.dropna(subset=num + fac).reset_index(drop=True)
    Z, znames = design_matrix(ph, num, fac)
    R = Residualiser(Z)
    y = ph[a.trait].to_numpy(float)
    yr = R.resid(y); yryr = yr @ yr
    row_of_iid = pd.Series(np.arange(len(ph)), index=ph["IID"])
    cnv = read_plink_cnv(a.cnv, keep_iids=ph["IID"])
    cnv = cnv[cnv.CHR.isin(chroms)].copy()
    cnv["row"] = row_of_iid.loc[cnv["IID"]].to_numpy()
    info = {"label": a.label, "trait": a.trait, "N": int(len(ph)), "n_covariate_columns": int(Z.shape[1]),
            "covariates": znames[:20] + (["..."] if len(znames) > 20 else []),
            "calls_DEL": int((cnv.TYPE == 1).sum()), "calls_DUP": int((cnv.TYPE == 3).sum()),
            "carriers_any": int(cnv.IID.nunique())}
    known = pd.read_csv(a.known_loci, sep="\t") if a.known_loci else None
    seg = read_bed_any(a.segdups) if a.segdups else None
    filt = [read_bed_any(f) for f in a.filtered_cnv_beds.split(",") if f]
    from scipy import stats as st
    df_res = R.N - R.k - 1
    for tname, tcode in TYPES.items():
        pr = load_probes(a.bp_dir, cnv, tcode, chroms)
        rows = []; sets = {}
        for c in chroms:
            ct = cnv[(cnv.CHR == c) & (cnv.TYPE == tcode)]
            pc = pr[pr.CHR == c].BP.to_numpy()
            if len(ct) == 0 or len(pc) == 0:
                continue
            ks = {}
            r = sweep(ct, pc, row_of_iid, R, yr, yryr, keep_sets=ks, set_cap=int(0.05 * R.N))
            for bp, nc, b, se, t, p in r:
                rows.append((c, bp, nc, b, se, t))
            sets.update({(c, k): v for k, v in ks.items()})
        S = pd.DataFrame(rows, columns=["CHR", "BP", "NCNVs", "BETA", "SE", "T"])
        S = S[S.NCNVs > 0].reset_index(drop=True)          # Omar step 05: probes with >= 1 call in the cohort
        S["P"] = np.where(np.isfinite(S["T"]), 2 * st.t.sf(np.abs(S["T"]), df_res), 1.0)
        S["N"] = R.N; S["FREQ"] = S.NCNVs / R.N
        # ---- Omar's probe filters (09_01_color.R) ----
        fail = S.FREQ > a.freq_max
        if known is not None:
            for _, L in known.iterrows():
                fail &= ~((S.CHR == L.CHR) & (S.BP >= L.START) & (S.BP <= L.END))
        if tname == "DUP" and a.dup_blacklist:
            bc, rg = a.dup_blacklist.split(":"); b1, b2 = map(int, rg.split("-"))
            fail |= (S.CHR == int(bc)) & (S.BP >= b1) & (S.BP <= b2)
        fail |= S.NCNVs < a.min_ncnv
        if filt:
            nf = sum(interval_hits(S.CHR.to_numpy(), S.BP.to_numpy(), b) for b in filt)
            S["filtered_calls"] = nf; fail |= nf > a.max_filtered
        if seg is not None:
            S["segdup"] = interval_hits(S.CHR.to_numpy(), S.BP.to_numpy(), seg); fail |= S.segdup > 0
        S["KEEP"] = ~fail
        # ---- independent tests (Omar calc_FWER) ----
        S["INDEP"] = False; S["GROUP"] = -1
        goff = 0
        for c in chroms:
            k = S[(S.CHR == c) & S.KEEP].sort_values("BP")
            if len(k) == 0:
                continue
            ct = cnv[(cnv.CHR == c) & (cnv.TYPE == tcode)]
            cs = {bp: sets[(c, bp)] if (c, bp) in sets else frozenset(ct.CNVID[(ct.BP1 <= bp) & (ct.BP2 >= bp)]) for bp in k.BP}
            ind, grp = independent_tests(k[["BP", "P"]], cs)
            S.loc[k.index, "INDEP"] = ind; S.loc[k.index, "GROUP"] = grp + goff; goff += grp.max() + 1
        ntest = int(S.INDEP.sum()); thr = 0.05 / ntest if ntest else np.nan
        info[f"{tname}_probes_tested"] = int(len(S)); info[f"{tname}_probes_kept"] = int(S.KEEP.sum())
        info[f"{tname}_independent_tests"] = ntest; info[f"{tname}_gws_threshold"] = thr
        S.to_csv(os.path.join(a.out, f"sumstats_{tname.lower()}.tsv.gz"), sep="\t", index=False, float_format="%.6g")
        # ---- lead probes per merged group ----
        sig = S[S.KEEP & (S.P < 1e-4)]
        if len(sig):
            lead = sig.loc[sig.groupby("GROUP").P.idxmin()].copy()
            span = sig.groupby("GROUP").agg(START=("BP", "min"), END=("BP", "max"), n_probes=("BP", "size"))
            lead = lead.merge(span, left_on="GROUP", right_index=True).sort_values("P")
            lead["GWS"] = lead.P < thr; lead["TYPE"] = tname
            if known is not None:
                lead["known_locus"] = [";".join(known[(known.CHR == r.CHR) & (known.START <= r.BP) & (known.END >= r.BP)].NAME) for r in lead.itertuples()]
            lead.to_csv(os.path.join(a.out, f"lead_{tname.lower()}.tsv"), sep="\t", index=False, float_format="%.6g")
            info[f"{tname}_gws_groups"] = int(lead.GWS.sum())
    json.dump(info, open(os.path.join(a.out, "run_info.json"), "w"), indent=1, default=str)
    print(json.dumps(info, indent=1, default=str))


if __name__ == "__main__":
    main()
