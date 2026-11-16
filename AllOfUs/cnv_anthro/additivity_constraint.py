"""Is sub-additivity across CNVs explained by a concave dose-response in the constrained-gene (HC) count? (Jon 2026-11-15)

Per person and dosage: HC count = n(genes with s_het > 0.029) hit by any of their calls of that dosage (each gene once), total genes
hit (n_all), the unsigned TWAS |z|>3 count (consensus mirror term), and the number of gene-hitting calls.
1. Dose-response: covariate-residualised trait mean by HC-count bin (0, 1, 2, 3-4, 5-9, 10+) per dosage, in carriers of one vs >= 2
   gene-hitting CNVs (and the matched comparison of equal total HC count reached by one vs several CNVs).
2. Nested OLS (carriers + non-carriers, covariate-residualised): trait ~ n_all + mirror + HC (per arm) + multi + multi x HC (per arm)
     linear   HC counts
     concave  + presence of >= 1 HC gene + log(1 + HC)   (per arm)
   Does the multi x HC interaction survive once the HC dose-response may be concave?
Aggregate outputs only.
"""
import argparse
import json
import os
import numpy as np
import pandas as pd
from scipy import stats
from common import read_plink_cnv, design_matrix
from conga_axes_qt import load_genes, incidence
from conga_grid_qt import twas_z_both, measure

BINS = [(0, 0, "0"), (1, 1, "1"), (2, 2, "2"), (3, 4, "3-4"), (5, 9, "5-9"), (10, 10 ** 6, "10+")]


def ols(yr, X, Q, names):
    # Aliased columns (e.g. no hit-gene carriers on an arm in a restricted loci subset) are dropped and reported as NaN instead of
    # making A singular; a column is kept only if it raises the rank of the residualised design (in order).
    Xr = X - Q @ (Q.T @ X); keep = []
    for j in range(Xr.shape[1]):
        if np.linalg.matrix_rank(Xr[:, keep + [j]]) == len(keep) + 1:
            keep.append(j)
    K = Xr[:, keep]; A = K.T @ K; bk = np.linalg.solve(A, K.T @ yr)
    df = len(yr) - Q.shape[1] - len(keep); rss = yr @ yr - bk @ (K.T @ yr); sek = np.sqrt(np.diag(np.linalg.inv(A)) * rss / df)
    b = np.full(X.shape[1], np.nan); se = np.full(X.shape[1], np.nan); b[keep] = bk; se[keep] = sek
    return pd.DataFrame({"term": names, "beta": b, "SE": se, "T": b / se, "P": 2 * stats.t.sf(np.abs(b / se), df)}), rss, df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno", required=True); ap.add_argument("--trait", required=True); ap.add_argument("--cnv", required=True)
    ap.add_argument("--twas", required=True); ap.add_argument("--exons", required=True); ap.add_argument("--gene_cov", required=True)
    ap.add_argument("--clingen", required=True)
    ap.add_argument("--numeric", default="age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10")
    ap.add_argument("--factors", default="sex,batch"); ap.add_argument("--out", required=True)
    ap.add_argument("--hit_def", default="top5_hi", choices=["top5_hi", "phaplo", "union"],
                    help="single-hit genes: top-5%% s_het or ClinGen HI; Collins pHaplo >= 0.86; or the union of all three")
    ap.add_argument("--phaplo", default=None, help="Collins 2022 dosage sensitivity scores TSV (gene, pHaplo, pTriplo)")
    ap.add_argument("--exclude_loci", default=None, help="TSV NAME CHR START END: drop people carrying a call matching one of these intervals (50%% reciprocal overlap)")
    ap.add_argument("--exclude_dosage", default="both", choices=["both", "DEL", "DUP"])
    ap.add_argument("--loci_mode", default="exclude", choices=["exclude", "only"],
                    help="exclude = drop carriers of --exclude_loci; only = keep only those carriers (plus people with no CNV call)")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    ph = pd.read_csv(a.pheno, sep="\t", dtype={"IID": str})
    num = [c for c in a.numeric.split(",") if c]; fac = [c for c in a.factors.split(",") if c]
    ph = ph[np.isfinite(ph[a.trait])].dropna(subset=num + fac).reset_index(drop=True)
    Zc, _ = design_matrix(ph, num, fac); Q, _ = np.linalg.qr(Zc); N = len(ph)
    y = ph[a.trait].to_numpy(float); yr = y - Q @ (Q.T @ y)
    cnv = read_plink_cnv(a.cnv, keep_iids=ph.IID); cnv = cnv[cnv.CHR <= 22].reset_index(drop=True)
    excl = {}
    if a.exclude_loci:   # 2026-11-15 (Jon): asymmetric (monogenic-like) CNVs out -> is the rest additive?
        L = pd.read_csv(a.exclude_loci, sep="\t"); hit = set()
        L = L.rename(columns={c: c.upper() for c in L.columns if c.lower() in ("name", "chr", "start", "end")})   # loci_sacks34.tsv is lower-case, one row per dosage
        L = L.drop_duplicates(["NAME", "CHR", "START", "END"]).reset_index(drop=True)
        for l in L.itertuples():
            c_ = cnv[(cnv.CHR == l.CHR) & (cnv.BP2 >= l.START) & (cnv.BP1 <= l.END)]
            if a.exclude_dosage != "both":
                c_ = c_[c_.TYPE == (1 if a.exclude_dosage == "DEL" else 3)]
            ov = (np.minimum(c_.BP2, l.END) - np.maximum(c_.BP1, l.START)).clip(lower=0)
            m = (ov >= 0.5 * (l.END - l.START)) & (ov >= 0.5 * (c_.BP2 - c_.BP1))
            hit |= set(c_.IID[m]); excl[l.NAME] = int(m.sum())
        # exclude: drop carriers of these loci; only: keep only carriers of these loci plus people with no CNV call at all
        drop = hit if a.loci_mode == "exclude" else set(cnv.IID) - hit
        ph = ph[~ph.IID.isin(drop)].reset_index(drop=True)
        Zc, _ = design_matrix(ph, num, fac); Q, _ = np.linalg.qr(Zc); N = len(ph)
        y = ph[a.trait].to_numpy(float); yr = y - Q @ (Q.T @ y)
        cnv = cnv[cnv.IID.isin(set(ph.IID))].reset_index(drop=True)
        excl = {"mode": a.loci_mode, "people_dropped": len(drop), "locus_carriers": len(hit), "dosage": a.exclude_dosage, "calls_matched_per_locus": excl}
    cnv["row"] = pd.Series(np.arange(N), index=ph.IID).loc[cnv.IID].to_numpy()
    U, ex, tx = load_genes(a.exons); M = len(U)
    gc = pd.read_csv(a.gene_cov, sep="\t", dtype={"entrez": str}); gc.loc[gc.shet_imputed == 1, "shet"] = np.nan
    shet = U[["entrez"]].merge(gc[["entrez", "shet"]], on="entrez", how="left").shet.to_numpy()
    Z, _ = twas_z_both(a.twas, "gtex_nonbrain", U, gc)
    cc_ = cnv.copy(); cc_["row"] = np.arange(len(cc_)); ng = np.zeros(len(cnv))
    for typ, feats, mode in [(1, ex, "overlap"), (3, tx, "within")]:
        m = (cnv.TYPE == typ).to_numpy(); ng[m] = incidence(cc_[m], feats, mode, len(cnv), M).getnnz(axis=1)[m]
    n_calls = np.bincount(cnv.row[ng > 0], minlength=N)
    Gd = incidence(cnv[cnv.TYPE == 1], ex, "overlap", N, M); Gu = incidence(cnv[cnv.TYPE == 3], tx, "within", N, M)
    wH = measure("none_H", Z, shet, 0, np.zeros(M, bool)); wU = measure("g:all:c3:uns:scaled", Z, shet, 0, np.zeros(M, bool))
    f = {}
    for d, G in [("DEL", Gd), ("DUP", Gu)]:
        f[f"nall_{d}"] = np.asarray(G.sum(1)).ravel(); f[f"HC_{d}"] = np.asarray(G @ wH).ravel(); f[f"mir_{d}"] = np.asarray(G @ wU).ravel()
    multi = (n_calls >= 2).astype(float); single = n_calls == 1
    info = dict(trait=a.trait, N=int(N), single_carriers=int(single.sum()), multi_carriers=int((n_calls >= 2).sum()), exclusion=excl)
    # 1. dose-response bins
    rows = []
    for d in ["DEL", "DUP", "total"]:
        h = f["HC_DEL"] + f["HC_DUP"] if d == "total" else f[f"HC_{d}"]
        for grp, mask in [("single CNV", single), ("2+ CNVs", n_calls >= 2)]:
            for lo, hi, lab in BINS:
                m = mask & (h >= lo) & (h <= hi)
                if m.sum() >= 1:
                    v = yr[m]; rows.append(dict(dosage=d, group=grp, HC_bin=lab, n=int(m.sum()), mean_SD=v.mean(),
                                                SE=v.std(ddof=1) / np.sqrt(m.sum()) if m.sum() > 1 else np.nan))
    pd.DataFrame(rows).to_csv(os.path.join(a.out, "dose_response_bins.tsv"), sep="\t", index=False, float_format="%.6g")
    # 2. nested models
    base = ["nall_DEL", "nall_DUP", "mir_DEL", "mir_DUP"]
    Xb = [f[c] for c in base]
    lin = Xb + [f["HC_DEL"], f["HC_DUP"], multi, multi * f["HC_DEL"], multi * f["HC_DUP"]]
    nl = base + ["HC_DEL", "HC_DUP", "multi", "multi_x_HC_DEL", "multi_x_HC_DUP"]
    T1, rss1, df1 = ols(yr, np.column_stack(lin), Q, nl)
    con = lin + [(f["HC_DEL"] > 0).astype(float), (f["HC_DUP"] > 0).astype(float), np.log1p(f["HC_DEL"]), np.log1p(f["HC_DUP"])]
    nc = nl + ["anyHC_DEL", "anyHC_DUP", "log1pHC_DEL", "log1pHC_DUP"]
    T2, rss2, df2 = ols(yr, np.column_stack(con), Q, nc)
    k12 = max(df1 - df2, 1); Fc = ((rss1 - rss2) / k12) / (rss2 / df2); P_concave = stats.f.sf(Fc, k12, df2)   # k12 = 4 unless terms aliased
    # same models without the multi terms: does concavity alone fit?
    T0, rss0, df0 = ols(yr, np.column_stack(Xb + [f["HC_DEL"], f["HC_DUP"]]), Q, base + ["HC_DEL", "HC_DUP"])
    T0c, rss0c, df0c = ols(yr, np.column_stack(Xb + [f["HC_DEL"], f["HC_DUP"], (f["HC_DEL"] > 0).astype(float), (f["HC_DUP"] > 0).astype(float),
                                                    np.log1p(f["HC_DEL"]), np.log1p(f["HC_DUP"])]), Q,
                           base + ["HC_DEL", "HC_DUP", "anyHC_DEL", "anyHC_DUP", "log1pHC_DEL", "log1pHC_DUP"])
    k0 = max(df0 - df0c, 1); F0 = ((rss0 - rss0c) / k0) / (rss0c / df0c)
    # 3. single hit vs constraint count (Jon 2026-11-15): is the duplication-axis constraint effect a single hit? Hit = a top-5% s_het
    #    or ClinGen haploinsufficient (score 2-3) gene on that arm, as in the locus-level single-hit test (yaxis_monogenic, top5_or_hi).
    q95 = float(np.nanquantile(gc.shet.to_numpy(float), 0.95))
    cg = pd.read_csv(a.clingen, sep="\t", skiprows=5, dtype=str, quoting=3); cg.columns = [c.strip().lstrip("#") for c in cg.columns]
    hi = U.symbol.isin(set(cg.loc[cg["Haploinsufficiency Score"].isin(["2", "3"]), cg.columns[0]])).to_numpy()
    top5hi = (np.nan_to_num(shet, nan=-1) >= q95) | hi
    if a.hit_def == "top5_hi":
        wHit = top5hi.astype(float)
    else:                                                       # Collins et al. 2022 pHaplo >= 0.86 (the paper's cut-off)
        co = pd.read_csv(a.phaplo, sep="\t"); co.columns = [c.lstrip("#") for c in co.columns]
        ph = U.symbol.isin(set(co.loc[co.pHaplo >= 0.86, "gene"])).to_numpy()
        wHit = (ph if a.hit_def == "phaplo" else (ph | top5hi)).astype(float)
    info.update(hit_def=a.hit_def, hit_genes=int(wHit.sum()))
    for d, G in [("DEL", Gd), ("DUP", Gu)]:
        f[f"hit_{d}"] = (np.asarray(G @ wHit).ravel() > 0).astype(float)
        f[f"HCnohit_{d}"] = np.asarray(G @ (wH * (1 - wHit))).ravel()      # constrained genes that are not hit genes
    hb = base + ["HC_DEL", "HC_DUP"]
    Th0, _, _ = ols(yr, np.column_stack([f[c] for c in hb]), Q, hb)
    Th1, _, _ = ols(yr, np.column_stack([f[c] for c in hb + ["hit_DEL", "hit_DUP"]]), Q, hb + ["hit_DEL", "hit_DUP"])
    hn = base + ["HCnohit_DEL", "HCnohit_DUP", "hit_DEL", "hit_DUP"]
    Th2, _, _ = ols(yr, np.column_stack([f[c] for c in hn]), Q, hn)
    # 4. two-component check (Jon 2026-11-15): do the mirror (TWAS) and constraint slopes differ in carriers of a hit-gene CNV on that arm?
    for d in ["DEL", "DUP"]:
        f[f"mirXhit_{d}"] = f[f"mir_{d}"] * f[f"hit_{d}"]; f[f"HCXhit_{d}"] = f[f"HC_{d}"] * f[f"hit_{d}"]
    hx = hb + ["hit_DEL", "hit_DUP", "mirXhit_DEL", "mirXhit_DUP", "HCXhit_DEL", "HCXhit_DUP"]
    Th3, _, _ = ols(yr, np.column_stack([f[c] for c in hx]), Q, hx)
    info.update(hit_carriers_DEL=int(f["hit_DEL"].sum()), hit_carriers_DUP=int(f["hit_DUP"].sum()), q95_shet=q95)
    out = pd.concat([T1.assign(model="linear HC + multi interactions"), T2.assign(model="concave HC + multi interactions"),
                     T0.assign(model="linear HC, no multi terms"), T0c.assign(model="concave HC, no multi terms"),
                     Th0.assign(model="hit: HC count only"), Th1.assign(model="hit: HC count + hit presence"),
                     Th2.assign(model="hit: HC count without hit genes + hit presence"),
                     Th3.assign(model="hit: slopes by hit status (mirror x hit, HC x hit)")])
    out.to_csv(os.path.join(a.out, "additivity_models.tsv"), sep="\t", index=False, float_format="%.6g")
    info.update(F_concavity_given_multi=float(Fc), P_concavity_given_multi=float(P_concave),
                F_concavity_no_multi=float(F0), P_concavity_no_multi=float(stats.f.sf(F0, k0, df0c)))
    json.dump(info, open(os.path.join(a.out, "run_info.json"), "w"), indent=1)
    pd.set_option("display.width", 200)
    print(json.dumps(info, indent=1)); print(out[out.term.str.contains("HC|multi")].to_string())
    print(pd.DataFrame(rows).to_string())


if __name__ == "__main__":
    main()
