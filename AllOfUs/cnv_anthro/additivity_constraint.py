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
import bisect
import json
import os
import numpy as np
import pandas as pd
from scipy import stats
from common import read_plink_cnv, design_matrix, rcg
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


def qualifying_calls(cnv, n, max_freq, exon_file, segdup_bed=None):
    """Boolean per call: rare (fewer than max_freq*n people carry a same-type call at >= 50% reciprocal overlap; the PGC frequency
    rule at a lower threshold) and gene-hitting by the model's rules (DEL: any exon by >= 1 bp; DUP: a whole transcript inside the call).
    With segdup_bed, calls whose midpoint lies in a segmental duplication do not qualify."""
    lim = max_freq * n
    site = cnv.groupby(["CHR", "BP1", "BP2", "TYPE"]).IID.agg(lambda s: set(s)).reset_index()
    site["rare"] = True
    for (c, t), g in site.groupby(["CHR", "TYPE"]):
        g = g.sort_values("BP1"); idx = g.index.to_numpy(); s = g.BP1.to_numpy(); e = g.BP2.to_numpy(); car = g.IID.to_list()
        maxlen = (e - s).max()
        for k in range(len(idx)):
            carriers = set()
            for m in range(bisect.bisect_left(s, s[k] - maxlen), bisect.bisect_right(s, e[k])):
                ov = min(e[k], e[m]) - max(s[k], s[m])
                if ov > 0 and ov >= 0.5 * (e[k] - s[k]) and ov >= 0.5 * (e[m] - s[m]):
                    carriers |= car[m]
                    if len(carriers) >= lim:
                        site.at[idx[k], "rare"] = False; break
    rare = cnv.merge(site[["CHR", "BP1", "BP2", "TYPE", "rare"]], on=["CHR", "BP1", "BP2", "TYPE"], how="left").rare.to_numpy(bool)
    U, ex, tx = load_genes(exon_file); cc = cnv.copy(); cc["row"] = np.arange(len(cc)); gene = np.zeros(len(cc), bool)
    for typ, feats, mode in [(1, ex, "overlap"), (3, tx, "within")]:
        m = (cc.TYPE == typ).to_numpy(); gene[m] = incidence(cc[m], feats, mode, len(cc), len(U)).getnnz(axis=1)[m] > 0
    q = rare & gene
    if segdup_bed:
        sd = pd.read_csv(segdup_bed, sep="\t", header=None, usecols=[0, 1, 2], comment="#")
        sd[0] = sd[0].astype(str).str.replace("chr", "", regex=False)
        mid = ((cc.BP1 + cc.BP2) // 2).to_numpy(); insd = np.zeros(len(cc), bool)
        for c, g in sd.groupby(0):
            if not c.isdigit():
                continue
            g = g.sort_values(1); st = g[1].to_numpy(); en = np.maximum.accumulate(g[2].to_numpy())
            m = (cc.CHR == int(c)).to_numpy(); i = np.searchsorted(st, mid[m], side="right") - 1
            insd[m] = (i >= 0) & (en[np.maximum(i, 0)] > mid[m])
        q &= ~insd
    return q


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
    ap.add_argument("--exclude_loci", default=None, help="TSV NAME CHR START END: drop carriers of these intervals (Molly's recurrent CNV genotyping)")
    ap.add_argument("--exclude_dosage", default="both", choices=["both", "DEL", "DUP"])
    ap.add_argument("--loci_mode", default="exclude", choices=["exclude", "only"],
                    help="exclude = drop carriers of --exclude_loci; only = keep only those carriers plus the --only_ref reference group")
    ap.add_argument("--only_ref", default="nocall", choices=["nocall", "no_rare_gene", "no_rare_gene_nosd"],
                    help="reference group of the 'only' runs: nocall (arrays), no_rare_gene (WGS), no_rare_gene_nosd (WGS sensitivity)")
    ap.add_argument("--rare_freq", type=float, default=0.01,
                    help="rare = carried by < this fraction of the sample (same-type calls at >= 50%% reciprocal overlap, the PGC frequency rule)")
    ap.add_argument("--segdup_bed", default="", help="segmental duplications BED (needed for --only_ref no_rare_gene_nosd)")
    ap.add_argument("--ref_exclude_loci", default="",
                    help="'only' runs with a no_rare_gene* reference: carriers of these intervals (either dosage) never enter the reference "
                         "(loci_sacks34.tsv), so the reference is the same in every run (Jon 2026-11-18)")
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
        G = rcg.genotype_calls(cnv)   # Molly's recurrent CNV genotyping (Jon 2026-11-24); counts below are carriers, not calls
        for l in L.itertuples():
            h_ = rcg.carriers(G, l.NAME, a.exclude_dosage)
            hit |= h_; excl[l.NAME] = len(h_)
        # exclude: drop carriers of these loci; only: keep only carriers of these loci plus the reference group (--only_ref):
        #   nocall              people with no CNV call at all (array definition; UKBB)
        #   no_rare_gene        people with no rare (< --rare_freq) gene-hitting call (WGS definition, Jon 2026-11-17: everyone
        #                       sequenced carries common CNVs, so "no call" is empty in AoU)
        #   no_rare_gene_nosd   people with no rare gene-hitting call whose midpoint lies outside segmental duplications (WGS sensitivity)
        if a.loci_mode == "exclude":
            drop = hit
        elif a.only_ref == "nocall":
            drop = set(cnv.IID) - hit
        else:
            q = qualifying_calls(cnv, len(ph), a.rare_freq, a.exons, a.segdup_bed if a.only_ref == "no_rare_gene_nosd" else None)
            drop = set(cnv.IID[q]) - hit
            if a.ref_exclude_loci:
                # Jon 2026-11-18: under _nosd a locus carrier whose call midpoint is in a segdup has no qualifying call, so it fell into
                # the reference of runs that do not include its locus (49 people in AoU). Carriers of any of these intervals stay out.
                R = pd.read_csv(a.ref_exclude_loci, sep="\t")
                R = R.rename(columns={c: c.upper() for c in R.columns if c.lower() in ("name", "chr", "start", "end")})
                R = R.drop_duplicates(["NAME", "CHR", "START", "END"]); rx = set()
                rx = rcg.carriers(G, R.NAME, "both")
                n_rx = len(rx - hit - drop)
                drop |= rx - hit
        ph = ph[~ph.IID.isin(drop)].reset_index(drop=True)
        Zc, _ = design_matrix(ph, num, fac); Q, _ = np.linalg.qr(Zc); N = len(ph)
        y = ph[a.trait].to_numpy(float); yr = y - Q @ (Q.T @ y)
        cnv = cnv[cnv.IID.isin(set(ph.IID))].reset_index(drop=True)
        excl = {"mode": a.loci_mode, "people_dropped": len(drop), "locus_carriers": len(hit), "dosage": a.exclude_dosage, "carriers_per_locus": excl}
        if a.loci_mode == "only":
            excl.update(only_ref=a.only_ref, reference_N=int(N - ph.IID.isin(hit).sum()), locus_carriers_in_sample=int(ph.IID.isin(hit).sum()))
            if a.only_ref != "nocall":
                excl.update(rare_freq=a.rare_freq, segdup_bed=os.path.basename(a.segdup_bed) if a.only_ref == "no_rare_gene_nosd" else None,
                            ref_exclude_loci=os.path.basename(a.ref_exclude_loci) or None,
                            ref_excluded_locus_carriers=n_rx if a.ref_exclude_loci else 0)
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
