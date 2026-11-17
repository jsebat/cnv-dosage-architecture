"""Test 4 (Jon 2026-11-16): are the 34 recurrent intervals enriched for TWAS genes beyond their gene count and their constraint?

The recurrent intervals are mostly developmental-disorder loci, so they are rich in constrained genes, and constrained genes are often TWAS
genes. Two nulls separate "TWAS genes per se" (4a) from "constraint" (4b):
  genes       random windows of the same number of consecutive genes (same chromosome, not overlapping any of the 34 intervals);
  genes+HC    as 'genes', and the same number of constrained genes (s_het > 0.029, GeneBayes, imputed values excluded).
Statistic: number of TWAS genes (|z| > cut) in the union of the intervals' genes, per TWAS measure; also among constrained genes only (HC).
Gene rule: a gene belongs to an interval if any exon overlaps it. Genes with >= 50% of their exonic bases in segmental duplications are
dropped everywhere (as in the size-bin grid). Public annotations only.

  python twas_enrichment_intervals.py --exons ... --gene_cov ... --twas_height ... --twas_bmi ... --loci loci_sacks34.tsv --segdup_bed ... --out <dir>
"""
import argparse
import json
import os
import numpy as np
import pandas as pd
from conga_axes_qt import load_genes, SHET_CUT
from conga_grid_qt import twas_z_both, segdup_fraction


def main():
    ap = argparse.ArgumentParser()
    for k in ["exons", "gene_cov", "twas_height", "twas_bmi", "loci", "segdup_bed", "out"]:
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--segdup_frac", type=float, default=0.5); ap.add_argument("--nsets", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--loo_runs", default="", help="leave-one-locus-out units (loo_runs format: run, del_members, dup_members = interval ids); each unit's intervals are dropped from the observed count and from the null (same draws, so the null is matched)")
    ap.add_argument("--extra_runs", default="", help="extra removal units in the same format (e.g. all of 16p11.2), reported but not counted in the LOO tally")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(a.seed)
    U, ex, tx = load_genes(a.exons); M = len(U)
    gc = pd.read_csv(a.gene_cov, sep="\t", dtype={"entrez": str}); gc.loc[gc.shet_imputed == 1, "shet"] = np.nan
    shet = U[["entrez"]].merge(gc[["entrez", "shet"]], on="entrez", how="left").shet.to_numpy()
    hc = np.nan_to_num(shet, nan=-1) > SHET_CUT
    keep = segdup_fraction(ex, M, a.segdup_bed) < a.segdup_frac
    # gene order on each chromosome (the gene circle of load_genes), kept genes only
    U["keep"] = keep; G = U[U.keep].reset_index(drop=True)
    gidx = G.gidx.to_numpy(); chrom = G.c.to_numpy()
    # interval genes (any exon overlap), union over dosage rows with the same coordinates
    L = pd.read_csv(a.loci, sep="\t"); print("loci columns:", list(L.columns))
    L = L.drop_duplicates(["chr", "start", "end"]).reset_index(drop=True)
    namecol = next((c for c in ["locus", "name", "unit", "id"] if c in L.columns), L.columns[0])
    exk = ex[ex.gidx.isin(set(gidx))]
    sets, blocked = [], np.zeros(M, bool)
    for r in L.itertuples():
        m = (exk.c == int(r.chr)) & (exk.start <= r.end) & (exk.end >= r.start)
        gs = np.unique(exk.gidx[m].to_numpy().astype(int)); sets.append(gs)
        bm = (U.c.to_numpy() == int(r.chr)) & (U.pos.to_numpy() <= r.end) & (U.pos.to_numpy() >= r.start - 2_000_000)   # generous block
        blocked |= bm
    for gs in sets:
        blocked[gs] = True
    pos_in_G = {g: i for i, g in enumerate(gidx)}
    # windows: start index i, length k, same chromosome, no blocked gene
    blockedG = blocked[gidx]
    cs_block = np.r_[0, np.cumsum(blockedG)]; cs_hc = np.r_[0, np.cumsum(hc[gidx])]
    chr_end = {c: np.flatnonzero(chrom == c).max() for c in np.unique(chrom)}
    def candidates(k, nhc):
        if k == 0:
            z0 = np.array([0]); return z0, z0
        i = np.arange(len(gidx) - k + 1)
        ok = (chrom[i] == chrom[i + k - 1]) & (cs_block[i + k] - cs_block[i] == 0)
        okhc = ok & (cs_hc[i + k] - cs_hc[i] == nhc) if nhc is not None else None
        return i[ok], (i[okhc] if okhc is not None else None)
    ks = [len(s) for s in sets]; nhcs = [int(hc[s].sum()) for s in sets]
    cand = [candidates(k, h) for k, h in zip(ks, nhcs)]
    info = dict(intervals=len(sets), genes_in_intervals=int(len(np.unique(np.concatenate(sets)))), genes_kept=int(keep.sum()), genes_total=M,
                genes_excluded_segdup=int((~keep).sum()), interval_gene_counts=dict(zip(L[namecol].astype(str), ks)),
                interval_hc_counts=dict(zip(L[namecol].astype(str), nhcs)),
                n_candidate_windows_min_genes=int(min(len(c[0]) for c in cand)), n_candidate_windows_min_genes_hc=int(min(len(c[1]) for c in cand)))
    # if an interval has no window with exactly its HC count, relax to +-1, then +-2 (recorded)
    relaxed = {}
    for j, (k, h) in enumerate(zip(ks, nhcs)):
        d = 0
        while len(cand[j][1]) == 0 and d < 3:
            d += 1; i = np.arange(len(gidx) - k + 1)
            ok = (chrom[i] == chrom[i + k - 1]) & (cs_block[i + k] - cs_block[i] == 0) & (np.abs(cs_hc[i + k] - cs_hc[i] - h) <= d)
            cand[j] = (cand[j][0], i[ok])
        if d:
            relaxed[str(L[namecol][j])] = d
    info["hc_match_relaxed"] = relaxed
    ids = L["id"].astype(str).to_numpy() if "id" in L.columns else L[namecol].astype(str).to_numpy()
    units = []      # (kind, run, boolean mask of intervals kept)
    for kind, fn in [("loo", a.loo_runs), ("extra", a.extra_runs)]:
        if fn:
            for r in pd.read_csv(fn, sep="	", keep_default_na=False).itertuples():
                mem = set(x for x in (str(r.del_members) + "," + str(r.dup_members)).split(",") if x)
                drop = np.isin(ids, list(mem))
                if drop.any():
                    units.append((kind, str(r.run), ~drop))
    info["loo_units"] = {u[1]: [str(x) for x in ids[~u[2]]] for u in units}
    rows, loo_rows = [], []
    allg = np.unique(np.concatenate(sets))
    for trait, tf in [("height", a.twas_height), ("bmi", a.twas_bmi)]:
        for panels in (["gtex_nonbrain"] if trait == "height" else ["gtex_nonbrain", "all"]):
            Z, npan = twas_z_both(tf, panels, U, gc)
            for tissue in ["scaled", "mean"]:
                z = np.nan_to_num(Z[tissue], nan=0.0)
                for cut in [3.0, 4.0]:
                    tw = np.abs(z) > cut
                    for genes in ["all", "HC"]:
                        w = (tw & hc) if genes == "HC" else tw
                        obs = int(w[allg].sum()); per = [int(w[s].sum()) for s in sets]
                        for null in ["genes", "genes+HC"]:
                            csw = np.r_[0, np.cumsum(w[gidx])]; S = np.zeros((len(ks), a.nsets))
                            for j, k in enumerate(ks):
                                c = cand[j][0] if null == "genes" else cand[j][1]
                                ws = csw[c + k] - csw[c]; S[j] = ws[rng.integers(len(c), size=a.nsets)]
                            sims = S.sum(0)
                            for kind, run, kp in units:     # matched LOO null: the same window draws without the dropped intervals
                                so = int(np.asarray(per)[kp].sum()); sl = S[kp].sum(0); mu_, sd_ = sl.mean(), sl.std(ddof=1)
                                loo_rows.append(dict(trait=trait, panels=panels, tissue=tissue, cut=cut, genes=genes, null=null, kind=kind, run=run,
                                                     observed_sum=so, null_mean=mu_, null_sd=sd_, ratio=so / mu_ if mu_ > 0 else np.nan,
                                                     z=(so - mu_) / sd_ if sd_ > 0 else np.nan, P_upper=(1 + np.sum(sl >= so)) / (a.nsets + 1)))
                            # observed uses the union (overlapping intervals counted once); null windows are drawn independently,
                            # so also report the observed sum over intervals (same counting as the null)
                            osum = int(sum(per))
                            rows.append(dict(trait=trait, panels=panels, tissue=tissue, cut=cut, genes=genes, null=null,
                                             observed_union=obs, observed_sum=osum, null_mean=sims.mean(), null_sd=sims.std(ddof=1),
                                             ratio=osum / sims.mean() if sims.mean() > 0 else np.nan,
                                             z=(osum - sims.mean()) / sims.std(ddof=1) if sims.std(ddof=1) > 0 else np.nan,
                                             P_upper=(1 + np.sum(sims >= osum)) / (a.nsets + 1), n_panels=npan))
                            print(rows[-1], flush=True)
                        # per-interval counts (for the table)
                        if null == "genes+HC" and cut == 3.0 and genes == "HC":
                            pd.DataFrame(dict(interval=L[namecol].astype(str), n_genes=ks, n_hc=nhcs, n_twas=per)).to_csv(
                                os.path.join(a.out, f"per_interval_{trait}_{panels}_{tissue}_HC_c3.tsv"), sep="\t", index=False)
    pd.DataFrame(rows).to_csv(os.path.join(a.out, "twas_enrichment_intervals.tsv"), sep="\t", index=False, float_format="%.5g")
    if loo_rows:
        LR = pd.DataFrame(loo_rows); LR.to_csv(os.path.join(a.out, "twas_enrichment_loo.tsv"), sep="	", index=False, float_format="%.5g")
        R = pd.DataFrame(rows); key = ["trait", "panels", "tissue", "cut", "genes", "null"]
        lo = LR[LR.kind == "loo"]
        g = lo.groupby(key).agg(loo_n=("P_upper", "size"), loo_held=("P_upper", lambda p: int((p < 0.05).sum())), loo_Pmax=("P_upper", "max"),
                                weakest=("P_upper", lambda p: lo.loc[p.idxmax(), "run"])).reset_index()
        for run in LR[LR.kind == "extra"].run.unique():
            x = LR[(LR.kind == "extra") & (LR.run == run)][key + ["observed_sum", "null_mean", "z", "P_upper"]]
            x.columns = key + [f"x_{run}_{c}" for c in ["observed", "null_mean", "z", "P"]]; g = g.merge(x, on=key, how="left")
        R.merge(g, on=key, how="left").to_csv(os.path.join(a.out, "twas_enrichment_intervals_loo_summary.tsv"), sep="	", index=False, float_format="%.5g")
    json.dump(info, open(os.path.join(a.out, "run_info.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
