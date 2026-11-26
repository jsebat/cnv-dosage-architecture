# Person-level ConGA (context-stratified gene-set analysis) of a quantitative trait: score test of each context gene set's CNV
# burden (genes of the set hit by the person's deletions, or duplications) beyond the total number of genes hit, against the
# genome-shift null with the lab's shared shifts; locus-removal runs with matched nulls. Portable (UK Biobank on Expanse, All of Us in
# the Workbench): every path comes in through arguments; run_pipeline.py step 'ctx' drives it from the config's "ctx" block.
# null model: genome shift of the gene annotation (set membership rotated around the refGene gene circle by the shared shifts of
#   resources/shared_shifts/<trait>_shifts.npy, restricted to >= 10 Mb), the same shifts in every cohort and every removal run
#   (replicate i = row i of the shift file everywhere, so cohorts can be meta-analysed replicate by replicate: meta_ctx.py)
"""Port of cnv_gwas/meta_aou_2026-11-24/conga_ctx_qt.py (the UK Biobank run behind Supplementary Fig. 1 of the CNV GWAS paper), with the
same statistics; only the defaults (Expanse paths) were removed. It is the quantitative-trait version of the ConGA-CNV fresh-run engine
(score_lib.R, score_sets: zE, the exact score test of the set count f in y ~ X0 + gene_count + f). For a linear model
    U = f' e,  V = s2 (f'f - |Q'f|^2),  z = U / sqrt(V)
with e the residual of y on [X0, gene_count], Q an orthonormal basis of [X0, gene_count], s2 = e'e / (N - rank). One fit per trait x
CNV type x run; every set and every null replicate is scored from it in gene space: f = G S, with G the person x gene incidence of the
CNV type (incidence() of conga_axes_qt.py: DEL = any exon >= 1 bp, DUP = a whole transcript) and S the gene x set membership.
Null: S is rotated around the gene circle (load_genes order) by each shared shift k; shifts shorter than the densest 10-Mb window are
dropped. Empirical two-sided P = (1 + #|z_null| >= |z_obs|) / (R + 1).
Locus-removal runs drop CNV calls before the fit:
  recurrent units: calls labelled by Molly Sacks's recurrent CNV genotyping (rcg.label_calls via conga_axes_qt.assign_loci) with the
                   unit's loci_sacks34 ids (both CNV types);
  span units (non-recurrent loci): calls overlapping the locus span (hg38) by >= 1 bp (both CNV types);
  minor_all: every remaining recurrent id together; no_34: every recurrent id.
Outputs (<out>/): conga_ctx_<runs>.tsv.gz (per set x run: z, z_emp, P_emp, n_genes, n_carriers; FDR), null_z_<run>.npy (R x sets,
float32, the score z of every null replicate; columns in the table's set order), info.json. n_carriers is a person count: the export
step drops it (export_aggregate.py), so only test statistics leave the Workbench.

  python conga_ctx_qt.py --pheno P --trait height --trait_col height_rint --cnv C --dosage DEL --sets <membership.tsv.gz> \
      --exons <refGene exons> --loci34 <loci_sacks34.tsv> --shifts <shared_shifts/{trait}_shifts.npy> --out <dir> \
      [--runs full | --runs loo --set_list <sets.txt> --loo_units <units.tsv>] [--numeric ... --factors ...]
"""
import argparse
import hashlib
import json
import os
import sys
import time
import numpy as np
import pandas as pd
from scipy import sparse

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from common import read_plink_cnv, design_matrix   # noqa: E402
from conga_axes_qt import load_genes, incidence, assign_loci   # noqa: E402


def build_runs(cnv, loci34, units):
    lab = np.asarray(assign_loci(cnv, loci34)).astype(str)
    ch, b1, b2 = cnv.CHR.to_numpy(), cnv.BP1.to_numpy(), cnv.BP2.to_numpy()
    runs = {"full": np.ones(len(cnv), bool)}; claimed = set()
    for u in units.itertuples():
        if u.kind == "recurrent":
            ids = u.ids.split("+"); claimed |= set(ids); drop = np.isin(lab, ids)
        else:
            c, se = u.ids.split(":"); s, e = map(int, se.split("-")); drop = (ch == int(c)) & (b1 <= e) & (b2 >= s)
        runs[f"loo_{u.run}"] = ~drop
    allids = sorted(set(loci34.id))
    runs["minor_all"] = ~np.isin(lab, [i for i in allids if i not in claimed])
    runs["no_34"] = ~np.isin(lab, allids)
    return runs, lab


class Fit:
    """Shared null fit of one trait x CNV type x run; scores sets for any rotation of the membership."""

    def __init__(self, y, Z, G):
        self.G = G.tocsr(); gc = np.asarray(self.G.sum(1)).ravel()
        W = np.column_stack([Z, gc]); Q, Rr = np.linalg.qr(W); keep = np.abs(np.diag(Rr)) > 1e-8 * np.abs(np.diag(Rr)).max()
        Q = Q[:, keep]; self.e = y - Q @ (Q.T @ y); self.s2 = float(self.e @ self.e) / (len(y) - Q.shape[1])
        rows = np.flatnonzero(gc > 0); Gc = self.G[rows].tocsr()                  # only carriers have f != 0
        # identical gene-hit patterns (recurrent CNVs) are scored once and weighted by their count
        key = pd.Series([Gc.indices[Gc.indptr[i]:Gc.indptr[i + 1]].tobytes() for i in range(Gc.shape[0])])
        code, uniq = pd.factorize(key); first = pd.Series(np.arange(len(code))).groupby(code).first().to_numpy()
        self.Gu = Gc[first].tocsr(); self.w = np.bincount(code).astype(np.float64)
        self.u = np.asarray(self.G.T @ self.e).ravel()                              # gene scores G'e
        self.QtG = np.asarray(self.G.T @ Q).T                                       # p x M = Q'G
        self.n_carriers = len(rows); self.n_patterns = self.Gu.shape[0]

    def score(self, S):
        """S: M x K dense membership (already rotated). Returns z, U, n carriers per set."""
        F = self.Gu @ S                                                             # patterns x K, dense
        ff = self.w @ (F * F)
        QF = self.QtG @ S                                                           # p x K = (Q'G) S
        U = self.u @ S
        V = self.s2 * (ff - (QF ** 2).sum(0))
        with np.errstate(invalid="ignore", divide="ignore"):
            z = np.where(V > 1e-8 * np.maximum(ff, 1e-12) * self.s2, U / np.sqrt(np.maximum(V, 1e-300)), np.nan)
        return z, U, (self.w @ (F > 0)).astype(int)


def rotate(S, k, M):
    """gene g takes the membership of gene (g - k) mod M"""
    return np.roll(S, int(k), axis=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno", required=True); ap.add_argument("--cnv", required=True)
    ap.add_argument("--trait", required=True, help="trait name used for the shift file ({trait}) and labels, e.g. height")
    ap.add_argument("--trait_col", default="", help="phenotype column (default <trait>_rint)")
    ap.add_argument("--dosage", required=True, choices=["DEL", "DUP"])
    ap.add_argument("--sets", required=True, help="TSV(.gz) set_id, entrez"); ap.add_argument("--out", required=True)
    ap.add_argument("--exons", required=True); ap.add_argument("--loci34", required=True)
    ap.add_argument("--shifts", required=True, help="shared shift file; may contain {trait}")
    ap.add_argument("--extra_rep", type=int, default=18000); ap.add_argument("--seed", type=int, default=20261124)
    ap.add_argument("--runs", default="full", choices=["full", "loo", "floor"]); ap.add_argument("--select", default="",
                    help="loo/floor: summary TSV of the full run (loo: FDR_BH < --loo_fdr; floor: sets at the replicate floor)")
    ap.add_argument("--loo_fdr", type=float, default=0.10); ap.add_argument("--loo_units", default="")
    ap.add_argument("--set_list", default="", help="loo: one set_id per line; overrides --select")
    ap.add_argument("--nrep", type=int, default=2000); ap.add_argument("--min_genes", type=int, default=10)
    ap.add_argument("--families", default="", help="TSV set_id, family (BH within trait x dosage x family)")
    ap.add_argument("--numeric", default="age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"); ap.add_argument("--factors", default="sex,batch")
    a = ap.parse_args(); t0 = time.time()
    if os.path.exists(os.path.join(a.out, "info.json")):
        raise SystemExit(f"{a.out} holds a finished run; choose a new output folder")
    os.makedirs(a.out, exist_ok=True)
    ph = pd.read_csv(a.pheno, sep="\t", dtype={"IID": str}); tr = a.trait_col or f"{a.trait}_rint"
    num = [x for x in a.numeric.split(",") if x]; fac = [x for x in a.factors.split(",") if x]
    ph = ph[np.isfinite(ph[tr])].dropna(subset=num + fac).reset_index(drop=True)
    Z, _ = design_matrix(ph, num, fac); y = ph[tr].to_numpy(float); N = len(y)
    cnv = read_plink_cnv(a.cnv, keep_iids=ph.IID); cnv = cnv[cnv.CHR <= 22].reset_index(drop=True)
    cnv["row"] = pd.Series(np.arange(N), index=ph.IID).loc[cnv.IID].to_numpy()
    U, ex, tx = load_genes(a.exons); M = len(U)
    ch = U.c.to_numpy(); pp = U.pos.to_numpy(); win = np.zeros(M, int)
    for c in np.unique(ch):
        idx = np.flatnonzero(ch == c); win[idx] = np.searchsorted(pp[idx], pp[idx] + 1e7, "left") - np.arange(len(idx))
    KMIN = int(win.max()) + 1
    shift_file = a.shifts.format(trait=a.trait)
    offs = np.load(shift_file).astype(np.int64); offs = offs[(offs >= KMIN) & (offs <= M - KMIN)][:a.nrep]
    if a.runs == "floor":   # sets at the replicate floor get fresh >= 10-Mb offsets; the same seed gives the same offsets in every cohort
        offs = np.random.default_rng(a.seed).integers(KMIN, M - KMIN + 1, a.extra_rep)
    sm = pd.read_csv(a.sets, sep="\t", dtype={"entrez": str})
    gi = dict(zip(U.entrez, U.gidx)); sm["g"] = sm.entrez.map(gi); sm = sm.dropna(subset=["g"]); sm["g"] = sm.g.astype(int)
    sm = sm.drop_duplicates(["set_id", "g"])
    n_on = sm.groupby("set_id").size(); keep = n_on[n_on >= a.min_genes].index
    if a.runs == "floor":
        sel = pd.read_csv(a.select, sep="\t"); sel = sel[(sel.trait == a.trait) & (sel.dosage == a.dosage) & (sel.run == "full")]
        sel = sel[sel.P_emp <= 1.0 / (sel.n_rep + 1) + 1e-12]
        keep = [s for s in keep if s in set(sel.set_id)]
        if not len(keep):
            json.dump(dict(note="no set at the floor"), open(f"{a.out}/info.json", "w")); print("no sets"); return
    if a.runs == "loo" and a.set_list:   # explicit list (every set that needs a robustness mark, both dosages)
        lst = set(pd.read_csv(a.set_list, sep="\t", header=None, dtype=str)[0])
        keep = [s for s in keep if s in lst]
    elif a.runs == "loo":
        sel = pd.read_csv(a.select, sep="\t"); sel = sel[(sel.trait == a.trait) & (sel.dosage == a.dosage) & (sel.FDR_BH < a.loo_fdr)]
        keep = [s for s in keep if s in set(sel.set_id)]
        if not len(keep):
            json.dump(dict(note="no set below loo_fdr"), open(f"{a.out}/info.json", "w")); print("no sets"); return
    sm = sm[sm.set_id.isin(keep)]; sets = sorted(sm.set_id.unique()); si = {s: i for i, s in enumerate(sets)}
    S = sparse.csr_matrix((np.ones(len(sm)), (sm.g.to_numpy(), sm.set_id.map(si).to_numpy())), shape=(M, len(sets))).toarray()
    loci34 = pd.read_csv(a.loci34, sep="\t")
    units = pd.read_csv(a.loo_units, sep="\t", dtype=str) if a.loo_units else pd.DataFrame(columns=["run", "kind", "ids"])
    runs, lab = build_runs(cnv, loci34, units)
    run_names = ["full"] if a.runs in ("full", "floor") else [r for r in runs if r != "full"]
    typ = 1 if a.dosage == "DEL" else 3
    info = dict(trait=a.trait, trait_col=tr, dosage=a.dosage, N=N, M=M, KMIN_10Mb=KMIN, nrep=int(len(offs)), sets=len(sets), runs=run_names,
                mode=a.runs, seed=a.seed if a.runs == "floor" else None, shift_file=os.path.basename(shift_file),
                shifts_sha1=hashlib.sha1(np.load(shift_file).astype(np.int64).tobytes()).hexdigest()[:16],
                offsets_sha1=hashlib.sha1(offs.astype(np.int64).tobytes()).hexdigest()[:16],
                calls_removed={r: int((~m).sum()) for r, m in runs.items()})
    out = []
    for rn in run_names:
        keepc = runs[rn] & (cnv.TYPE.to_numpy() == typ); cc = cnv[keepc]
        G = incidence(cc, ex if a.dosage == "DEL" else tx, "overlap" if a.dosage == "DEL" else "within", N, M)
        fit = Fit(y, Z, G)
        z0, U0, nc = fit.score(S)
        Zn = np.empty((len(offs), len(sets)), np.float32)
        print(rn, f"fit ready: carriers={fit.n_carriers} patterns={fit.n_patterns} {time.time() - t0:.0f}s", flush=True)
        for i, k in enumerate(offs):
            Zn[i] = fit.score(rotate(S, k, M))[0]
            if i in (9, 99) or (i + 1) % 500 == 0:
                print(rn, f"rep {i + 1}/{len(offs)} {time.time() - t0:.0f}s", flush=True)
        ge = (np.abs(Zn) >= np.abs(z0)[None, :] - 1e-12).sum(0); ok = np.isfinite(Zn).sum(0)
        mu = np.nanmean(Zn, 0); sd = np.nanstd(Zn, 0)
        d = pd.DataFrame(dict(trait=a.trait, dosage=a.dosage, run=rn, set_id=sets, n_genes=S.sum(0).astype(int),
                              n_carriers=nc, z=z0, z_emp=(z0 - mu) / np.where(sd > 0, sd, np.nan), P_emp=(1 + ge) / (ok + 1), n_rep=ok,
                              null_mean=mu, null_sd=sd, n_ge=ge))
        out.append(d); np.save(f"{a.out}/null_z_{rn}.npy", Zn)
        print(rn, f"carriers={fit.n_carriers}", f"min P_emp={d.P_emp.min():.2e}", f"{time.time() - t0:.0f}s", flush=True)
    D = pd.concat(out, ignore_index=True)
    if a.families:
        fam = pd.read_csv(a.families, sep="\t"); D = D.merge(fam, on="set_id", how="left")
    else:
        D["family"] = "all"
    from scipy.stats import false_discovery_control as fdc
    D["FDR_BH"] = np.nan
    for _, idx in D.groupby(["run", "family"]).groups.items():
        p = D.loc[idx, "P_emp"].to_numpy(); m = np.isfinite(p)
        if m.any():
            D.loc[np.asarray(idx)[m], "FDR_BH"] = fdc(p[m], method="bh")
    D.to_csv(f"{a.out}/conga_ctx_{a.runs}.tsv.gz", sep="\t", index=False, float_format="%.5g")
    info["seconds"] = round(time.time() - t0); json.dump(info, open(f"{a.out}/info.json", "w"), indent=1)


if __name__ == "__main__":
    main()
