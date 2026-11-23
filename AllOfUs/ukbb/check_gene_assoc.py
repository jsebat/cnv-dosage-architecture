"""Checks of the gene-based association output (run where the individual data are; prints aggregate numbers only).
1. For 20 genes per dosage (spread over carrier counts), refit trait ~ covariates + carrier indicator directly (numpy lstsq on the full design,
   carriers recomputed from the CNV table with the standard rules) and compare beta, SE and t with gene_assoc.tsv.gz.
2. One-cohort meta = the exported table. 3. Two-copy meta (a self-test, not a result): where both copies have the gene, beta equal and SE / sqrt2.
null model: none (analytic OLS).
  python check_gene_assoc.py --config <pipeline config> --table <results>/all/gene_assoc/<trait>/gene_assoc.tsv.gz --trait height --meta_one D --meta_two D
"""
import argparse, json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cnv_anthro"))
from common import read_plink_cnv, design_matrix
from conga_axes_qt import load_genes, incidence

ap = argparse.ArgumentParser(); ap.add_argument("--config", required=True); ap.add_argument("--table", required=True); ap.add_argument("--trait", required=True)
ap.add_argument("--meta_one", default=""); ap.add_argument("--meta_two", default=""); ap.add_argument("--cnv", default="", help="call table of --table if not the config one (e.g. cnv_rare.cnv)")
a = ap.parse_args(); C = json.load(open(a.config)); T = C["traits"][a.trait]; bad = 0
R = pd.read_csv(a.table, sep="\t", dtype={"entrez": str})
ph = pd.read_csv(C["pheno_table"], sep="\t", dtype={"IID": str})
ph = ph[np.isfinite(ph[T["column"]])].dropna(subset=T["numeric"] + T["factors"]).reset_index(drop=True); N = len(ph)
Zc, _ = design_matrix(ph, T["numeric"], T["factors"]); y = ph[T["column"]].to_numpy(float)
cnv = read_plink_cnv(a.cnv or C["cnv_table"], keep_iids=ph.IID); cnv = cnv[cnv.CHR <= 22].reset_index(drop=True)
cnv["row"] = pd.Series(np.arange(N), index=ph.IID).loc[cnv.IID].to_numpy()
U, ex, tx = load_genes(C["grid"]["exons"]); M = len(U); gi = dict(zip(U.entrez.astype(str), U.gidx))
worst = 0.0
for dos, typ, feats, mode in (("DEL", 1, ex, "overlap"), ("DUP", 3, tx, "within")):
    G = incidence(cnv[cnv.TYPE == typ], feats, mode, N, M).tocsc()
    d = R[(R.dosage == dos) & (R.n_carriers >= 5)].sort_values("n_carriers"); pick = d.iloc[np.linspace(0, len(d) - 1, 20).astype(int)]
    for r in pick.itertuples():
        x = np.asarray(G[:, gi[r.entrez]].todense()).ravel(); X = np.column_stack([Zc, x])
        b, res, rk, _ = np.linalg.lstsq(X, y, rcond=None); e = y - X @ b; s2 = e @ e / (N - X.shape[1])
        se = np.sqrt(s2 * np.linalg.pinv(X.T @ X)[-1, -1])
        dd = max(abs(b[-1] - r.beta) / max(abs(r.beta), 1e-6), abs(se - r.SE) / r.SE, abs(int(x.sum()) - r.n_carriers)); worst = max(worst, dd)
print(f"1. direct regression vs table, 40 genes: worst relative difference (beta, SE) or carrier-count difference = {worst:.2e}"); bad += worst > 1e-4
print(f"   table: N {N}, rows {len(R)} (DEL {(R.dosage == 'DEL').sum()}, DUP {(R.dosage == 'DUP').sum()}); rows with >= 5 carriers {(R.n_carriers >= 5).sum()}, >= 21 {(R.n_carriers >= 21).sum()}")
for dos in ("DEL", "DUP"):
    d = R[(R.dosage == dos) & (R.n_carriers >= 21)]; lam = np.median(d.t ** 2) / 0.4549364
    print(f"   {dos}: genes with >= 21 carriers {len(d)}; lambda (median t^2) {lam:.2f}; P < 2.5e-6: {(d.P < 2.5e-6).sum()}; top: " + ", ".join(f"{s} ({t:+.1f})" for s, t in zip(d.sort_values('P').symbol.head(6), d.sort_values('P').t.head(6))))
if a.meta_one:
    m = pd.read_csv(os.path.join(a.meta_one, a.trait, "gene_assoc_meta.tsv.gz"), sep="\t", dtype={"entrez": str}).merge(R[R.n_carriers >= 5], on=["entrez", "dosage"], suffixes=("_m", ""))
    w = float(np.nanmax(np.abs(m.beta_m - m.beta) / np.maximum(np.abs(m.beta), 1e-6))); print(f"2. one-cohort meta vs table: rows {len(m)}, worst relative beta difference {w:.2e}"); bad += w > 1e-4
if a.meta_two:
    m = pd.read_csv(os.path.join(a.meta_two, a.trait, "gene_assoc_meta.tsv.gz"), sep="\t", dtype={"entrez": str}); two = m[m.k_cohorts == 2]
    rb = (two.beta / two.beta_UKBB); rs = (two.SE / two.SE_UKBB)
    print(f"3. two-copy self-test: rows {len(m)}, in both copies {len(two)} (the min_cell-20 copy drops genes with < 21 carriers); beta ratio {rb.min():.4f}-{rb.max():.4f} (expect 1), SE ratio {rs.min():.4f}-{rs.max():.4f} (expect 0.7071)")
    bad += abs(rs.median() - 0.7071) > 1e-3
print("OK" if not bad else "FAILED"); sys.exit(1 if bad else 0)
