# Tissue sensitivity of the TWAS dosage panel (Fig. panel A, height)

MS 2026-11-03. Panel A uses each gene's TWAS z averaged over all TWAS30 panels (45 tissues). Here the panel A statistic (all
non-overlapping CNVs except Williams, n = 22; Spearman rho, one-sided circular permutation P, 10,000, seed 2026) is recomputed
with z from (a) one tissue only and (b) all tissues except one. Baseline reproduces rho 0.69, P 0.0008.

1. `tissue_drivers.py <outdir>` (run with `.` from this folder) -> tissue_drivers_rho.csv (single / leave-one-out / tissue-group
   rho and P for all, mirror, asymmetric), tissue_drivers_cov_decomp.csv (each panel's additive share of cov(summed z, dosage
   effect)), tissue_contrib_per_cnv.csv.
2. `tissue_gene_counts.py` -> tissue_gene_counts.csv (CNV genes with a model in each tissue / in all but that tissue).
3. `fig_tissue_sensitivity.py` -> fig_tissue_sensitivity.{png,pdf}, tissue_sensitivity_table.csv. Bars = rho; * P < 0.05, ** P < 0.05/43
   (Bonferroni over tissues with models in these CNVs); (n) = CNV genes with an expression model.

Result: no single tissue drives it. Blood is the only group significant alone (rho 0.42, P 0.036; GTEx whole blood 0.47, P 0.02)
and carries ~60% of the covariance (YFS ~44%); dropping any one panel leaves P < 0.012. Brain alone rho 0.29 (CMC 0.33, P 0.07).
Blood's lead may reflect the large YFS/NTR eQTL panels (better-powered weights) rather than tissue biology.
