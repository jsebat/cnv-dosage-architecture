# BMI TWAS dosage panel: GTEx only, and GTEx + PsychENCODE, with tissue sensitivity

MS 2026-11-03. Per CNV, BMI dosage effect (DUP - DEL)/2 vs summed BMI TWAS z of its genes (each gene's z = mean over panels).
Same non-overlapping CNV set as twas_dosage_bmi.py (Jon); all CNVs mirror for BMI. Circular permutation P (10,000, seed 2026).
Genes keyed by current HGNC symbol (hgnc_map.tsv; previous/alias symbols mapped).

Modes (`python bmi_gtex_analysis.py <mode>`, then `python fig_bmi_gtex.py <mode>`):
- `gtex`     : TWAS30 GTEx v6 panels only (Mancuso 2017; GIANT BMI GWAS). 6 of 22 CNVs have no GTEx-modelled gene -> n = 16.
- `gtex_pec` : GTEx + PsychENCODE cortex (PECfull weights, Tan et al. 2024 BMI GWAS; ../../psychencode_twas/raw/BMI_TAN).
               NB different GWAS for the two sources. PsychENCODE covers 169 of the 224 modelled CNV genes, so it dominates. n = 22.
Outputs per mode: bmi_<mode>_{points,stats,tissue_rho,gene_counts,tissue_sensitivity_table}.csv,
fig_bmi_<mode>_scatter.{png,pdf}, fig_bmi_<mode>_tissue_sensitivity.{png,pdf}.
Sensitivity figure shows only tissues with >= 5 CNV genes with a model (MIN_GENES); Bonferroni over tissues shown; the hidden
tissues still contribute to the all-tissue mean. * one-sided P < 0.05, ** Bonferroni.

Results
- GTEx only (n 16): rho 0.12, P 0.33 (without Williams 0.05, P 0.43).
- GTEx + PsychENCODE (n 22): rho -0.20, P 0.80. PsychENCODE alone -0.11 (P 0.68); GTEx alone in this gene universe 0.03 (P 0.42).
- Single tissues: a few nominal positives (sun-exposed skin, spleen, heart LV), none Bonferroni.
- The sources disagree strongly for some CNVs: 16p11.2 BP1-BP3 summed z +30.9 PsychENCODE vs -10.1 GTEx; 2q11.2 -9.6 vs +8.5.
