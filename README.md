# cnv-dosage-architecture

Code and summary results for Sacks et al. (2026), "Gene dosage architecture of anthropometric traits": how copy-number
variants (CNVs) change height and body mass index (BMI), split into a mirror component (deletion and duplication of the same
genes act in opposite directions), a constraint component and a single-hit component.
Jonathan Sebat (UCSD), analyses carried out with Claude, July to November 2026.

The companion paper on body mass as a modifier of growth and autistic traits in CNV carriers has its own repository,
https://github.com/jsebat/cnv-growth-bmi. The two repositories were split from one on 2026-11-19; this one keeps the full
history of the shared code (the `cnv-gwas` branch of cnv-growth-bmi). The biobank association tests and meta-analysis are in
Molly Sacks's repository, https://github.com/mollysacks/Sacks2026.

## Layout

| Folder | What it holds | Figures |
|---|---|---|
| `fig1_main_effects/` | dose-response panels of recurrent CNV effects on height and BMI | Fig. 1 |
| `recurrentCNVgenotyping/` | recurrent-CNV genotyping tutorial | Fig. 1 |
| `main_effects/` | federated recurrent-CNV main effects on height, BMI and head circumference, all ages and < 18 y vs >= 18 y (run per cohort, aggregates only), meta-analysis, Fig. 1 style plots | Fig. 1 |
| `fig3_mirror_asymmetric/` | mirror line, perpendicular offset and its bimodality, mirror/continuous/single-hit model comparisons, extended two-component model grid | Figs. 2-3, Supplementary Note 5 |
| `concordance/` | common-variant TWAS concordance of CNV dosage effects (circular genomic permutation), exome-gene concordance | Fig. 2, Supplementary Note 4 |
| `fig4_pathway_analysis/` | pathway convergence of asymmetric CNVs and single-gene height associations | Fig. 4 |
| `growth_models/` | locus-level joint models of the 34 recurrent intervals: bivariate (deletion and duplication) model, mirror decomposition (TWAS share), single-hit presence vs count, family-level leave-one-out | Fig. 5A-D, H-I |
| `cnv_gwas/` | UK Biobank genome-wide CNV analyses: probe-level CNV GWAS, ConGA axes and joint models with genome-shift nulls and leave-one-locus-out, aggregate results | Fig. 5E-G, Supplementary Figs. |
| `AllOfUs/` | the same pipeline packaged for All of Us (WGS and arrays), meta-analysis and figure code (`figures/make_fig_models.py` = Fig. 5) | replication |

`growth_models/code/` scripts were run from the lab Dropbox (`two_axis_CNV_dose_response/v2_2026-11-08`); their absolute
paths point there. Inputs are per-locus summary statistics (`growth_models/inputs/`); gene annotations (RefSeq exons, GeneBayes
s_het, ClinGen, DDG2P) and the public TWAS tables (Mancuso et al. 2017; PsychENCODE) are not redistributed because of their size.

No participant-level data are included. Run `cnv_gwas/scripts/id_scan.sh AllOfUs` (and `cnv_gwas`) before every push.
