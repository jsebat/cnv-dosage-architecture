# AllOfUs: UK Biobank + All of Us replication and meta-analysis (Sacks 2026)

Everything needed to repeat, in All of Us, every genome-wide UK Biobank analysis of the Sacks 2026 drafts, for height and BMI:
- the CNV GWAS;
- the S9 alleles;
- the growth-locus effects;
- the genome-wide burden;
- the univariate grid and the joint models of the gene-dosage axes;
- additivity.

The folder also holds the code that meta-analyses the two cohorts and redraws the figures from the meta-analysis.

**Start here:** [`HANDOFF_Milad.md`](HANDOFF_Milad.md).

**Fig. 3 joint model of the dosage paper (added 2026-11-23):** [`HANDOFF_Milad_Fig3.md`](HANDOFF_Milad_Fig3.md).
- `cnv_anthro/run_fig3.py`: fits, in-cohort summaries and the meta-ready export, from one spec for every cohort
  (`designs/fig3_spec_2026-11-23.json`, `designs/fig3_joint_2026-11-23.json`); also reachable as `run_pipeline.py --steps fig3,fig3_export`.
- `ukbb/assemble_ukbb_fig3.py`, `ukbb/fig3_ukbb_export_meta_2026-11-23.sbatch`: the UK Biobank export from the finished runs.
- `meta/meta_fig3.py`, `ukbb/meta_fig3_expanse.sbatch`: the meta-analysis at UCSD; `meta/compare_fig3_tables.py` checks two table sets.
- `figures/figs_dosage_story.py D` with `GW2_FIG3_META=<meta folder>` redraws Fig. 3 from the meta-analysis.
- Rare-CNV filter (2026-11-23): `cnv_anthro/cnv_freq_filter.py`, set in the spec (`cnv_filter`: per-position carrier frequency, calls mostly
  above 1% dropped) and applied by `run_fig3.py` before every fit, in every cohort.
- Univariate tests and dose-response bins of every Fig. 3 stratum: spec fits `uni_*`, `cnv_anthro/fig3_dose_bins.py`; figure key `B`.
- Gene-based association (2026-11-23): `cnv_anthro/gene_assoc_qt.py` (`run_pipeline.py --steps genes,genes_export`),
  `cnv_anthro/export_gene_assoc.py`, `meta/meta_gene_assoc.py`.
- Milad's Workbench notebook and notes: `workbench_milad/`; his All of Us port of the recurrent CNV genotyping: `../recurrentCNVgenotyping/aou/`.

| folder | what it holds | where it runs |
|---|---|---|
| `cnv_anthro/` | **The pipeline**, the same code for both cohorts:<br>• single entry point `run_pipeline.py`, with one config per cohort<br>• `configs/aou_sacks2026.json`: the All of Us run<br>• `configs/ukbb_as_aou_sacks2026_test.json`: the plumbing test on Expanse<br>• `designs/`: the term sets of the univariate grid and the joint models<br>• `workbench/`: the notebook and the dsub launcher | AoU Workbench, Expanse |
| `ukbb/` | **UK Biobank side:**<br>• regenerate and verify the UKBB genome-shift lists that AoU must reuse<br>• lay out the existing UKBB results in the AoU export format<br>• build the resources bundle<br>• run the plumbing test | Expanse |
| `meta/` | `meta_sacks2026.py`, which meta-analyses every analysis and writes figure-ready tables with the UKBB file names | Workbench or laptop |
| `figures/` | The current figure scripts, with the data root as a parameter:<br>• Figs 2 and 3<br>• the model figure<br>• the joint-model heatmap | Laptop (the locus panels read the lab Dropbox) |
| `resources_small/` | The interval lists (34 Sacks intervals, 13 asymmetric and 21 mirror intervals, LOO units) | |

**Data rules:**
- No individual-level data are in this repository, and none may be added.
- UK Biobank individual data stay on Expanse.
- All of Us data stay in the Workbench.
- Only `export_aggregate/` (where no count from 1 to 20 appears) and the meta tables leave the Workbench.
