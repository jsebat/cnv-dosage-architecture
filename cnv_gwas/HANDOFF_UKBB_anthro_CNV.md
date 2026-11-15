# UKBB height/BMI rare-CNV analyses: handoff

Sebat lab, 2026-11-15. Companion to `report/UKBB_anthro_CNV_report.pdf` (2026-11-14). Started from the "Conga CNV paper V2" session.

## What it is, in one paragraph

Omar Shanta's PGC probe-level CNV GWAS, re-implemented for quantitative traits (linear model) and run on height and BMI in the UK Biobank calls from the PGC CNV Freeze_2024 (Task A). Allele-level tests at the 64 Table S9 alleles. ConGA gene-dosage axes (mirror, constraint; Engchuan et al., companion manuscript) as person-level burden in UKBB, with a genome-shift null and leave-one-locus-out (Task B). One code base runs UKBB on Expanse and All of Us in the Workbench through a config swap.

Headline: the number of highly constrained (HC) genes hit lowers height and raises BMI on both dosages (robust to all 18 LOO and to removing all 18 loci). The ConGA signed mirror is null for both traits; the height z>0-gene mirror is a secondary, quantified alongside (P 0.027).

## Where everything is

**Data rule.** Individual-level UKBB data stay on Expanse in `data/` (group-only, mode 750). Everything under `Jon/Claude/UKBB_anthro_CNV/` (Dropbox) is aggregate.

| What | Expanse (`/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/`) | Local (`Jon/Claude/UKBB_anthro_CNV/`) |
|---|---|---|
| Report | – | `report/UKBB_anthro_CNV_report.pdf` (+ `.html`); builder `report/make_report.py`, text in `report/narrative.json`, inputs `report/export_view/` |
| Code (group-readable) | `code/cnv_anthro/` (runner `run_pipeline.py`, configs, `README.md`) | `code/cnv_anthro/` (master copy; push with `jobs/sync_code.sh`) |
| Code resources | `code/resources/` (`omar_alleles/`, `conga/`, `twas/`, `segDups.bed`) | – |
| Phenotype table | `data/ukbb_pheno_std.tsv` (346,202 people after QC) | – |
| CNV calls | `data/ukbb_pgc_xd.cnv` (PLINK .cnv, PGC XD dataset `ukbb`) | – |
| UKB fields used | `data/ukb42720_anthro_fields.tsv` (from the j3guevar archive `ukb42720.tab`) | – |
| Full results | `results/pipeline_ukbb/` (taskA_gwas, taskA_alleles, taskB_axes; full sumstats) | – |
| Aggregate export | `results/export_ukbb/` | `results_aggregate/pipeline_ukbb/` |
| Replication of Omar | `~/ukbb_anthro/replication/out/` | `replication/*.R`, `results_aggregate/replication/` |
| Fast vs reference check | `results/reference_check/` | `results_aggregate/reference_check/` |
| AoU code path test (UKBB re-encoded) | `results/aou_path_test*`, `data/aou_path_test/` | – |
| Job scripts + frozen snapshots | `jobs/*.sbatch`, `jobs/snap_*` | `jobs/` |
| Logs | `logs/` | – |
| New-loci locus effects (2026-11-15) | `results/locus_effects_newloci_2026-11-15/` | `../two_axis_CNV_dose_response/growth_newloci_2026-11-15/` |
| Shared notes | – | side-session notes `ukbb_anthro_cnv_findings.md` |

**External inputs (read-only, not ours):**
- PGC XD .phe/.cnv: `/expanse/projects/sebat1/oshanta/home/IndividualDisorder/Freeze_2024/inputs/data_10perc/ALL/XD/all_platforms_10kb_10probes_XD_ALL.{phe,cnv}`
- ID link: `/expanse/projects/sebat1/msacks/UKBB/2024/genomatrix_UKBB_2025.csv`
- Omar's current code (reproduced): `/home/oshanta/IndividualDisorder/Freeze_2024/scripts/ALL/<DIS>/noFam_noSex` + `META_noSex`; allele level `scripts_figures/scripts_trumpet_noSex`
- PGC breakpoints: `.../Freeze_2024/inputs/bp_files`; filter beds: `/expanse/projects/sebat1/oshanta/home/CrossDisorder/Freeze_2024/inputs/data_10perc/filteredCNVs/`
- TWAS30 (copied into `code/resources/twas/`): `Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/concordance/twas30/`
- Table S9 workbook: `G:/My Drive/Manuscripts/Shanta_PGC/Supplementary_Tables_2026_S5updated_2026-11-14.xlsx`

## How to rerun

```
cd /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/code/cnv_anthro
python run_pipeline.py --config configs/ukbb_expanse.json --steps prep,gwas,plots,alleles,axes,export
```
Submit through a frozen snapshot (`jobs/ukbb_full.sbatch`; many tasks: `jobs/node2_tasks.sbatch`, 2-node ind-compute). Wall time: GWAS + plots + alleles 1 min 44 s per trait; axes 35-40 min per trait x panel set (2,000 nulls, 20 runs). Report: `python report/make_report.py report/export_view results_aggregate/replication report/UKBB_anthro_CNV_report.html`, then headless Chrome to PDF.

All of Us: `configs/aou_workbench.json` (fill the SV VCF and `ancestry_preds.tsv` paths from the Data Dictionary), `workbench/cnv_anthro_aou.ipynb`, `workbench/dsub_run.sh`. Only `export_aggregate/` leaves the Workbench. Then `compare_allele_freq.py` for array vs WGS allele frequencies.

## Jobs (all finished)

55739525, 55739749, 55739752, 55739815, 55739842, 55739843 (axes array, ind-shared, before the 2-node rule), 55739865, 55739867, 55739871 (reference check), 55740332 (workers identity check), 55742064 (new-loci locus effects, 25 s).

## Caveats to carry forward

1. **Lead-probe and selected-locus effects are inflated (winner's curse).** At the 17 Sacks locus x dosage x trait results that reach genome-wide significance in UKBB, the UKBB lead-probe effect is a median 1.43x the Sacks T3 effect; the UKBB effect over the Sacks interval (50% reciprocal overlap, same carriers almost exactly) is the same size (e.g. 16p11.2 BP4-BP5 DEL height -1.12 lead, -1.11 interval, -0.94 Sacks). Across all 34 Sacks intervals, unselected, UKBB and Sacks agree (slope 0.91-1.20, r 0.68-0.91). So the excess at GWS loci comes from selecting on UKBB significance, not from the lead probe vs interval definition. Effects of the new UKBB-discovered loci carry the same inflation (`growth_newloci_2026-11-15/plus_ukbb/leadprobe_vs_sacks.csv`).
2. **Table S3** in the 2026-11-14 workbook is Omar's Nov-2024 `META` run; it does not reproduce from current files (unknown-sex removal, Oct 2025). Table S9 reproduces exactly.
3. **No current UKB withdrawal list** was found on Expanse; apply the latest before publication.
4. **TWAS30 holds 45 panels.** Primary = GTEx non-brain (32); brain (GTEx brain + CMC) and all panels are sensitivity analyses.
5. **PGC psychiatric cases are included** (29,445 UKBB people with a case flag); **relatedness is not pruned**.
6. **Not on GitHub.** Jon to choose the repo.

## Open items

- **Genome-wide joint model (Jon 2026-11-15).** The genome-wide burden analysis needs one joint model carrying all axes together, instead of separate fits per variant (primary / unsigned / split). To specify with Jon.
- **New height/BMI loci in the growth locus-level analysis (in progress).** 26 UKBB genome-wide-significant loci outside the Sacks set, appended to the Sacks T3 effects (Jon's choice); union list for both traits. Work in `two_axis_CNV_dose_response/growth_newloci_2026-11-15/` (`ref34/` = reproduction of v3 on the 34 intervals; `plus_ukbb/` = 34 + new loci).
- AoU run in the Workbench; array vs WGS allele frequency comparison at the S9 loci.
