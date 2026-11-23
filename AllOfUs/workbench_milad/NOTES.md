# All of Us Workbench run (Milad, Nov 2026)

`make_nb.py` generates `sacks2026_aou_workbench.ipynb` (`python make_nb.py sacks2026_aou_workbench.ipynb`); the notebook is committed without outputs. Molly's recurrent CNV genotyping port lives in `../../recurrentCNVgenotyping/aou/` and is embedded into notebook section 15 at build time.

## Pipeline code

No changes to `AllOfUs/cnv_anthro/*.py`. The runs used upstream commits only (sex strata fix 4b64136; `dup_burden` small strata and "only"-run reference b97150f / ce1ce67; grid/joint 05fa65b); on this repo that is `main` at 946fb1f. Everything Workbench-specific is in the notebook:

- **SV → CNV per chromosome** (stage `sv`). CDR v9 ships one GATK-SV VCF per chromosome; the 22 are converted in parallel with `aou_sv_to_cnv.py` and the config's filters, then concatenated. The PGC 10% rule is computed within chromosome over the same samples, so this equals a single-file run.
- **Staged background driver**: stage markers, a bucket backup after each stage, and a restore on a new app. A paused VM kills background jobs; the run moved from a Spark (Dataproc) app to a plain JupyterLab app midway, so paths are relative to `$HOME`.
- **Array inputs from Molly's PGC calls** (section 12): `results/4_filtered/GDA-8v1-0/*.plt`, consensus calls only (`algorithm` = IP,UP: PennCNV + QuantiSNP; DELDUP dropped); breakpoints = `consStart/consEnd`; sample QC `Freq_Calc == 1`; genotyping batch from `sample_info`; SITES = probes in `probes.txt` inside the call. No call merging (fragmented calls are kept, as in the PGC freeze). Set `MOLLY_GS` to the bucket of the PGC CNV workspace.
- **Array config (Jon 2026-11-20)**: "only"-run reference `no_rare_gene` (+ `nocall` sensitivity, `_nocall`); burden `single_def = rare_gene`. A second burden version with `all_calls` runs as its own stage (`burden_allcalls` → `<stratum>/burden_all_calls/`), because `run_pipeline.py` writes one burden version per run.
- **JupyterLab app**: `google-cloud-bigquery` was newer than `google-api-core` (import error); section 3 upgrades both when the import fails.
- **QC and export**: 8.6 (1q21.1 TAR DUP breakpoints, WGS), 8.7 (calls split within a locus, array), and an export audit (9.1) on top of `export_aggregate.py`'s masking.

## WGS vs array

| | WGS (srWGS, GATK-SV v9) | Array (PGC calls, GDA-8v1-0) |
|---|---|---|
| People | srWGS participants | participants **without** srWGS (no overlap) |
| Call filters | PASS/MULTIALLELIC, ≥ 10 kb, autosomes, PGC 10% rule | IP,UP consensus, ≥ 10 kb, ≥ 10 probes, PGC 10% rule |
| Covariates | age, sex, 16 PCs | + genotyping batch (factor) |
| Calls per person (≥ 10 kb) | median 57; ~55% of midpoints in segdups | median 8; 0.1% with no call (UKBB: 1, ~21%) |
| "Only"-run reference | `no_rare_gene` (+ `_no_rare_gene_nosd`) | `no_rare_gene` (+ `_nocall`: ~460 people) |
| Burden panel | `rare_gene` | `rare_gene` + `all_calls` |
| Extra run | B6: no GWAS segdup filter | — |

Known call-set features:

- **WGS**: 1q21.1 TAR BP1-BP2 DUP excess from two segdup-anchored read-depth sites; TAR BP2-BP3 not callable; calls at 1q21.1 distal are larger than the locus, so they fail the reciprocal rule.
- **Array**: calls are fragmented at segdup-flanked loci. DUP carrier frequencies vs UKBB: 22q11.21 A-D about 0.4×, 1q21.1 about 0.3×, 17q12 and TAR BP1-BP3 lower still; TAR BP1-BP2 DUP is high (pieces of BP1-BP3 duplications). Deletions and 16p11.2 match UKBB.

## Recurrent CNV genotyping (`recurrentCNVgenotyping/aou/`)

A port of `GenotypeRecurrentCNVs_tutorial.ipynb`: same locus sheet (saved as TSV; no Google access in the Workbench), same 50% reciprocal / one-way rules, same segment patterns. It uses pandas instead of bedtools, and its input is the pipeline's PLINK-CNV table, so one script serves WGS and array calls.

`crosscheck_molly.py` runs the notebook cell by cell with its bedtools binary on 4,000 simulated people. The results are identical except:

- (a) people with both a DEL and a DUP at one locus, where "first call" depends on bedtools' output order;
- (b) the notebook reads each bedtools output with `pd.read_csv`'s default `header=0`, which drops the first hit of each intersect. The port keeps it.

On AoU, the port and `alleles.py` agree at all the main loci in both datasets. They differ only at:

- loci `alleles.py` doesn't define (16p12.2 BP1-BP2, 22q11.21 A-C / E-H / D-H, 15q11.2-13.1 BP2-BP3, 16p11.2 BP3-BP4);
- WGS 1q21.1 distal (the reciprocal rule);
- array 22q11.21 sub-segments and 15q13.3 BP4-BP5, where exact-pattern genotyping gives mixed fragment patterns no allele.
