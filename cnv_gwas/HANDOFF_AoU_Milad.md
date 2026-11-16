# CNV GWAS of height and BMI in All of Us: handoff for Milad

> **Superseded 2026-11-16 by [`AllOfUs/HANDOFF_Milad.md`](../AllOfUs/HANDOFF_Milad.md)**, which covers this CNV GWAS plus every other genome-wide analysis of the Sacks 2026 drafts, the meta-analysis and the figures. Use that one.

Sebat lab, 2026-11-15. Code: branch `cnv-gwas` of https://github.com/jsebat/cnv-growth-bmi (folder `cnv_gwas/`).

## What this is, in one paragraph

We ran a rare-CNV analysis of height and BMI in UK Biobank (PGC CNV Freeze_2024 array calls, N = 345k) and want the same analysis in All of Us (srWGS GATK-SV calls, ~97k genomes) for a two-cohort meta-analysis. You run one pipeline in the AoU Researcher Workbench. It produces summary statistics overall, by genetic ancestry and by sex. You then export **aggregate results only**, following the AoU small-cell rule. The pipeline is already tested end to end on 5,000 UKBB people re-encoded as AoU inputs, with the same code, configs and steps.

**What you run:**
1. **Probe-level CNV GWAS**: deletions and duplications separately, linear model. It is a port of Omar Shanta's PGC CNV GWAS. It tests AoU's own breakpoints plus the UKBB probe positions (for the meta-analysis).
2. **Allele-level tests** at the 64 Table S9 CNV alleles (Shanta et al.).
3. **Per-locus DEL and DUP effects** at 60 growth loci: 34 Sacks recurrent intervals and 26 loci discovered in the UKBB CNV GWAS.
4. **Genome-wide gene-dosage grid**:
   - person-level burden of TWAS-gene counts (signed and unsigned, several thresholds) and of constrained-gene counts, on deletions and duplications;
   - a call-size split (>5 / >10 genes);
   - two joint models;
   - each with a genome-shift null that uses **the same shifts as UKBB**, so the two cohorts can be meta-analysed replicate by replicate.
5. **Export**: aggregate tables only, the small-cell rule, cross-stratum suppression and a manifest.
6. **Meta-analysis**: optional, but it is the best way to keep the download small. Upload the UKBB aggregate export into the workspace and run the two meta scripts there; only small meta tables then need to leave.

## 1. Set up the workspace

- **Workspace:** a Controlled Tier workspace with genomic data access, on the current CDR (v8 or later).
- **Cloud environment (Jupyter, Python):** we suggest 64–96 vCPUs, ≥ 128 GB RAM and a ≥ 200 GB persistent disk. Each process peaks at about 2 GB. The grid runs many variants in parallel.
- **Code** (in a notebook terminal):
  ```
  cd ~ && git clone -b cnv-gwas https://github.com/jsebat/cnv-growth-bmi.git
  cd ~/cnv-growth-bmi/cnv_gwas/cnv_anthro
  ```
  Requirements: Python 3 with numpy, scipy, pandas, matplotlib and google-cloud-bigquery. All are preinstalled in the Workbench image. No R and no bedtools.
- **Resources:** one 11 MB bundle of public annotations and UKBB aggregates (no individual data).
  - Expanse: `/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/aou_bundle/cnv_gwas_resources_2026-11-15.tar.gz`, checksum in `.sha256` next to it.
  - Copy it to your computer, upload it to the workspace bucket (Workbench file upload, or `gsutil cp` into `$WORKSPACE_BUCKET`), then:
    ```
    cd ~/cnv-growth-bmi/cnv_gwas && gsutil cp $WORKSPACE_BUCKET/cnv_gwas_resources_2026-11-15.tar.gz . && tar xzf cnv_gwas_resources_2026-11-15.tar.gz
    ```
    This creates `cnv_gwas/resources/`, which is where `configs/aou_workbench.json` expects it (`../resources`).
- **Two paths to fill** in `configs/aou_workbench.json`, both marked `FILL_FROM_DATA_DICTIONARY`. Take them from the Workbench's genomic data organisation page for your CDR version.
  - The joint-called **GATK-SV VCF** (srWGS structural variants).
  - **`ancestry_preds.tsv`** (research_id, ancestry_pred, pca_features). Copy it locally first.
  - Use `gsutil -u $GOOGLE_PROJECT` for the requester-pays genomic buckets.

## 2. Run

Everything goes through one runner. Steps can be run separately. Strata are set in the config:
- `all`;
- `eur`, `afr`, `amr`, `eas`, `sas`, `mid` (AoU `ancestry_pred`);
- `female`, `male`.

```
cd ~/cnv-growth-bmi/cnv_gwas/cnv_anthro
python run_pipeline.py --config configs/aou_workbench.json --steps prep                  # SV VCF -> CNV table; OMOP -> phenotypes
python run_pipeline.py --config configs/aou_workbench.json --steps gwas,plots,alleles,loci  # all strata
python run_pipeline.py --config configs/aou_workbench.json --steps grid --jobs 24           # grid strata: all, eur, afr, amr, female, male
python run_pipeline.py --config configs/aou_workbench.json --steps export
```
- **Long runs:** use `nohup ... &` in a terminal, or `workbench/dsub_run.sh` to send steps to a batch VM.
- **Splitting the run:** `--strata eur,afr` or `--only_trait height` restrict a run, so it can be split across sessions.

**Phenotypes and CNVs** (handled by `prep`):
- **Height and BMI:** median of all adult measurements (age ≥ 18).
  - Units are converted, and implausible values are dropped (height 120–220 cm, BMI 12–75).
  - Each trait is rank-inverse-normal transformed within sex, again within every stratum.
  - Covariates: age and age² at the median measurement, sex, and 16 ancestry PCs.
- **CNVs:** GATK-SV DEL and DUP sites plus multi-allelic CNV sites (CN < 2 = DEL, CN > 2 = DUP).
  - Filters: PASS (or MULTIALLELIC for CNV sites), ≥ 10 kb, autosomes, and the PGC 10% frequency rule.
  - Only people with SV calls enter the phenotype table, so everyone analysed has been sequenced.

**Expected time:**
- `prep`: under 1 h (the VCF stream is the slow part).
- `gwas`, `plots`, `alleles`, `loci`: a few minutes per stratum and trait.
- `grid`: the main cost. It is 102 variants per trait in each of 6 strata, each with 1,000–2,000 null replicates over 11 matched runs. Expect several hours on 64–96 vCPUs. Run `all` first, then the ancestry and sex strata.

## 3. Check before export (QC)

- `results/<stratum>/taskA_gwas/<trait>/run_info.json`: N, calls and carriers per stratum, plus probes tested and kept.
  - With AoU's per-genome calls (no 10-probe rule on WGS), expect many more small calls than on arrays.
- **Allele frequencies at the S9 loci:** compare with UKBB using `compare_allele_freq.py`. Recurrent CNVs (e.g. 16p11.2, 22q11.21, 15q11.2 BP1–BP2) should have similar carrier frequencies in European ancestry. A 10× difference means a call-set or overlap problem, not biology.
- **Probe-level genomic inflation:** the Miami plots show it. Large inflation in one stratum usually means a covariate or ancestry problem.

## 4. Export and egress (read this before downloading anything)

AoU rules, from the AoU support pages (Dissemination Policy, egress alert, large downloads; checked 2026-11-15):
- **No participant count from 1 to 20 may leave the Workbench.** This includes values from which such a count can be worked out (percentages, differences between tables, other formulas). 0 is allowed. Show small cells as "≤20" or suppress them.
- **No participant-level data** may be downloaded, and no screenshots of it.
- **Egress alerts:** downloading files from Jupyter, large outputs, or large printouts can trigger an alert.
  - When one fires, the VM is shut down and the security team sends a questionnaire.
  - High-volume egress can disable the account until it is reinstated.
  - The size threshold is not published. **Ask support@researchallofus.org for the current limit before downloading.**
  - Files above the limit that comply with the policy need a **Large Download Exemption Request**: a 48-hour window, approved per workspace and user.
- **Genomic data:** read the in-Workbench article "Downloading and Disseminating Genomic Data". It is not public.
  - The AoU FAQ says **allele counts below 40** must not leave either.
  - Our tables report carrier counts (participants), and the export withholds counts 1–20.
  - If that article requires carrier counts ≥ 40 for genomic tables, rerun the export with `--min_cell 39` (set `"min_cell": 39` in the config). Tell Jon if so.
- **Before any publication or presentation** that uses AoU results: file the Publication and Presentation Reporting Form at least 2 weeks ahead. Include the AoU acknowledgement and a data-availability statement (Controlled Tier, CDR version).

**What the export step does** (`export_aggregate.py`, run by `--steps export`):
- **Copies only aggregate files.** Phenotype tables, CNV tables and anything with genotype, pheno or .cnv in the name are never copied.
- **Primary rule:** every per-probe, per-allele and per-locus row whose carrier count is 1–20 loses its count ("<=20") and its effect statistics ("suppressed"). The SE of a carrier effect would otherwise reveal the count. Withheld probe rows are dropped from the probe tables.
- **Secondary rule:** across strata, a small count must not be recoverable by subtraction, e.g. male = all − female, or the participants without an ancestry call = all − Σ ancestries.
  - For every row key, if any child stratum or that remainder has 1–20 carriers, the row is suppressed in every child stratum of that partition.
- **Writes `EXPORT_MANIFEST.txt`**: total size, file list and the number of withheld rows.

**Download procedure:**
1. `python run_pipeline.py --config configs/aou_workbench.json --steps export`
2. Read `export_aggregate/EXPORT_MANIFEST.txt`. Check the total size against the limit support gave you.
3. Spot-check the tables: `zcat .../sumstats_del.tsv.gz | head`. No count column should show a value from 1 to 20.
4. Optional, recommended (meta-analysis in the Workbench):
   - Upload our UKBB aggregate export (on Expanse, `results/export_ukbb_strata/`, available once the UKBB stratified run finishes) into the workspace.
   - Then:
     ```
     python meta_sumstats.py --cohort UKBB=<ukbb_export>/all --cohort AoU=export_aggregate/all --trait height --out meta/height_all --ref_label UKBB
     python meta_grid.py --cohort UKBB=<ukbb_export>/all/grid/grid_univariate/height --cohort AoU=export_aggregate/all/grid/grid_univariate/height --out meta/grid_univariate_height_all.tsv
     ```
     Repeat for each trait and stratum, and for the `grid_callsize` and `joint_2026-11-15` designs.
   - `meta_grid.py` checks that both cohorts used identical shifts (`shared_shifts` column).
   - The meta tables are small and contain no counts below 21.
5. Package what you will download: `tar czf aou_cnv_export_<date>.tar.gz export_aggregate meta`. Download it once, via the Jupyter file browser.
   - The grid null-replicate arrays (`null_runs.npz`, about 0.2 MB per variant and stratum) are the bulk of the export.
   - If you run the meta in the Workbench, leave them out of the download (`--exclude='*.npz'`).
6. Put the download in `Jon/Claude/UKBB_anthro_CNV/aou_results/` (Dropbox) and tell Jon.

## 5. What a reviewer will want (all produced per stratum: all, each ancestry, female, male; each trait)

| table | file (under `export_aggregate/<stratum>/`) | content |
|---|---|---|
| probe GWAS, full | `taskA_gwas/<trait>/sumstats_{del,dup}.tsv.gz` | CHR BP NCNVs BETA SE T P N FREQ, filters, independent-test group, REF_PROBE (UKBB position) |
| lead loci | `taskA_gwas/<trait>/lead_{del,dup}.tsv`, `run_info.json` | groups at P < 1e-4, genome-wide flag, threshold (0.05 / independent tests), N, call counts |
| Miami plot | `taskA_gwas/<trait>/miami_<trait>.png` | |
| S9 alleles | `taskA_alleles/<trait>/allele_assoc.tsv` | 64 alleles × DEL/DUP: carriers, frequency, beta (SD), SE, P |
| growth loci | `locus_effects/locus_effects_<trait>.tsv` | 60 loci × DEL/DUP: carriers, beta, SE, P, rule, median call size |
| grid | `grid/<design>/<trait>/<variant>/grid_summary.tsv` | observed t, empirical z/P, held across the 8 major loci, all-minor-loci and all-loci removals |
| grid nulls | `.../null_runs.npz`, `observed_terms.tsv`, `per_run_z.tsv` | for the replicate-matched meta-analysis |
| QC | `run_info.json` files, `EXPORT_MANIFEST.txt` | sample sizes, call counts, suppression counts |

## 6. Design notes (for the methods section)

- **Probes:** AoU breakpoints, unlike the PGC array breakpoints, so probe positions differ from UKBB. `--ref_probes` also tests the UKBB positions, which makes a probe-level meta possible (`meta_sumstats.py --ref_label UKBB`).
- **Filters:**
  - The array rule "≥ 10 probes" has no WGS analogue and is not applied.
  - Minimum carriers per probe is 21, so the small-cell rule is built into the tests.
- **Grid null:**
  - Genome-shift of gene annotations (TWAS z and s_het rotated together), using the UKBB shift list (`resources/shared_shifts/<trait>_shifts.npy`).
  - Leave-one-out uses the 34 Sacks recurrent intervals (24 families). The 8 major families, chosen in UKBB as those whose removal moves any readout by |Δt| ≥ 1, each get their own matched run. The remaining families are dropped together in one run, plus an all-loci run.
- **Meta-analysis:**
  - Probes, alleles and loci: inverse-variance fixed and random effects, with Cochran's Q and I².
  - Grid: Stouffer of t weighted by √N, computed replicate by replicate under the shared shifts.

Questions: Jon. The UKBB side (code, results, handoff) is in `HANDOFF_UKBB_anthro_CNV.md` next to this file.
