# Sacks 2026 in All of Us: handoff for Milad

Sebat lab, 2026-11-16. Code: branch `cnv-gwas` of https://github.com/jsebat/cnv-growth-bmi, folder `AllOfUs/`.
Questions: Jon.

## In short

**What it is.** The Sacks 2026 drafts describe the CNV dose response of height and BMI. Every genome-wide panel in them comes from UK Biobank (UKBB) today. You will:
1. run the same analyses in All of Us (AoU), for both traits;
2. meta-analyse the two cohorts;
3. redraw those panels from the meta-analysis.

**Same code, same summary statistics.** The code is the code that produced the UKBB results; only the config file changes. Every AoU output has the same form as its UKBB counterpart. Every genome-shift null in AoU reuses the exact UKBB shift list, so the two cohorts can be combined replicate by replicate.

**Three parts:**

| Part | Where | What | Time |
|---|---|---|---|
| A | Expanse | Build and verify the UKBB inputs AoU needs (shift lists, UKBB reference export, resources bundle), then a plumbing test of the AoU config on 5,000 UKBB people | ~½ day, mostly waiting on jobs |
| B | AoU Workbench | Prep, CNV GWAS, burden, additivity, axis scores, univariate grid, joint models, export | 1–2 days of VM time |
| C | Workbench or laptop | Meta-analysis (`meta/meta_sacks2026.py`), then the final figures (`figures/`) | < 1 h |

**What you produce, for height and BMI:**

| Analysis | Step | Figure panels |
|---|---|---|
| Probe-level CNV GWAS (port of Omar Shanta's PGC GWAS), Miami plot, lead loci | `gwas`, `plots` | Supplement; meta Miami |
| Allele tests at the 64 Table S9 CNV alleles | `alleles` | Supplement |
| DEL and DUP effects at 60 growth loci (34 Sacks intervals + 26 UKBB-discovered loci) | `loci` | Supplement, replication |
| Genome-wide burden: signed burden bins, U-shape test, removal of single-hit duplication carriers | `burden` | Fig. 2A-B, 3O-P |
| Constrained-gene dose response; sub-additivity; asymmetric vs mirror duplications; single hit vs constraint count | `additivity` | Fig. 2I-J, 3H-I, model figure |
| Cross-fitted mirror and constraint scores and their additivity across CNVs | `scores` | Fig. 2L |
| Univariate grid of TWAS and constraint measures, with removal of the 34 recurrent intervals | `grid` (`grid_univariate`) | Fig. 2G, 2K, model figure E |
| Mirror by call size (> 5 vs ≤ 5 genes) | `grid` (`grid_callsize`) | Fig. 2H |
| Joint models: all CNVs, 24 single-unit removals, all minor units together, no recurrent loci | `grid` (`joint_sacks34`) | Model figure F-G |
| Joint models: all CNVs and 4 carrier removals (single-gene DUPs, single-gene CNVs, 13 asymmetric loci, both) | `joint` | Joint-model heatmap, UK Biobank columns |

**What we did not port.** Jon decided there is no tissue or cell-type grid of TWAS panels: it is overkill and showed nothing. Each trait uses one TWAS panel set:
- **height:** GTEx non-brain;
- **BMI:** all 45 panels for the joint models (Jon, 2026-11-16) and GTEx non-brain elsewhere, as in the UKBB runs behind the current figures.

There is one extra run: the BMI univariate grid on all 45 panels, in both cohorts. See section 6, item 1.

## 1. Before you start

- **Access.** You need:
  - the AoU Controlled Tier, with genomic data;
  - an Expanse account on `ddp195`;
  - to be listed on the lab's UK Biobank application. Part A reads UKBB phenotypes and CNV calls in `/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/data`. If you are not listed, ask Jon; Part A can then be run by him and you start at Part B with the bundle.
- **Code.** Clone the repository on Expanse and, later, in the Workbench:
  ```bash
  git clone -b cnv-gwas https://github.com/jsebat/cnv-growth-bmi.git
  ```
- **Python.** Python 3 with numpy, scipy, pandas and matplotlib; nothing else. No R and no bedtools. The figures in Part C need R with data.table, ggplot2 and patchwork, on Jon's laptop.
- **Your own output folder on Expanse.** Do not write into `jsebat/`:
  ```bash
  OUT_ROOT=/expanse/projects/sebat1/$USER/sacks2026_aou
  ```

## 2. Expanse rules (Jon's; please follow them)

### Use ind-compute well, and leave room for the lab

From the lab guide `james-guevara/sebat-lab-guides`, `slurm-expanse.md`, section "Limits, and getting past 16 nodes", commit 31b0bb7, and Jon's rules from 2026-11.

**The limits** (QOS, checked 2026-11-05):

| Partition | Binding limit |
|---|---|
| `ind-shared` | The **account** may hold at most 16 nodes. Many small jobs land on many nodes, so one person's array can fill all 16 and everyone else in `ddp195` waits (`MaxNodePerAccount`). |
| `ind-compute` | 16 running jobs **per user**. The account may use 32 nodes and 4,096 cores. |

**Large batches go to ind-compute as 2-node jobs.**
- Twelve 2-node jobs give 24 nodes (3,072 cores) and still leave nodes and job slots for the lab:
  ```bash
  #SBATCH --partition=ind-compute
  #SBATCH --qos=ind-normal
  #SBATCH --nodes=2
  #SBATCH --ntasks=2
  #SBATCH --ntasks-per-node=1
  #SBATCH --cpus-per-task=128
  #SBATCH --mem=240G            # per node
  #SBATCH --nice=100
  srun --ntasks=2 --ntasks-per-node=1 --cpus-per-task=128 bash -c 'worker.sh'
  ```
- The two workers must not do the same work. Each walks the same task list and claims a line atomically (`mkdir claims/<line>` succeeds for exactly one worker).
- Template: `cnv_gwas/jobs/node2_tasks.sbatch` in this repository. To submit a slice: `sbatch --export=ALL,TASKS=<file>,FROM=1,TO=56,SLOTS=32 node2_tasks.sbatch`.
- **Use many 2-node jobs, not one.** For a large project, split the task list into FROM..TO slices and submit one 2-node job per slice, up to about 12, sized so the whole project finishes fast.

**Etiquette.**
- Stay below the account cap (Jon's ConGA run used 24 of 32 ind-compute nodes and at most 6 of 16 ind-shared nodes).
- Always use `--nice`.
- **Never fill ind-shared with arrays of small tasks.**
- Add nothing while another `ddp195` user has a job pending on an account cap in the same partition. Check first:
  ```bash
  squeue -A ddp195 -p <partition> -t PD -o "%u %i %r"
  ```
  Look for the reasons `MaxNodePerAccount` and `AssocGrp*`.
- ind-shared jobs must state `--nodes=1 --ntasks=1`. Without `--ntasks`, sbatch rejects the job silently.

The Part A jobs are small (1 node, 2–32 cores, ≤ 6 h), so they go to ind-shared. Nothing in this project needs ind-compute unless you rerun a UKBB grid (section 6, item 1).

### Job hygiene

- Submit from a **frozen code snapshot**, which every sbatch here copies into its output folder, so editing the code cannot break a running job.
- **Smoke test first**, on a slice that fails in seconds.
- Run `sbatch --test-only` before every submit.
- **One output folder per job variant.** Two configurations must never write to the same path; the sbatch files here refuse an existing folder.
- **scp scripts to the cluster and run them there. Never pipe code through ssh**: nested quotes and `$` get eaten.
- Cancel jobs only by numeric job id (`scancel <id>`), never `scancel -u` or `-n`. Never cancel someone else's running job, or one of Jon's, without asking.
- Report a job as a job id plus the command that checks it. "Submitted" is not "finished".
- Nothing in `/tmp`: it is not shared between login and compute nodes.

### Save space for other users

From `filesize-expanse.md` and Jon's practice in this project.

**The file system.** `/expanse/projects/sebat1` is NFS, about 92% full, with **no per-user quota**. Whatever you write comes out of everyone's space.

**Rules:**
- **Do not copy what you can link.** `assemble_ukbb_reference.py` builds the UKBB tree from symlinks; only the small aggregate export is a real copy.
- **Delete person-level intermediates as soon as the job is done.** The sbatch files here remove the `pheno_*.tsv` carrier subsets and `stratum_pheno.tsv` at the end; `make_union_strata.py` rebuilds them in a minute.
- **Do not keep what you will not use.** Delete smoke-test outputs once the full run is verified. The null arrays (`null_runs.npz`) are kept only where a meta-analysis needs them, and the export stores them as float32.
- **Measure in a job, not on the login node.** Use `diskus -j <cores>` inside an ind-shared job, and `lfs quota` on Lustre. `du` over a large tree on the login node is slow and gets killed.
- **Before planning a cleanup, check hard links** (`nlink > 1`). Deleting one link of a hard-linked file frees nothing.
- The output of this project is small: shift lists 0.2 MB; UKBB reference export about 100–200 MB; plumbing test a few GB, deleted after review.

## 3. Part A: Expanse (UKBB side)

All commands are run from a folder of yours on the Expanse login node. Each job writes its log in the submit folder.

```bash
REPO=/expanse/projects/sebat1/$USER/cnv-growth-bmi/AllOfUs
OUT_ROOT=/expanse/projects/sebat1/$USER/sacks2026_aou
PY=<python3 with numpy/scipy/pandas>
```

### A1. UKBB shift lists (about 1 h)

```bash
sbatch --test-only --export=ALL,REPO=$REPO,OUT_ROOT=$OUT_ROOT,PY=$PY $REPO/ukbb/make_ukbb_shifts.sbatch
sbatch --export=ALL,REPO=$REPO,OUT_ROOT=$OUT_ROOT,PY=$PY $REPO/ukbb/make_ukbb_shifts.sbatch
```

**Why this exists.** The genome-shift null rotates the gene annotation (TWAS z, s_het and ClinGen HI together) along the genome. Each replicate shift is drawn from the seed and kept only if it preserves the CNV burden of the analysed people. So the list depends on:
- the trait;
- the TWAS panel set;
- the carrier subset.

The UKBB runs did not save their lists. The job regenerates the 11 lists those runs used:

| Trait | Panel set | Carrier subsets |
|---|---|---|
| height | GTEx non-brain | all, noUnionDUP, noAsym13, noUnionDUP_noAsym13, noUnionAny |
| BMI | GTEx non-brain | all |
| BMI | all 45 panels | all, noUnionDUP, noAsym13, noUnionDUP_noAsym13, noUnionAny |

**Check.** `shared_shifts_sacks2026_2026-11-16/shared_shifts/SHIFTS_MANIFEST.tsv` must say `yes` in the `match` column on every row that has a reference run. The job compares each regenerated list with the `perms_sha1` recorded by the original UKBB run, so `yes` means identical.
- A `NO` stops the job. Tell Jon; do not ship those files.
- "no sha1 recorded" (older runs) or "ref run not found" is acceptable only if the same list matched `yes` against another run.

### A2. UKBB reference export (about 15 min)

```bash
sbatch --export=ALL,REPO=$REPO,OUT_ROOT=$OUT_ROOT,PY=$PY $REPO/ukbb/assemble_ukbb_reference.sbatch
```

This lays out the existing UKBB results exactly as the AoU pipeline writes its own, then exports them with the same exporter. It recomputes nothing and uses symlinks. The output is `ukbb_reference_sacks2026/ukbb_reference_sacks2026_export.tar.gz` (+ `.sha256`).

Read `ASSEMBLY_MANIFEST.tsv`. Every row should be `ok` once the UKBB grid rerun has finished (section 6, item 1).

### A3. Resources bundle (light; login node is fine)

```bash
OUT_ROOT=$OUT_ROOT bash $REPO/ukbb/make_aou_bundle_sacks2026.sh
```

The output is `aou_bundle/cnv_resources_sacks2026_2026-11-16.tar.gz` (+ `.sha256`), about 15 MB. It holds:
- the public annotations (refGene exons, GeneBayes s_het, ClinGen, Collins 2022 dosage scores, TWAS30, segmental duplications, Omar's allele definitions, growth-locus lists);
- the UKBB probe positions;
- the verified shift lists.

It holds no individual data.

### A4. Plumbing test (about 2–4 h; do not skip it)

```bash
sbatch --export=ALL,REPO=$REPO,OUT_ROOT=$OUT_ROOT,PY=$PY $REPO/ukbb/aou_path_test_sacks2026.sbatch
```

This runs every AoU step on 5,000 UKBB people re-encoded as AoU inputs, through the same code path and config as the Workbench but with small nrep. It then runs the meta-analysis against the UKBB export.

**Check:**
- The log ends with the meta table count and the first lines of `joint/all/height/grid_summary.tsv`.
- `shared_shifts` is `True` in the meta grid tables.
- The run reports the minutes per step. Scale them by the AoU sample size to plan Part B.

When A4 passes:
- copy the two tarballs (bundle and UKBB export) to your computer;
- delete `$OUT_ROOT/aou_path_test_*` (section 2, "Save space").

## 4. Part B: AoU Workbench

### B1. Workspace

- **Workspace:** Controlled Tier, CDR v8 or later.
- **Jupyter environment:** 64–96 vCPUs, ≥ 128 GB RAM, ≥ 200 GB persistent disk. Each process peaks at about 2 GB, and the grid runs many variants in parallel.

In a notebook terminal:
```bash
cd ~ && git clone -b cnv-gwas https://github.com/jsebat/cnv-growth-bmi.git
gsutil cp $WORKSPACE_BUCKET/cnv_resources_sacks2026_2026-11-16.tar.gz ~/cnv-growth-bmi/AllOfUs/
cd ~/cnv-growth-bmi/AllOfUs && sha256sum cnv_resources_sacks2026_2026-11-16.tar.gz && tar xzf cnv_resources_sacks2026_2026-11-16.tar.gz
```

Upload the two tarballs to the workspace bucket first (Workbench file upload). They contain no individual data. The bundle creates `AllOfUs/resources/`, which is where `configs/aou_sacks2026.json` expects it.

**Fill two paths** in `cnv_anthro/configs/aou_sacks2026.json`, both marked `FILL_FROM_DATA_DICTIONARY`. Take them from the Workbench page "How the All of Us genomic data are organized" for your CDR:
- the srWGS **GATK-SV joint-called VCF**;
- **`ancestry_preds.tsv`**, copied locally first.

Use `gsutil -u $GOOGLE_PROJECT` for the requester-pays buckets. The notebook `cnv_anthro/workbench/cnv_anthro_aou.ipynb` walks through the same steps.

### B2. Run

All commands run in `~/cnv-growth-bmi/AllOfUs/cnv_anthro`. Run long steps in a terminal with `nohup ... &`; `workbench/dsub_run.sh` sends a step to a batch VM.

```bash
C=configs/aou_sacks2026.json
python run_pipeline.py --config $C --steps prep                                  # SV VCF -> CNV table; OMOP -> phenotypes (< 1 h)
python run_pipeline.py --config $C --steps gwas,plots,alleles,loci,burden        # every stratum (all, 6 ancestries, female, male)
python run_pipeline.py --config $C --steps additivity,scores --strata all        # minutes
nohup python run_pipeline.py --config $C --steps grid,joint --strata all --jobs 24 > ~/cnv_anthro_run/grid_joint.log 2>&1 &
python run_pipeline.py --config $C --steps export
```

**Phenotypes and CNVs** (unchanged from the CNV GWAS handoff):
- **Height and BMI:** median of all adult measurements (age ≥ 18).
  - Units are converted, and implausible values are dropped (height 120–220 cm, BMI 12–75).
  - Each trait is rank-inverse-normal transformed within sex, again inside every stratum.
  - Covariates: age and age² at the measurement, sex, and 16 ancestry PCs.
- **CNVs:** GATK-SV DEL and DUP sites, plus multi-allelic CNV sites (CN < 2 = DEL, CN > 2 = DUP).
  - Filters: PASS (or MULTIALLELIC), ≥ 10 kb, autosomes, and the PGC 10% frequency rule.
  - Only people with SV calls enter the phenotype table.

**What the grid and joint steps run** (stratum `all`):

| Output | Variants | Nulls | Removal runs | UKBB shift list |
|---|---|---|---|---|
| `grid/grid_univariate/<trait>/` | 44 per trait (adjusted for total genes) | 1,000 | 24 single-unit removals + no_loci | `<trait>_gtex_nonbrain_all` |
| `grid/grid_callsize/<trait>/` | 12 per trait | 1,000 | same | same |
| `grid/joint_sacks34/<trait>/` | `height_selected`, `bmi_selected` | 2,000 | 24 single-unit removals + minor_all (8 fixed major units kept) + no_loci | same |
| `grid/grid_univariate_bmi_allpanels/bmi/` | 44 | 1,000 | 24 + no_loci | `bmi_all_all` |
| `joint/<subset>/<trait>/` | 5 carrier subsets × 2 traits | 2,000 | 24 + minor_all + no_loci | `<trait>_<panels>_<subset>` |

`make_union_strata.py` builds the carrier subsets from AoU's own calls:
- **noUnionDUP:** no duplication of a "union" gene (top-5% s_het or ClinGen haploinsufficient);
- **noUnionAny:** no CNV of a union gene;
- **noAsym13:** no CNV matching one of the 13 asymmetric intervals at 50% reciprocal overlap;
- **noUnionDUP_noAsym13:** both of the first and third.

A missing shift file stops the run (`require_shared_shifts`). That is deliberate: AoU's own shifts cannot be meta-analysed replicate by replicate.

**Time.** The grid and joint steps are the main cost. In UKBB (N = 345k), one joint run with 2,000 nulls took 45–75 min on 8 threads. Use the minutes from the plumbing test, scaled to the AoU N, for your own estimate.

### B3. QC before export

1. **Sample sizes.** Check `results/<stratum>/taskA_gwas/<trait>/run_info.json`, `burden/<trait>/run_info.json` and `joint/strata/union_strata_info.json` for N, carriers, union-gene carriers and asymmetric-interval carriers.
2. **Allele frequencies.** At the S9 loci, compare frequencies with UKBB using `compare_allele_freq.py`. In European ancestry, recurrent CNVs (16p11.2, 22q11.21, 15q11.2 BP1-BP2) should have similar carrier frequencies. A tenfold difference is a call-set or overlap problem, not biology.
3. **The "only" runs** (`additivity/<trait>/onlyasym_*`, `onlymirror_*`). These keep carriers of the listed intervals plus **people with no CNV call at all**.
   - On arrays (UKBB) that is a large reference group. With WGS, most people carry at least one call ≥ 10 kb, so the reference group may be small.
   - Check N in their `run_info.json` and tell Jon the number before relying on Fig. 3H-I (section 6, item 2).
4. **Small binary terms.** In `additivity_models.tsv`, terms such as `hit_DUP` are presence terms. Check that `hit_carriers_DEL` and `hit_carriers_DUP` in `run_info.json` are > 20; the export removes counts from 1 to 20 but not a beta that rests on them. If one is ≤ 20, tell Jon and leave that run out of the download.
5. **Genomic inflation.** The Miami plots show probe-level inflation. Large inflation in one stratum usually means a covariate or ancestry problem.

### B4. Export and egress (read before downloading anything)

AoU rules (Dissemination Policy, egress alert and large-download pages; checked 2026-11-15):

**No participant count from 1 to 20 may leave the Workbench.**
- This includes anything from which such a count can be worked out.
- 0 is allowed.
- The FAQ also mentions allele counts below 40. Our tables report carrier counts (participants). If the in-Workbench article "Downloading and Disseminating Genomic Data" asks for ≥ 40, set `"min_cell": 39` in the config, rerun `--steps export`, and tell Jon.

**No participant-level data and no screenshots of it.**

**Egress alerts.**
- Large downloads trigger an alert: the VM is stopped and the security team sends a questionnaire.
- Ask support@researchallofus.org for the current size limit before downloading.
- Above the limit, file a Large Download Exemption Request.

**Publications.** Before any publication or talk that uses AoU results, file the Publication and Presentation Reporting Form at least 2 weeks ahead.

**What the export step does** (`export_aggregate.py`):
- It copies aggregate files only. Phenotype tables, CNV tables and anything with `pheno`, `genotype` or `.cnv` in the name are never copied.
- **Primary rule:** every row with a count from 1 to 20 loses its count (written as `<=20`) and its statistics (written as `suppressed`). This covers burden bins, dose-response bins, probes, alleles and loci.
- **Secondary rule:** a small count must not be recoverable by subtracting strata (e.g. male = all − female).
- It writes `EXPORT_MANIFEST.txt`.

**Download procedure:**
1. Run the meta-analysis in the Workbench (Part C1).
2. `tar czf aou_sacks2026_export_<date>.tar.gz --exclude='*.npz' export_aggregate meta`. The null arrays stay in the Workbench once the meta-analysis has used them.
3. Read `EXPORT_MANIFEST.txt`, spot-check that no count column shows a value from 1 to 20, and download once via the Jupyter file browser.
4. Put it in Dropbox `Jon/Claude/UKBB_anthro_CNV/aou_results/` and tell Jon.

## 5. Part C: meta-analysis and final figures

### C1. Meta-analysis (in the Workbench, so only small tables leave)

Upload `ukbb_reference_sacks2026_export.tar.gz` to the bucket, then:
```bash
mkdir -p ~/ukbb_export && tar xzf ukbb_reference_sacks2026_export.tar.gz -C ~/ukbb_export     # gives ~/ukbb_export/all
cd ~/cnv-growth-bmi/AllOfUs/cnv_anthro
python ../meta/meta_sacks2026.py --ukbb ~/ukbb_export/all --aou ~/cnv_anthro_run/export_aggregate/all --out ~/cnv_anthro_run/meta
```

How each analysis is combined:

| Analysis | Method |
|---|---|
| CNV GWAS, alleles, loci | Inverse-variance fixed and random effects, Cochran Q and I², on the UKBB probe grid; meta Miami plot |
| Bins and models | Inverse-variance per bin or term, carriers summed |
| Sub-additivity (Fig. 3I) | Meta of each class, then the asymmetric − mirror difference |
| Score additivity (Fig. 2L) | Meta of each cohort's ratio of slopes, carriers of 2+ CNVs vs one; the scores are cohort-specific, the ratio is not |
| Grid and joint models | √N-weighted Stouffer of t, computed for every null replicate under the shared shifts; then empirical z, P, and how many removals hold. `shared_shifts` must be `True` on every row |

The output uses the same folder and file names as the UKBB aggregates (`dup_burden_2026-11-15/`, `growthloci_2026-11-15/`, …), plus `joint/<subset>/<trait>/grid_summary.tsv` and `cnv_gwas/<trait>_all/`.

### C2. Final figures (Jon's laptop)

The locus panels read the lab Dropbox. Download the meta folder to `Jon/Claude/UKBB_anthro_CNV/aou_results/meta`, then, in `AllOfUs/figures`:
```bash
GENOMEWIDE_ROOT=<meta> GENOMEWIDE_LABEL="UK Biobank + All of Us" FIG_OUT=<out> python make_fig2_fig3_reorg.py
GENOMEWIDE_ROOT=<meta> GENOMEWIDE_LABEL="UK Biobank + All of Us" FIG_OUT=<out> python make_fig_models.py
GENOMEWIDE_META=<meta> FIG_OUT=<out> Rscript fig_joint_model_heatmap.R
```
Without these variables, the same scripts draw the current UKBB-only figures. We checked that both modes run and that the meta mode reproduces the UKBB values when both cohorts are given the same tables.

## 6. Open points (Jon decides; flag them when you hit them)

1. **UKBB grid rerun (Jon, 2026-11-16).**
   - The UKBB univariate grid, call-size and joint_sacks34 runs predate the code that saves `null_runs.npz`, so they could not be meta-analysed.
   - Jon approved a rerun with the same inputs and the verified shift lists: `results/growthloci_v2_2026-11-16`. It includes the BMI univariate grid on all 45 panels.
   - `assemble_ukbb_reference.py` reads that folder. Wait for it before running A2.
2. **Reference group of the "only" runs on WGS** (B3, item 3). If it is small, the definition may need to change in both cohorts, e.g. to "no gene-hitting call".
3. **Platform.** UKBB calls come from arrays; AoU calls come from WGS (GATK-SV). WGS finds many more small calls, so burden bins, the "single CNV" groups and gene counts per person are not identical in meaning. The meta tables keep each cohort's estimate and the heterogeneity test; look at I² before interpreting a pooled value.
4. **Shift matching in AoU.** AoU reuses the UKBB shifts, which were burden-matched to the UKBB CNVs. The null stays a valid rotation in AoU, but the burden match is approximate there. This is the price of replicate-matched meta-analysis.

## 7. Design notes for the methods

- **CNV → gene.** A deletion hits a gene if it overlaps any exon. A duplication hits a gene only if it spans a whole transcript.
- **Gene annotations:**
  - Constrained genes: GeneBayes s_het > 0.029 (imputed values excluded).
  - Union single-hit genes: top-5% s_het (≥ 0.2272) or ClinGen haploinsufficiency 2–3.
  - TWAS: FUSION TWAS30. Per-gene z is either Σz/√(number of panels) ("scaled") or the mean over panels ("mean").
- **Joint models.**
  - Height: total genes, summed |TWAS z| of constrained genes, and the constrained-gene count on each arm; union presence on each arm (all-CNV model).
  - BMI: total genes and the count of constrained genes with |mean TWAS z| > 4 on each arm; union presence on each arm; the constrained-gene count on duplications.
  - Readouts: each arm, and the mirror, (t_DUP − t_DEL)/√2.
- **Null.** A circular shift of the gene annotation, with TWAS z, s_het and HI rotated together. The shift lists are those of UKBB.
- **Removals:**
  - 24 recurrent units (the 34 Sacks intervals with overlapping intervals merged), each with its own matched null;
  - minor_all: the 16 minor units removed together, with the 8 major units kept;
  - no_loci: all units removed.

  "Held" means |z| > 1.96 with the same sign.
- **Meta-analysis.** As in C1.
