# cnv_anthro — rare-CNV analyses of quantitative traits (UK Biobank, All of Us)

One code base, one JSON config per cohort. Two analyses:

- **Task A:** Omar Shanta's PGC probe-level CNV GWAS, adapted to a linear model. Also runs the allele-level tests at the Table S9 loci.
- **Task B:** person-level burden of the two ConGA gene-dosage axes, mirror and constraint (Engchuan et al., companion manuscript). It includes a genome-shift null, leave-one-locus-out over the 18 recurrent loci, and a run with all 18 loci removed.

```
python run_pipeline.py --config configs/<cohort>.json --steps prep,gwas,plots,alleles,axes,export \
    [--only_trait height] [--only_panels gtex_nonbrain]
```

| file | role |
|---|---|
| `prep_ukbb.py` | UKBB adapter. Takes the PGC Freeze_2024 UKBB calls, links IDs, extracts UKB fields, applies QC, and writes the standard phenotype table. |
| `aou_sv_to_cnv.py` | AoU adapter 1. Converts the GATK-SV srWGS VCF to a per-person CNV table in PLINK .cnv layout. |
| `aou_pheno_extract.py` | AoU adapter 2. Builds the standard phenotype table from OMOP (BigQuery) plus the ancestry PCs. |
| `cnv_gwas_qt.py` | Task A probe GWAS. Uses Omar's probes and filters, his independent-test count and his GWS threshold, with an exact sparse OLS in place of the logistic model. |
| `alleles.py` | Port of Omar's allele calling, followed by the trait association. |
| `conga_axes_qt.py` | Task B. ConGA consensus2 axes, joint OLS, burden-matched genome-shift null, 18 LOO runs plus a no-loci run, all with matched nulls. |
| `plots.py` | Miami plot. |
| `export_aggregate.py` | Copies aggregates only, enforcing the small-cell rule (no count from 1 to 20). |
| `configs/` | `ukbb_expanse.json`, `aou_workbench.json`, and `ukbb_as_aou_test.json` (the AoU code path run on UKBB data). |
| `workbench/` | AoU Jupyter notebook and dsub launcher. |
| `tests/` | Synthetic AoU unit test, plus the builder for the UKBB-as-AoU test. |

Requirements: Python 3 with numpy, scipy, pandas and matplotlib. On AoU, add google-cloud-bigquery. All of these are preinstalled in the Workbench image. Nothing else is needed: no bedtools and no R.

## Standard tables (the adapter contract)

- **Phenotype table:** one row per analysed person, with columns `IID sex(1/2) age age2 [batch] PC1..PCk height bmi height_rint bmi_rint`. The traits are rank-inverse-normal transformed within sex (Blom).
- **CNV table:** PLINK .cnv layout with columns `FID IID CHR BP1 BP2 TYPE SCORE SITES`. TYPE 1 = DEL, 3 = DUP. CHR is numeric.
- **People without sequencing or calls:** they must not appear in the phenotype table, because anyone in it with no calls is treated as a non-carrier.

## Task A choices

These follow Omar's current code: `scripts/ALL/<DIS>/noFam_noSex` and `META_noSex`, from Nov 2025.

- **Probes:** the PGC master breakpoint files, restricted to probes overlapped by at least one cohort call of that type. With `bp_dir: null`, the cohort's own breakpoints are used instead.
- **Model:** OLS of trait_rint on carrier status plus covariates (age, age², sex, genotyping batch/array, PCs). DEL and DUP are tested separately, on chr1–23 (UKBB) or chr1–22 (AoU).
- **Probe filters:**
  - carrier frequency > 2% (known loci are exempt);
  - fewer than 12 carriers (AoU: fewer than 21, the small-cell rule);
  - more than 20 QC-filtered calls at the probe (PGC filter beds, UKBB only);
  - overlap with a segmental duplication;
  - for DUP, the blacklist region chr11:134824089-134923878.
- **Threshold:** 0.05 / n_independent, where independent tests come from Omar's greedy Jaccard merge (≥ 0.04) of the call sets at adjacent probes.

## Task B choices

These follow ConGA `conga_joint.R --design consensus2`.

- **Mirror:** count_A = n(z > 3) − n(z < −3) over the genes a person's CNVs hit.
  - The TWAS is the GTEx-based FUSION TWAS30 for height or BMI. Per-gene z = Σz / √(number of panels).
  - Primary: GTEx non-brain panels. Sensitivity: brain panels (GTEx brain + CMC), and all panels (the panel set of the earlier locus-level analysis).
- **Constraint:** none_H = n(genes with GeneBayes s_het > 0.029). Imputed s_het values are excluded.
- **CNV → gene:** a DEL hits a gene if it overlaps any exon by ≥ 1 bp. A DUP hits a gene only if it spans a whole transcript (any isoform).
- **Model:** trait ~ covariates + [n_all, count_A, none_H] for DEL + the same three terms for DUP.
  - Readouts: DEL arm, DUP arm, and mirror net = (t_DUP − t_DEL)/√2 for count_A; DEL and DUP arms for none_H.
  - Secondary variants, quantified alongside the primary: unsigned n(|z| > 3), and split n(z > 3) with n(z < −3).
- **Null:** ConGA true_perm, i.e. one circular shift of the autosomal gene circle per replicate.
  - Each shift is redrawn until the median burden ratio is within 10%. There are 2,000 replicates.
  - The same shift is used in all 20 runs: full, 18 LOO and no_loci.
- **LOO:** `loo2_runs.tsv` (18 runs). Calls are assigned to loci with the `assign_loci.awk` rule.
  - "Held" means the empirical |z| > 1.96 with the same sign as the full run.
  - A min-over-LOO statistic is also tested against the same minimum taken in each null replicate.

## All of Us

- **CNV source:** the srWGS GATK-SV callset, covering about 97k participants (CDRv8/v9). AoU releases no curated array CNV callset.
  - DEL/DUP sites: a person is a carrier if their genotype carries the non-reference allele.
  - Multiallelic CNV sites: CN < 2 is a DEL, CN > 2 a DUP.
  - Filters: PASS, ≥ 10 kb, site frequency ≤ 10%, autosomes only.
  - The array rule "≥ 10 probes" has no WGS analogue and is not applied.
  - Joint-called sites give all carriers the same breakpoints, so each site yields two probes with identical carrier sets. The Jaccard merge collapses them into one independent test.
- **Before running:**
  - Fill the two marked paths in `configs/aou_workbench.json` (the SV VCF and `ancestry_preds.tsv`) from the Workbench Data Dictionary.
  - Run prep in the notebook VM, which has `gsutil`.
  - Run gwas/alleles/axes in the notebook or with `workbench/dsub_run.sh`.
- **What leaves the Workbench:** only `export_aggregate/`.

## Speed: fast code vs reference

`reference_check.py` keeps the naive implementation runnable as a reference: a full dense OLS refit per test, which is what a per-probe `lm`/`glm` loop does. The fast code computes the same OLS Wald t-test exactly; no statistic is changed. The check below is the UKBB chr22 height run (N = 345,169, 119 covariate columns, 4 CPUs, job 55739871, 2026-11-14).

**Equivalence.** Probe GWAS, 400 random chr22 probes:
- carrier counts are identical;
- the largest relative difference is 4.4e-6 for β and SE, and the largest absolute difference in −log10 P is 1.8e-6.

These differences are entirely the 6-significant-digit rounding of the saved sumstats. For the axes model, t differs by at most 2.7e-12 between the fast and reference fits, on the observed annotation and on 3 shifted ones.

**Timing.**

| step | reference (naive refit per test) | fast | how |
|---|---|---|---|
| one probe test | 3.0 s | ~2 µs per probe (sweep) | Covariate QR reused (FWL); carrier sums updated incrementally along the chromosome (only calls entering or leaving) |
| chr22 GWAS (9.4k probes) | ~7.9 h (extrapolated) | 18 s wall, 2.2 GB | as above |
| whole-genome GWAS + plots + 154 allele tests, one trait | ~17 CPU-days (extrapolated, ~500k probes) | 1 min 44 s | as above |
| one axes fit (6 dosage terms) | 3.49 s | 0.017 s | Gene-space sufficient statistics (Q′G precomputed per run; F built from carriers only) |
| axes: 2,000 nulls × 3 variants × 20 runs, one trait × panel set | ~5 days (extrapolated) | 35-40 min on 4 CPUs | as above |

`--workers N` (fork-parallel over the 20 runs) produces bitwise-identical output to `--workers 1`; checked on height and BMI, 2,000 nulls each (job 55740332). It was **slower** as run: 57 min with 20 workers on a 128-core node, against 35 min with 1 worker on 4 CPUs. Each worker started a full 128-thread BLAS pool, which oversubscribes the node. Set `OPENBLAS_NUM_THREADS = cores / N` before using it. The default stays at 1 worker, and no speedup is claimed for it.

## Validation

All checks were run on Expanse on 2026-11-14.

- **Probe-level re-implementation:** reproduces Omar's META_N_noSex and META_N probe P values exactly at the 22q11.21, 16p11.2, 3q29 and 15q11.2 lead probes for SCZ.
- **Allele level:** reproduces Table S9 exactly at 3q29 DEL, 22q11.21 A-D DEL, 15q11.2 BP1-BP2 DEL and 16p11.2 BP4-BP5 DUP.
- **Allele-call port:** all 154 alleles are identical to Omar's UKBB (XD Axiom) genotype matrix, with 0 discordant calls among 346,202 people.
- **GWAS sweep:** on chr22, the probe set and the carrier count per probe are identical to Omar's XD Axiom sumstats (3,434 DEL and 5,993 DUP probes).
