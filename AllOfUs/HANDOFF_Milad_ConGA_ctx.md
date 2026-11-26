# All of Us: person-level ConGA of tissue, cell-type and signal-group gene sets (handoff for Milad, 2026-11-25)

## Why

The CNV GWAS paper (UK Biobank + All of Us) has a ConGA figure and Supplementary Fig. 1 that are UK Biobank only:

- CNVs of highly constrained (HC) genes in growth signal groups lower height;
- CNVs of HC genes in neural groups raise BMI;
- deletions and duplications act in the same direction.

ConGA needs individual-level CNV calls, so it has to run inside the Workbench. The test is the same code UK Biobank used. Only test
statistics leave the Workbench. At UCSD we then combine UK Biobank and All of Us replicate by replicate, under the same genome shifts.

## What you need

1. **Code.** Pin the commit Jon sends (branch `main` of `jsebat/cnv-dosage-architecture`). New and changed files:
   - `cnv_anthro/conga_ctx_qt.py` (new): the ConGA engine.
   - `run_pipeline.py`: new step `ctx`.
   - `export_aggregate.py`: exports the step's tables without carrier counts.
   - `configs/aou_sacks2026*.json`: new `ctx` block.
2. **Add-on resources bundle.** `cnv_resources_ctx_2026-11-25.tar.gz` (+ `.sha256`), 8.2 MB, in
   `/expanse/projects/sebat1/jsebat/sacks2026_aou/aou_bundle/`. Upload it to the workspace bucket and unpack it in `AllOfUs/` next to the
   2026-11-16 bundle. It adds:
   - `resources/conga_ctx/all_sets_membership.tsv.gz`: 1,650 gene sets, which are the context sets (GTEx tissues, HPA cell types,
     embryonic limb, brain) and the height/BMI signal groups, each split into all / HC / LC genes;
   - `resources/conga_ctx/all_sets_families.tsv`, `loo_units_meta12.tsv` (the 12 locus-removal units), `loo_set_list.txt` (the sets
     that get removal runs) and `loci_sacks34.tsv`;
   - `resources/shared_shifts/height_shifts.npy` and `bmi_shifts.npy`: the 2,000 UK Biobank genome shifts. They are **not** in the
     2026-11-16 bundle. `SHIFTS_MANIFEST_ctx.tsv` confirms they are the files the UK Biobank runs used.

   The exon file comes from the 2026-11-16 bundle (`resources/conga/hg38_refGene_20200708.exon.txt`). Nothing in either bundle is
   individual data.

## What you run

Use the existing results folders and configs. Only the new step runs; nothing else is refit.

```
cd ~/cnv-growth-bmi/AllOfUs/cnv_anthro            # your clone at the pinned commit
python run_pipeline.py --config configs/aou_sacks2026.json --steps ctx --strata all
python run_pipeline.py --config configs/aou_sacks2026_array.json --steps ctx --strata all
```

- **What it runs.** For height and BMI, deletions and duplications, the step runs:
  - the full test of every set against the 2,000 shared shifts;
  - the 14 locus-removal runs (12 loci, the remaining recurrent intervals, all 34 intervals) for the sets in `loo_set_list.txt`.

  The output goes to `<results_dir>/all/ctx/{full,loo}_<trait>_<DEL|DUP>/`. A folder that already has `info.json` is skipped, so an
  interrupted run can simply be restarted.
- **Runtime.** UK Biobank (345,169 people, 16 cores, 120 GB) took about 30 minutes per trait × CNV type for the full run, and about
  33 minutes for the removal runs. With `"jobs": 4` the eight runs take about 1–1.5 hours per dataset on a 16-core, 128-GB machine.
  All of Us has fewer people, so expect less. Memory scales with the number of carriers; use at least 64 GB.
- **Check.** Each run writes `info.json` with `shifts_sha1`:
  - height: `157a247259f02719`;
  - BMI: `be761a32dd575298`.

  These must match. A missing shift file stops the run, deliberately: All of Us's own shifts could not be combined with UK Biobank
  replicate by replicate.

## Export

Export into a new folder; an export is never overwritten.

```
# copy the config, set "export_dir" to a new folder, e.g. ${HOME}/cnv_anthro_run/export_ctx_20261125
python run_pipeline.py --config <copy> --steps export --strata all
```

What leaves the Workbench from each `ctx/<run>/` folder:

- `conga_ctx_<full|loo>.tsv.gz`. Per set and run: z, z_emp, P_emp, null mean and SD, number of genes in the set, FDR. The per-set
  carrier count (`n_carriers`) is removed by the exporter.
  - That count could otherwise be recovered by subtraction between the full run and a removal run.
  - z and P are test statistics and may be shown (Jon 2026-11-25).
- `null_z_<run>.npy`: z of every null replicate (2,000 × sets, float32). These are test statistics only.
- `info.json`: N, the shift checksum, and the number of calls removed per run (values 1–20 are written as `<=20`).

Expected size: about 15 MB per trait × CNV type for the full runs, about 12 MB for the removal runs, about 100 MB in all per dataset.
Check the size against the egress limit before downloading.

## What to send back

- Package `export_ctx_<date>/all/ctx` for WGS and for arrays as `aou_wgs_ctx_<date>.tar.gz` and `aou_array_ctx_<date>.tar.gz`, each with
  a `.sha256`.
- Copy them to `/expanse/projects/sebat1/jsebat/sacks2026_aou/aou_results/`.
- Send the two `info.json` checksum lines (height and BMI) from any one run.

## At UCSD (Jon or the paper's agent)

```
python AllOfUs/cnv_anthro/meta_ctx.py \
   --cohort UKBB=<conga_ctx_56044342>,<conga_supp_56063298>,<conga_fig16_56072429> \
   --cohort AoU_WGS=<unpacked>/ctx --out <new meta dir>
python cnv_gwas/meta_aou_2026-11-24/conga_contrasts.py --stage B --ctx_run <meta dir> --loo_run <meta dir> \
   --groups <conga_supp_56063298>/groups --out <figure dir> --cohort_label "UK Biobank + All of Us"
```

- **How the cohorts are combined.** `meta_ctx.py` combines the cohorts as √N-weighted Stouffer of the set z. It does this for the
  observed z and for every null replicate of every run, which gives a matched null.
- **Readouts and figure.** `conga_contrasts.py` then computes:
  - the four readouts (DEL, DUP, mirror, dose-aligned) in LC and HC genes;
  - FDR and locus-robust marks;
  - the figure and tables, in the same form as the UK Biobank Supplementary Fig. 1.

## Tests (2026-11-25, Expanse job 56073458: all checks PASS)

- `AllOfUs/ukbb/test_ctx_port_2026-11-25.sbatch` runs two checks:
  1. **UK Biobank reproduction.** The ported engine reproduces z, z_emp and P_emp of the UK Biobank runs, for the full run and the
     removal runs.
  2. **Synthetic cohorts.** Two synthetic cohorts (`tests/make_synthetic_ctx.py`) are run through `run_pipeline.py --steps ctx,export`,
     then the export checks, the meta checks (A+A gives √2·z and the same z_emp; A+B uses shared shifts), and the figure drawn from the
     meta output.
- Checks: `tests/test_ctx_port.py`. Results:
  - the UK Biobank full runs (height DEL, BMI DUP; 6 sets each) and the 14 removal runs (height DEL) reproduce z, z_emp and P_emp exactly (max difference 0);
  - the export has no carrier count and float32 nulls;
  - meta A+A gives z = √2·z_A with z_emp unchanged;
  - meta A+B used shared shifts;
  - `conga_contrasts.py` drew the figure from the meta output.
  - Output: `/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/results/ctx_port_test_56073458`.
