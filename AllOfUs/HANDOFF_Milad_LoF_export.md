# All of Us: full probe P values and gene-level LoF (DEL + SNV) export (handoff for Milad, 2026-11-25)

## Why

We want to show where exome/genome loss-of-function (LoF) SNV signals and CNV deletion signals converge on the same genes. The current
exports drop every probe and gene with 1–20 carriers, so genes with few deletions are invisible in All of Us. One example is CPZ on 4p16.1:
in UK Biobank, 25 deletion carriers have height −0.84 SD (P 1e-5), but there is no All of Us row at all.

Jon (2026-11-25): "Don't exclude gene based association p-values of any kind. but we should be careful with raw counts data."

**New export rule.** Every probe and gene row is exported with its **t statistic and P value**. When a row has 1–20 carriers:
- the count, beta, SE and frequency are withheld;
- the count is written as `<=20` (probes) or as NA with `n_carriers_masked = 1` (genes);
- beta and SE are withheld because together they give the count back (SE ≈ 1/√n); t and P alone do not.

Counts that could be recovered by subtraction stay withheld as before: across strata, and now also across the gene classes DEL / PLOF / LOF.

## What you run

Pin the new commit of `jsebat/cnv-dosage-architecture` (branch `main`; Jon will send the hash). All three steps reuse results you
already have. **No CNV GWAS or other fit is rerun.**

### 1. Probe level: re-export the CNV GWAS (WGS and arrays; minutes)

Your results folders already hold every probe (`results/<stratum>/taskA_gwas/...`). Only the export changes.
- Copy each config.
- Set `"export_dir"` to a new folder (an export is never overwritten), e.g. `${HOME}/cnv_anthro_run/export_aggregate_allP_20261125`.
- Point `"results_dir"` at the folder that holds your CNV GWAS.

```
python run_pipeline.py --config <new config> --steps export
```

- Withheld probe rows now stay in `sumstats_*.tsv.gz`, with `T` and `P` filled, `BETA`/`SE`/`FREQ` = `suppressed`, `NCNVs` = `<=20` and
  `stats_withheld = 1`.
- The files get larger: roughly every probe with at least 1 carrier.
- The log line still ends with the suppression counts. Package the folder as before.

### 2. Gene level, deletions and duplications: re-export (WGS and arrays; minutes)

The in-Workbench gene tables already contain every gene with at least 1 carrier. In the same new config:
- set `"genes": {"min_carriers": 1, ...}` (already the default in the repo configs);
- set `"genes": {"export_dir": "<new folder>"}`.

```
python run_pipeline.py --config <new config> --steps genes_export
```

- Check that the log ends with `AUDIT OK`.
- New columns:
  - `n_carriers_masked`: 0 = shown; 1 = 1–20 carriers; 2 = more than 20 carriers but withheld because subtraction would give a small count;
  - `freq_max`: an upper bound on frequency, 21/N, for masked = 1.

### 3. Gene level, LoF SNVs (WGS dataset; about an hour, mostly Hail)

**a. Extract rare pLoF carriers per gene from the short-read WGS SNV/indel callset (Hail cluster).**

```
python cnv_anthro/plof_carriers_aou.py --vat <VAT path of the current CDR> --mt <srWGS exome split Hail MatrixTable> \
    --people <the pheno table of the WGS config> --out ${HOME}/plof/plof_carriers_wgs.tsv.gz
```

- **Variants:** canonical-transcript `stop_gained`, `frameshift_variant`, `splice_acceptor_variant`, `splice_donor_variant` with VAT
  `gvs_all_af` < 1%. These are the classes of the Kosmicki 2026 and Akbari 2021 exome pLoF masks, without LOFTEE, which the VAT does not
  carry.
- **Carriers:** het or hom-alt genotype, `FT` PASS.
- **Output:** `IID`, `symbol`. It is person-level and stays in the Workbench.
- **Before running:**
  - Check the VAT column names (`vid`, `gene_symbol`, `consequence`, `is_canonical_transcript`, `gvs_all_af`) against your CDR. The script
    stops if one is missing.
  - Fill in the VAT and MatrixTable paths from the genomics data page.
- **This script has not been run.** It may need small fixes for your CDR (column names, `FT` type, the `GT` vs `LGT` field). Send Jon the
  error if it fails.
- It prints the share of analysed people who have srWGS genotypes. **Use the output only where that share is 100%.** Otherwise people
  without srWGS are counted as non-carriers. This should hold for the WGS dataset. For the array dataset, run step 3 only if it also prints
  100%.

**b. Re-run the gene step with the pLoF table, then export.**

In the WGS config:
- set `"genes": {"plof": "${HOME}/plof/plof_carriers_wgs.tsv.gz", ...}`;
- set `"results_dir"` to the folder that holds your gene tables.

The gene step rewrites `gene_assoc/<trait>/gene_assoc.tsv.gz` in place. DEL and DUP rows come out identical; `PLOF` and `LOF` rows are
added.

```
python run_pipeline.py --config <WGS config> --steps genes
python run_pipeline.py --config <WGS config> --steps genes_export      # into another new genes export_dir
```

- `PLOF` = carries a rare pLoF SNV in the gene.
- `LOF` = carries a pLoF SNV **or** an exonic deletion (DEL = any exon overlapped by at least 1 bp, the same rule as before).
- The export withholds a count whenever two shown counts would give a 1–20 count by subtraction: LOF − DEL, LOF − PLOF, or DEL + PLOF − LOF.
  The audit checks it.

## What to send

1. The three export folders (probe; gene DEL/DUP for WGS and arrays; gene with PLOF/LOF for WGS), packaged and copied to Expanse as before.
2. The printed lines of `plof_carriers_aou.py`: carrier × gene pairs, genes, and the srWGS share.

## Please confirm before egress

- **Policy check.** Is releasing t and P for rows with 1–20 carriers, with no count, beta, SE or frequency, within the AoU Data and
  Statistics Dissemination Policy? Jon believes it is, because no count can be recovered from them. If you or the AoU team disagree, run
  step 1 with `--drop_small_probe_rows`, set `min_carriers` to 21 for the genes, and tell Jon.
- **Changes to the earlier export rule.** These replace the 2026-11-23 rule: the 21-carrier gene floor, and dropping withheld probe rows.
  Nothing else in the earlier handoffs changes.
