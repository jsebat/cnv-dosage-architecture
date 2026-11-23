# Sacks 2026 in All of Us: Fig. 3 joint model (handoff for Milad, 2026-11-23)

This adds one analysis to the run you already finished: the joint model behind the new Fig. 3 of the dosage-architecture paper.
It reuses your CNV and phenotype tables, your resources bundle and your notebook's `launch()` / `status()`. Nothing you ran before is
rerun or changed. The earlier handoff (`HANDOFF_Milad.md`) still holds for everything else.

## In short

1. Pin commit `ed0c509` of `jsebat/cnv-dosage-architecture` (branch `main`) in your notebook and rerun sections 1, 2, 4 and 6.
2. Add the two stages in section A below and launch `fig3` on the WGS dataset (8 joint-model fits, 6 univariate fits, 2 gene-mirror runs
   and 2 binned dose-response runs; several hours).
3. Run the QC cell (section B) and paste its output to Jon. It prints only z statistics, P values and sample sizes.
4. Launch `fig3_export`, check that the log ends with `AUDIT OK`, package the folder and copy the tarball to Expanse (section C).
5. Repeat 2 to 4 on the array dataset when its main run is done (section D).

Do not change the specification (`cnv_anthro/designs/fig3_spec_2026-11-23.json`). UK Biobank ran exactly this file, and the meta-analysis
stops if the two cohorts differ.

## What the model is

One linear model per trait (height, BMI) and per set of CNVs. Each term is a count of genes hit by a person's deletions (DEL) or
duplications (DUP).

| fit | CNVs | terms, for DEL and for DUP | pairs fitted |
|---|---|---|---|
| `size2t` | single-gene CNVs and CNVs hitting 2 or more genes, as separate term sets | total genes hit; genes with TWAS z > +3 ("up") and z < −3 ("down"), each split into less-constrained (LC, s_het ≤ 0.029) and highly constrained (HC); other HC genes (weak or no TWAS). Other LC genes are the reference. | height CNV × height TWAS, BMI × BMI, and the two cross-trait pairs |
| `all2t` | all CNVs pooled | the same terms | height × height, BMI × BMI |
| `totdose` | all CNVs pooled | total genes hit only | height × height, BMI × BMI |
| gene mirror | all CNVs, and single-gene CNVs | per-gene score z for DEL and DUP; correlation of the two across genes | height, BMI |
| `uni_N1_inf`, `uni_N1_1`, `uni_N2_inf` | all CNVs; single-gene; ≥ 2-gene | **univariate**: one model per class (TWAS up, down, weak or no TWAS, each LC and HC; and total genes), with only that class's DEL and DUP gene counts | height, BMI |
| dose-response bins | all CNVs; single-gene; ≥ 2-gene | **univariate**: per class, the trait in fixed bins of the number of genes deleted or duplicated | height, BMI |

TWAS is GTEx non-brain (32 panels), scaled z = sum of z / √(number of panels). Covariates are the ones in your config for each trait.

**Univariate tests and dose-response plots** (Jon 2026-11-23). These come before the joint model: each gene class is tested on its own
(covariates only, no other class terms) in each CNV-size stratum. Two outputs per class and stratum:

- the univariate DEL z and DUP z of the class's gene count, against the same genome-shift null and the same locus-removal runs as the joint
  model (total genes: analytic t);
- the numbers behind a dose-response plot like Fig. 2a–h of the paper: the adjusted difference in the trait (SD) from people with no gene
  of that class, per bin of genes deleted or duplicated, with and without the 34 recurrent intervals. The bins are fixed in the code (`1`,
  `2`, `3-4`, `5-9`, `10+` for the TWAS and HC classes; `1-2`, `3-5`, `6-10`, `11-20`, `21+` for weak-or-no-TWAS LC genes and for total
  genes), so they are the same in every cohort and are meta-analysed bin by bin. As in Fig. 2, people who carry gene-hitting CNVs of both
  dosages are left out. Bin points and their 95% CIs are analytic (OLS) and descriptive.

**Null model.** The joint-model terms are tested against a genome shift of the gene annotation: 2,000 shifts, and they must be the
UK Biobank shifts. You already have them: `resources/shared_shifts/height_gtex_nonbrain_all.npy` and `bmi_gtex_nonbrain_all.npy` from
the 2026-11-16 bundle are byte-identical to the lists the UK Biobank Fig. 3 runs loaded (sha1 `157a247259f02719` and
`be761a32dd575298`). `run_fig3.py` checks the sha1 before and after every fit and stops if it differs. Each fit runs on all CNVs, with
each of 8 major loci removed (16p11.2 BP4–BP5, 16p11.2 BP1–BP3, 22q11.21 A–D, 15q11.2 BP1–BP2, 1q21.1, 16p13.11 BP1–BP2, 17q12, 17p12),
with all minor loci removed (`minor_all`) and with all 34 recurrent intervals removed (`no_loci`): 11 runs, each with its own matched null.
The gene mirror uses a gene-circle shift (10,000 shifts, fixed seed, so the offsets are the same in every cohort). Total-genes terms do not
depend on the annotation and are tested by their analytic t.

**Locus carriers.** The locus removals use the same assignment as your grid and joint runs (`assign_loci` in `conga_axes_qt.py`,
`resources/growth_loci/loci_sacks34.tsv`: a call is assigned to a locus of its own dosage at ≥ 50% reciprocal overlap). This is what UK
Biobank used for Fig. 3, so the two cohorts match. If you genotyped recurrent CNVs another way for any step (for example with
`recurrentCNVgenotyping/GenotypeRecurrentCNVs_tutorial.ipynb`), tell Jon: nothing here reads those genotypes.

## A. Run (WGS dataset)

In section 1 set `PIN_COMMIT = 'ed0c509'`, then rerun sections 1, 2, 4 and 6. Section 6 rebuilds the filled config from the template,
which now has a `"fig3"` block; nothing else in the template changed. Check what changed in the code first:

```
cd ~/repos/cnv-growth-bmi && git diff --stat 4b64136 ed0c509 -- AllOfUs/cnv_anthro
```

Scripts your finished stages used are unchanged except `conga_grid_qt.py` (new gene classes and one extra output column;
old designs give the same numbers) and `run_pipeline.py` (two new steps).

Add this to the `STAGES` dict in section 7 (and to each dataset's `stages` list if you want them in `status()`), rerun that cell, then launch:

```python
STAGES['fig3']        = '--steps fig3 --strata all --jobs 5'      # 14 fits + 2 gene-mirror + 2 dose-bin runs, then the in-Workbench summaries
STAGES['fig3_export'] = '--steps fig3_export --strata all'        # aggregate-only export -> <run>/export_aggregate_fig3
launch(['fig3'], tag='', ds='wgs')
```

- Always pass `--strata all` as above. Without it the pipeline rebuilds every stratum's phenotype file first.
- Output: `~/cnv_anthro_run/results/all/fig3/` (`size2t_mirror/`, `all2t_mirror/`, `totdose/`, `uni_N1_inf/`, `uni_N1_1/`, `uni_N2_inf/`,
  `gene_mirror/`, `dose_bins/`, `logs/`).
- The stage is restartable: a finished fit (it has `run_info.json`) is skipped. Your `launch()` backs the folder up as usual.
- Time and memory. In UK Biobank (N = 345,000, about 1 call per person) one `size2t` fit took 1 to 1.5 h, and 5 fits running at once
  peaked at 12 GB of memory. Your WGS table has far more calls, so expect several hours, as for your `gridjoint` stage. `--jobs 5`
  runs 5 fits at once; if `status()` shows memory pressure, stop, set `--jobs 3` and relaunch (finished fits are kept).
- The nosegdup dataset does not need this run: the segdup filter only affects the probe-level GWAS, and these fits read the same CNV and
  phenotype tables as the WGS dataset.

## B. QC (safe to paste; z, P and N only)

```python
import numpy as np
DS = 'wgs'; F3 = f"{DATASETS[DS]['run']}/results/all/fig3"
spec = json.load(open(f'{CODE}/designs/fig3_spec_2026-11-23.json'))
for f in sorted(glob.glob(f'{F3}/*/grid/*/*/run_info.json')):
    i = json.load(open(f)); t = i['trait'].replace('_rint', '')
    print(os.path.relpath(os.path.dirname(f), F3), '| N', i['N'], '| nrep', i['nrep'], '| shifts',
          'OK' if i['perms_sha1'] == spec['perms_sha1'][t] else 'MISMATCH', '| major loci', len(i['major_loci']))
    o = pd.read_csv(f.replace('run_info.json', 'observed_terms.tsv'), sep='\t')
    print('   (run, term) cells with 1-20 carriers, withheld at export:', int(pd.to_numeric(o.n_carriers_term, errors='coerce').between(1, 20).sum()),
          'of', int(o.n_carriers_term.notna().sum()))
for f in sorted(glob.glob(f'{F3}/uni_*/grid/*/univariate/grid_summary.tsv')):      # univariate DEL z / DUP z per class and CNV-size stratum
    d = pd.read_csv(f, sep='\t'); d = d[~d.readout.str.endswith('|mirror')]
    d['z'] = np.where(d.readout.str.contains('n_all'), d.t_analytic, d.z_emp)       # total genes: analytic t (its shift null is degenerate)
    print('---', os.path.relpath(os.path.dirname(f), F3)); print(d[['readout', 'z', 'P_emp', 'loo_held', 'loo_n', 'noloci_z_emp']].round(3).to_string(index=False))
for name in ['size2t_mirror', 'all2t_mirror', 'totdose']:
    for tab in ['mirror4_cross_summary.tsv', 'dose_summary.tsv']:
        p = f'{F3}/{name}/{tab}'
        if os.path.exists(p) and os.path.getsize(p) > 5:
            d = pd.read_csv(p, sep='\t'); d = d[d.stat.isin(['M4', 'armDEL', 'armDUP', 'D4', 'doseOth', 'doseTot', 'mirTot'])]
            print('---', name, tab); print(d[['trait', 'size', 'stratum', 'stat', 'z', 'P', 'loo_held', 'loo_n', 'noloci_z']].round(3).to_string(index=False))
for t in ['height', 'bmi']:
    d = pd.read_csv(f'{F3}/gene_mirror/{t}/gene_mirror_summary.tsv', sep='\t')
    print('--- gene mirror', t); print(d[['size', 'stratum', 'stat', 'obs', 'z', 'P', 'n_genes', 'loo_held', 'n_loo', 'noloci_z']].round(3).to_string(index=False))
for t in ['height', 'bmi']:      # dose-response bins: already masked when written (bins with 1-20 people have no n and no estimate)
    i = json.load(open(f'{F3}/dose_bins/{t}/run_info.json')); d = pd.read_csv(f'{F3}/dose_bins/{t}/fig3_dose_bins.tsv', sep='	')
    d = d[(d.run == 'full') & (d.side != 'none')]; est = d.beta.notna()
    print('--- dose bins', t, '| N', i['N'], '| single-dosage or non-carriers', i['full_N_single'], '| bins with an estimate', int(est.sum()), 'of', len(d),
          '| withheld', int((d.n.astype(str) == '<=20').sum()), '| empty', int((d.n.astype(str) == '0').sum()))
    print(d[est & (d.partition == 'uni')].groupby(['size', 'cls', 'side']).size().unstack('side').to_string())
```

`loo_held`, `loo_n` and `n_loo` are numbers of leave-one-locus-out runs (at most 9), and `n_genes` is a number of genes. None is a
participant count, so the 1–20 rule does not apply to them; your section 9.1 audit would flag them by name, which is why the Fig. 3 export
has its own audit.

Things to look at: `single-dosage or non-carriers` is the sample the dose-response bins use. On WGS, where most people carry gene-hitting
deletions and duplications, it may be a small part of N: report the number to Jon (open point 3). Every line says `shifts OK`; N equals your joint-run N for the trait; the number of withheld cells is small. If many
single-gene terms are withheld, tell Jon before exporting.

## C. Export and egress

```python
launch(['fig3_export'], tag='', ds='wgs')
```

When it is done, `~/cnv_anthro_run/logs/fig3_export.log` must end with `AUDIT OK`. The export is in `~/cnv_anthro_run/export_aggregate_fig3/`
(about 50 MB in UK Biobank). It holds, per fit, `observed_terms.tsv` (beta, SE and t of every term in every run), `null_runs.npz` (the t
of every term in each of the 2,000 shift replicates of every run), `run_info.json`; per trait `gene_z.npz` (per-gene score z) and the dose-response
bin tables (`fig3_dose_bins.tsv`, `fig3_dose_shape.tsv`: beta, SE, carriers and mean genes per bin); the design,
the spec and `EXPORT_MANIFEST.txt`. No summary table is exported: the meta-analysis recomputes them.

Small-cell rule, applied by the export step with your `MIN_CELL`:

- a term carried by 1–20 people in a run loses its beta, SE and t in that run, its null column and its mirror readout;
- a gene hit in 1–20 people (per run, CNV size and dosage) loses its z;
- gene carrier counts are never exported: the file holds 0 (not usable) or 21 (usable) only;
- counts of 1–20 calls or carriers in `run_info.json` are written as `<=20`;
- a dose-response bin with 1–20 people has no count and no estimate, and the bins of the run without the 34 intervals carry no counts
  at all (the difference between the two runs' counts could reveal a small number).

The export step refuses to overwrite an existing export folder: move the old one away if you need to export again.

Package and send, as in B4 of the first handoff (ask support for the egress size limit first if in doubt):

```python
DRUN = DATASETS[DS]['run']; T = f"{DATASETS[DS]['tar']}_fig3_export_{time.strftime('%Y%m%d')}.tar.gz"
sh(f"tail -2 {DRUN}/logs/fig3_export.log; head -3 {DRUN}/export_aggregate_fig3/all/fig3/EXPORT_MANIFEST.txt")
sh(f'cd {DRUN} && tar czf {T} export_aggregate_fig3 && sha256sum {T} > {T}.sha256 && ls -la {T}* && gsutil cp {T} {T}.sha256 {BUCKET_GS}/')
```

Then `scp` the tarball and its `.sha256` to `/expanse/projects/sebat1/jsebat/sacks2026_aou/aou_results/` and tell Jon. That is all we need back:

- `aou_sacks2026_fig3_export_<date>.tar.gz` + `.sha256` (WGS)
- `aou_array_sacks2026_fig3_export_<date>.tar.gz` + `.sha256` (arrays, section D)
- the pasted output of the QC cell for each dataset

## D. Array dataset

After the array main run is finished, rerun the array config cell of section 14 (it rebuilds `aou_sacks2026_array_filled.json` from the
template, which now has the `"fig3"` block), then:

```python
launch(['fig3'], tag='ar', ds='array')        # then QC with DS = 'array', then:
launch(['fig3_export'], tag='ar', ds='array')
```

Same QC, export and packaging with `DS = 'array'`.

## What happens at UCSD (for information)

- UK Biobank export (done): `/expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23/ukbb_export_55986701/ukbb_fig3_export.tar.gz`
  (51 MB; sha256 `fd9a29a3f51eb193442baadadbf28a572bd2d1b6d65425fb87b0d08f8c44be66`). It stays at UCSD.
- Tests on Expanse: a meta-analysis of the UK Biobank export alone reproduces the UK Biobank tables and the current figure, pixel for pixel.
  The code you will run (`run_pipeline.py --steps fig3,fig3_export`) was also refitted from scratch on UK Biobank data (1 h 52 min on
  16 cores for the joint fits and gene mirror) and gave the same tables as the finished runs. The univariate fits and bins took 1 h more.
- Meta-analysis: `AllOfUs/ukbb/meta_fig3_expanse.sbatch` → `AllOfUs/meta/meta_fig3.py`. Per term and run, Z = Σ √N·t / √(Σ N) across
  cohorts, formed in the observed data and in each of the 2,000 shared shift replicates; the mirror (M4), dose-aligned (D4), arm and total
  contrasts, the 9 locus-removal runs and `no_loci` are then recomputed from the meta t. The gene mirror combines the per-gene z the same way
  and recomputes the correlation and its null. The univariate z are combined like the joint-model terms; the dose-response bins by
  inverse-variance weighting, bin by bin.
- Figures: `GW2_FIG3_META=<meta folder> python AllOfUs/figures/figs_dosage_story.py D` (Fig. 3) and `... B` (dose-response curves of every
  stratum, one figure per trait).

## Open points (Jon decides)

1. **Call set on WGS.** The Fig. 3 fits use every call in your CNV table, as your grid and joint runs did. On WGS (median 57 calls per
   person) the total-genes term is then dominated by common calls. A version restricted to rare (< 1%) gene-hitting calls, like the
   `no_rare_gene` rule of the additivity runs, is not in this handoff. Do not add it unless Jon asks.
2. **Gene inclusion in the gene mirror.** UK Biobank uses genes with ≥ 5 deletion and ≥ 5 duplication carriers. Under the 1–20 rule an All
   of Us gene needs ≥ 21 carriers to contribute. Jon decides whether UK Biobank should use the same limit.
3. **Dose-response bins on WGS.** The bins use people with deletions only, duplications only, or neither. On WGS that may leave few
   people. If so, Jon decides between rare (< 1%) gene-hitting calls only, or the array dataset alone for these plots.
4. **Which All of Us dataset enters the main meta-analysis** (arrays, as for the other figures, with WGS reported separately) is unchanged
   from the first handoff.
