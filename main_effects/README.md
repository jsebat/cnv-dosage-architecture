# main_effects: recurrent CNV main effects on height, BMI and head circumference (federated)

Written 2026-11-20. Estimates, inside one cohort, the effect of every recurrent CNV (deletion and duplication at each
locus of Supplementary Table 1) on height, BMI and head circumference, as in Fig. 1 of Sacks et al. 2026, and the same
effects split by age (< 18 y and >= 18 y). The site that holds the data runs the scripts and sends back **aggregate
tables only**; the tables from all sites are then combined by fixed-effects meta-analysis.

The purpose of the two additions to the original Fig. 1 analysis: head circumference can be added to the figure (marked
as coming from paediatric cohorts), and age-dependent effects can be shown by comparing < 18 y with >= 18 y.

## What you need

- R (>= 4.0) with the packages `MASS` and `jsonlite` (`data.table` speeds up reading; `ggplot2` only for the figures;
  `metafor` only for one optional test).
- This folder, complete: `main_effects.R`, `00_scoring.R` and `growth_refs_lms.csv` must stay together.
- **Genotype matrix** from `../recurrentCNVgenotyping/GenotypeRecurrentCNVs_tutorial.ipynb`
  (`genomatrix.recurrentCNV.<date>.csv`): one row per sample, an `identifier` column, one column per locus coded
  1 = deletion, 2 = no CNV, 3 = duplication, empty = sample failed QC. Locus columns are found automatically (cytoband-like
  names with only 1/2/3 values; the `_oneway` and helper columns of the notebook are ignored). A cohort whose carriers
  are known from clinical testing can supply any table in this format.
- **Phenotype table** (tab- or comma-separated, or .rds), one row per measurement visit:

  | column | required | notes |
  |---|---|---|
  | id | yes | must match the genotype matrix `identifier` |
  | sex | yes | any coding; give the male / female values in the config |
  | age at measurement | yes | years (or months / days: set `age_unit`) |
  | height (cm) | per trait | |
  | weight (kg) or BMI | per trait | BMI is computed from weight and height when no BMI column is given |
  | head circumference (cm) | per trait | |
  | PC1 ... PC10 | optional | genetic principal components; used as covariates when present |
  | site | optional | added as a covariate (factor) |
  | cohort | optional | if given, each cohort value is analysed separately |

  Repeated visits are welcome: put each visit on its own row (see the visit rule below).
- **Config**: copy `config_example.json` and fill in the two file paths and your column names.

## Commands

```
# 1. at each site
Rscript main_effects.R --config my_config.json
# 2. at the coordinating site, on the tables returned by the sites
Rscript meta_main_effects.R --in siteA/main_effects_A.tsv siteB/main_effects_B.tsv \
        --age_in siteA/age_group_test_A.tsv siteB/age_group_test_B.tsv --out meta_out
# 3. figures (from the meta-analysis, or from one site's main_effects_<cohort>.tsv)
Rscript plot_fig1_style.R --main meta_out/meta_main_effects.tsv --age meta_out/meta_age_group_test.tsv --out figures
# test on synthetic data (about 1 minute; writes to a temporary folder)
Rscript tests/test_main_effects.R
```

`main_effects.R` also accepts `--geno`, `--pheno`, `--out`, `--cohort` and `--min_carriers` to override the config.

## What to send back

The four files in `out_dir`, and nothing else:

| file | content |
|---|---|
| `main_effects_<cohort>.tsv` | the effects (schema below) |
| `age_group_test_<cohort>.tsv` | < 18 y vs >= 18 y difference per locus, allele and trait |
| `scoring_info_<cohort>.tsv` | per trait, age group and sex: how many non-carriers were scored, their mean and SD Z, the adult Box-Cox exponent and its source |
| `run_log_<cohort>.json` | settings, counts of excluded rows, list of loci |

## Data rules

- No individual-level data leave the site. The outputs contain no ids, no dates and no individual measurements;
  please look at them before sending.
- Small cells: any locus x allele x trait x stratum with fewer than `min_carriers` carriers (default 5) is suppressed:
  the count is written as `<5` and the effect, SE, P and mean age are empty. All of Us requires 20; the All of Us export
  in this repository hides counts of 1 to 20, so use `"min_carriers": 21, "complementary": "suppress"` there.
- Complementary cells: when a hidden count could be recovered by subtraction (for example all ages minus < 18 y), a
  second count in the same line is masked as `>=5`; its effect estimate is kept. With `"complementary": "suppress"`
  the estimate is removed too.
- Never add input or output files of real cohorts to this repository; the tests use simulated data with `FAKE_` ids.

## Output schema: `main_effects_<cohort>.tsv`

One row per locus x allele x trait x age stratum x sex stratum x estimator.

| column | meaning |
|---|---|
| cohort, locus | cohort label; locus name as in the genotype matrix |
| allele | DEL or DUP |
| trait | height, bmi, hc (head circumference) |
| stratum_age | all, lt18 (< 18 y), ge18 (>= 18 y) |
| stratum_sex | all, male, female |
| estimator | `within-cohort` or `reference-based` (see below) |
| estimator_role | primary; sensitivity; fallback (no non-carriers) |
| n_carriers, n_noncarriers | persons (one visit each); `<5` = suppressed, `>=5` = masked |
| beta_sd, se, t, df, p | effect in SD units of the Z score, standard error, t statistic, degrees of freedom, two-sided P |
| mean_age_carriers | mean age (years) of the carriers at the visit used |
| scoring | how the trait was turned into a Z score in this stratum |
| covariates | covariates in the model |
| status | ok; suppressed; no carriers; not estimable (with the reason) |
| min_carriers | the minimum carrier count that was applied |

`age_group_test_<cohort>.tsv`: per locus x allele x trait x sex stratum, `diff_lt18_minus_ge18`, `se_diff`, `p_diff`,
the two stratum estimates, the number of carriers in each stratum and in both, for two tests: the Wald test on the two
stratum estimates and the interaction term of one model.

## Known choices

- **Scoring, children (< 18 y).** Age- and sex-normed Z scores from growth references: WHO 2006 (< 5 y) and WHO 2007
  (5 to 18 y) for height and BMI, UK 1990 for head circumference; |Z| > 5 is set to missing. The scorer and the LMS table
  are copied from the companion repository (`00_scoring.R`, `growth_refs_lms.csv`; exported from the R package
  childsds 0.9.11). Measurements outside wide plausible ranges (height 30-250 cm, weight 1-400 kg, BMI 5-100, head
  circumference 20-80 cm; `plausible` in the config) are set to missing first.
- **Scoring, adults (>= 18 y).** There is no growth reference, so adults are scored relative to their own cohort,
  with the adult scorer of the companion BMI paper (`score_adult_bc()` in `00_scoring.R`, called unchanged so the two
  repositories agree). Within sex: Box-Cox transform; then `lm(transformed y ~ age + age^2)` fitted in the cohort's
  non-carrier adults; Z = residual / non-carrier residual SD (|Z| > 10 set to missing). Carriers are scored with the
  non-carrier model. At least `min_adult_ref` (default 50) non-carrier adults per sex are needed; otherwise the adult
  rows read "not estimable".
  - `"adult_boxcox": "fixed_ukbb"` (default): the Box-Cox exponent is fixed at the value estimated once in UK Biobank
    non-carriers: height 0.65 (men) and 0.60 (women), BMI -0.75 (men) and -1.05 (women). Estimating it in a small or
    clinical cohort is unstable (in the lab's clinical table the height estimate ran to the edge of the search grid, 2.0).
  - `"adult_boxcox": "estimate"`: the exponent is estimated in the cohort's own non-carrier adults (profile likelihood,
    grid -2 to 2). Only sensible for biobank-size cohorts; a warning is printed when the estimate is on the grid edge.
  - `scoring_info` reports the exponent used and where it came from.
  - Difference from the companion repository: adults start at 18 y here (the < 18 / >= 18 split of this module); there
    height and BMI use the adult scale from 19 y.
- **Adult head circumference.** No adult reference exists, and UK Biobank has no head circumference, so it is not
  transformed (exponent 1 in both modes). It is estimated on the cohort-relative scale above when the cohort has
  enough non-carrier adults with head circumference; otherwise the row says so. Adult head-circumference
  effects are therefore relative to the cohort, not to a reference.
- **Few adults are expected and fine.** In a paediatric cohort the >= 18 y strata will mostly read "suppressed" or
  "not estimable". The script does not stop; the < 18 y rows are the ones that matter, and the age-group comparison is
  then made at the meta-analysis stage against cohorts that do have adults.
- **Non-carriers.** A non-carrier is a person with no CNV at any locus in the genotype matrix (a carrier of another
  recurrent CNV is not used as a control). Persons whose genotype is missing at a locus are left out for that locus.
- **Two estimators, labelled on every row.**
  - `within-cohort` (primary): carriers against the cohort's own non-carriers, both scored the same way,
    `lm(Z ~ genotype (reference = no CNV) + age + sex + PC1..PC10 [+ site])`, deletion and duplication in one model.
    Age and sex are included because clinical carriers and controls can differ in age and sex even after norming; set
    `"adjust_age_sex": false` to get the manuscript's model `Z ~ genotype + PCs` exactly. The all-ages model also has an
    age-group term (children and adults are on different Z scales).
  - `reference-based` (sensitivity; children only): the mean Z of the carriers against the growth reference, one-sample
    t-test. It is the only estimate available for a clinic with no non-carriers (role "fallback"). It assumes the
    cohort's children would average Z = 0 without the CNV; `scoring_info` shows how far the cohort's own non-carriers
    are from that. Set `"reference_based": false` to leave these rows out.
  The meta-analysis pools the two estimators separately. A row that combines them is labelled `mixed` and has
  `mixed_estimators = TRUE`.
- **Visit rule.** One visit per person for each trait and stratum: the latest visit with a valid value of that trait
  (`"visit_rule": "earliest"` for the first). The < 18 y stratum uses the person's latest visit before 18, the >= 18 y
  stratum the latest visit from 18 on, and the all-ages stratum the adult visit when there is one. A person measured
  both as a child and as an adult therefore appears in both age strata; `n_carriers_in_both_strata` reports how many.
- **Age-group test.** (a) Wald test of (effect < 18 y) minus (effect >= 18 y) from the two stratum estimates; it treats
  the strata as independent, which is conservative when the same persons are in both. (b) The genotype x age-group
  interaction term of one model in which each person appears once. The meta-analysis repeats (a) on the pooled
  estimates and pools the per-cohort terms of (b).
- **Meta-analysis.** Fixed-effects, inverse-variance weights, Cochran Q and I2, the same numbers as
  `metafor::rma(method = "FE")` (checked in the test).
- **Figures.** Height and BMI are drawn from the all-ages stratum, head circumference from the < 18 y stratum, with an
  asterisk and the footnote "paediatric clinical cohorts".

## Files

| file | purpose |
|---|---|
| `main_effects.R` | per-cohort estimates (run at each site) |
| `meta_main_effects.R` | meta-analysis across cohorts and pooled age-group comparison |
| `plot_fig1_style.R` | Fig. 1 style effect plot and < 18 y vs >= 18 y plot |
| `00_scoring.R`, `growth_refs_lms.csv` | growth-reference scorer and LMS table (WHO 2006, WHO 2007, UK 1990) |
| `config_example.json` | config template |
| `tests/make_synthetic.R`, `tests/test_main_effects.R` | simulated cohorts with known effects; end-to-end test |
