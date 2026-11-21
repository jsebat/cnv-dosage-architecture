# Handoff for Kelly Lee (and her Claude): replication for the two Sacks 2026 papers

Written 2026-11-20 by Jonathan Sebat's Claude session. You (Kelly Lee and Zhandong Liu) have a clinical cohort of CNV
carriers with age-matched non-carrier controls, mostly children. Jon would like you to rerun a defined set of analyses on
your data and send back summary tables. Nothing individual-level leaves your site.

There are two repositories, one per paper. Both are public; clone them, no access request is needed.

| Paper | Repository | Your part |
|---|---|---|
| Gene dosage architecture of anthropometric traits | `github.com/jsebat/cnv-dosage-architecture` (this one) | Part A below |
| BMI as a modifier of CNV effects on development and autistic traits | `github.com/jsebat/cnv-growth-bmi` | Part B below (`federated/` folder there, with its own `HANDOFF_Kelly.md`) |

## Part A (this repository): CNV main effects and the 22q11.2 subregions

**Goal.** Main effects of every recurrent CNV on height, BMI and head circumference, for all ages and split at 18 years,
as in Fig. 1 of the paper. Head circumference will be added to Fig. 1 from paediatric cohorts, and the under-18 versus
18-and-over comparison will illustrate age-dependent effects. Your 22q11.2 estimates also go into Fig. 6 (panel C, effects
of the nested 22q11.2 intervals) and into a pooled version of panel A.

**Steps**

1. **Genotype the recurrent CNVs.** `recurrentCNVgenotyping/` (Molly Sacks): `GenotypeRecurrentCNVs_tutorial.ipynb` and
   `README.txt`. Input: your QC'd CNV calls and a sample manifest. Output: one row per person, one column per locus
   (1 = deletion, 2 = no CNV, 3 = duplication). If your carriers come with a clinical genotype instead of array calls,
   build the same matrix from the clinical genotype; say so when you send results.
   For 22q11.2, keep the breakpoint class (A-B, A-C, A-D, B-D, C-D, distal) as separate loci. Fig. 6 needs them apart.
2. **Run the main effects.** `main_effects/README.md` has the exact commands.
   - Copy `main_effects/config_example.json`, set your column names, run `main_effects.R`.
   - Children (under 18) are scored against WHO norms for height and BMI and the UK 1990 reference for head
     circumference. Adults are Box-Cox scaled within your cohort with fixed exponents. This is done for you.
   - The primary estimate is carriers against your own controls, adjusted for age and sex.
   - Few adults is expected. Sparse strata come out as "suppressed" or "not estimable"; the run does not stop.
   - Run the tests first (`main_effects/tests/test_main_effects.R`) to check your R setup. They use synthetic data.
3. **Send back** the output folder of `main_effects.R` (aggregate TSV files only) and one line on how CNVs were
   genotyped, how many visits per person you had, and which visit was used (default: the latest with the outcome).

**What happens next (done at UCSD).** `main_effects/meta_main_effects.R` combines sites by fixed-effects meta-analysis;
`main_effects/plot_fig1_style.R` draws the Fig. 1 panels. The Fig. 6 test (is the A-D effect smaller than the sum of the
A-B and B-D effects?) is computed from the per-interval estimates you send.

**Data rules**
- Cells with fewer than 5 carriers are suppressed by the script (`min_carriers` in the config; raise it if your site
  requires).
- Do not send IDs, dates or any row-level table. The scripts write none.
- If something in your data does not fit the input format, change the config or your input table, not the model. If a
  model must change, tell us what and why, so all sites stay comparable.

## Part B: BMI modifier paper

See `HANDOFF_Kelly.md` and `federated/README.md` in `github.com/jsebat/cnv-growth-bmi`. In short: one script
(`federated/site_export.R`) runs the paper's non-genetic analyses (Figs. 1, 2, 3 and the growth and autistic-trait parts
of Fig. 5; not Fig. 4, which needs polygenic scores) and writes aggregates, including growth trajectories as spline
coefficients and binned means so that the pooled figure can be drawn without individual values.

## Things you should know before starting

- **Scoring is fixed.** Adults always Box-Cox; children WHO norms (UK 1990 for head circumference). Please do not
  substitute CDC charts or site-specific norms; a second scoring can be sent as a labelled sensitivity run.
- **Controls matter.** On 2026-11-20 we found that the child controls of our own growth analysis had been autistic
  SSC probands, not unaffected siblings, and corrected it. Your controls should be non-carriers without the CNVs under
  study; tell us how they were recruited and whether any have a neurodevelopmental diagnosis.
- **Relatedness.** If you have a family identifier, put it in the config; the BMI-paper script then reports
  family-clustered standard errors.
- **Which loci.** The BMI paper uses only 16p11.2 BP4-BP5 and 22q11.2 A-D. The dosage paper uses every recurrent locus.
- **Questions** go to Jon (jsebat). If your Claude finds a bug in these scripts, open an issue or a pull request on the
  repository; do not patch silently, because the other sites run the same code.
