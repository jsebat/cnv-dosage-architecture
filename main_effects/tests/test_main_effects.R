# test_main_effects.R  (2026-11-20)
# What: end-to-end test of the main_effects module on SYNTHETIC cohorts (tests/make_synthetic.R): generates the data,
#   runs main_effects.R through its command line for each cohort, then meta_main_effects.R and plot_fig1_style.R, and
#   checks that (1) the simulated effects are recovered within sampling error, (2) small cells are suppressed,
#   (3) the two estimators are labelled and never pooled without a flag, (4) the age-group tests find the simulated
#   difference, (5) the meta-analysis equals metafor (if installed), (6) one visit per person is used.
# Inputs: none. Outputs: PASS / FAIL lines; everything is written to a temporary folder outside the repository
#   (or to the folder given as first argument).
# null model: none; recovered = |estimate - simulated effect| < 4 SE + 0.06 (the 0.06 allows for the Box-Cox step and
#   for the covariate-adjusted adult scale).
# Usage: Rscript tests/test_main_effects.R [work_dir]
a <- commandArgs(FALSE); f <- gsub("~+~", " ", sub("^--file=", "", a[grep("^--file=", a)][1]), fixed = TRUE)
TESTS <- dirname(normalizePath(f)); MOD <- dirname(TESTS)
work <- commandArgs(TRUE)[1]; if (is.na(work)) work <- file.path(tempdir(), "main_effects_test")
unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
if (startsWith(normalizePath(work), normalizePath(MOD))) stop("choose a work folder outside the repository")
RS <- file.path(R.home("bin"), "Rscript"); nfail <- 0
check <- function(ok, what) { ok <- isTRUE(ok); cat(if (ok) "PASS " else "FAIL ", what, "\n"); if (!ok) nfail <<- nfail + 1; invisible(ok) }
run <- function(...) { st <- system2(RS, shQuote(c(...))); if (st != 0) stop("command failed: ", paste(..., collapse = " ")) }
rd <- function(p) read.delim(p, check.names = FALSE, stringsAsFactors = FALSE, na.strings = "NA")

MAIN_EFFECTS_SOURCE_ONLY <- TRUE; MAIN_EFFECTS_DIR <- MOD
source(file.path(TESTS, "make_synthetic.R")); source(file.path(MOD, "main_effects.R")); source(file.path(MOD, "meta_main_effects.R"))
syn <- make_all(work, file.path(MOD, "growth_refs_lms.csv")); truth <- syn$truth

# ---- unit checks
set.seed(1); X <- cbind(1, rnorm(200), rbinom(200, 1, 0.3), rnorm(200)); colnames(X) <- c("i", "a", "b", "c"); y <- rnorm(200) + X[, 3]
o <- ols(X, y); m <- summary(lm(y ~ X - 1))$coefficients
check(max(abs(o$b - m[, 1])) < 1e-10 && max(abs(o$se - m[, 2])) < 1e-10, "ols() equals lm() (coefficients and SE)")
pv <- pick_visit(c("p1", "p1", "p1", "p2", "p2"), c(3, 9, 12, 5, 4), c(TRUE, TRUE, FALSE, TRUE, TRUE))
check(setequal(pv, c(2, 4)), "visit rule: latest visit with the outcome, one per person")
cl <- close33(matrix(c(40, 37, 3, 20, 18, 2, 20, 19, 1), 3, 3, dimnames = list(AGES, SEXES)), 5)
check(all(cl$primary["ge18", ]) && all(rowSums(cl$primary | cl$complementary) != 1) && all(colSums(cl$primary | cl$complementary) != 1), "small-cell closure: no line with exactly one hidden cell")

# ---- run every cohort through the command line
for (l in names(syn$configs)) run(file.path(MOD, "main_effects.R"), "--config", syn$configs[[l]])
R <- lapply(names(syn$configs), function(l) rd(file.path(work, paste0("out_", l), sprintf("main_effects_%s.tsv", l)))); names(R) <- names(syn$configs)
AG <- lapply(names(syn$configs), function(l) rd(file.path(work, paste0("out_", l), sprintf("age_group_test_%s.tsv", l)))); names(AG) <- names(syn$configs)
tv <- function(locus, allele, trait, age) { r <- truth[truth$locus == locus & truth$allele == allele, ]
  if (trait == "bmi") { if (age == "lt18") r$bmi_child else if (age == "ge18") r$bmi_adult else NA } else r[[trait]] }
recovered <- function(d, label, est = "within-cohort", ages = c("all", "lt18", "ge18"), bcol = "beta_sd", scol = "se") {
  d <- d[d$estimator == est & d$stratum_age %in% ages & is.finite(d[[bcol]]), ]
  d$truth <- mapply(tv, d$locus, d$allele, d$trait, d$stratum_age); d <- d[!is.na(d$truth), ]
  dev <- abs(d[[bcol]] - d$truth) / (4 * d[[scol]] + 0.06)
  check(nrow(d) > 0 && all(dev < 1), sprintf("%s: %d %s estimates recover the simulated effects (worst |est - truth| = %.2f of the tolerance; median SE %.3f)", label, nrow(d), est, max(dev), median(d[[scol]])))
}
recovered(R$A, "cohort A"); recovered(R$C, "cohort C"); recovered(R$K, "cohort K (paediatric clinic with non-carriers)", ages = "lt18")
recovered(R$K, "cohort K", est = "reference-based", ages = "lt18"); recovered(R$B, "cohort B (carrier-only)", est = "reference-based", ages = "lt18")
check(sum(is.finite(R$A$beta_sd[R$A$estimator == "within-cohort"])) >= 6 * 3 * 9 - 9, "cohort A: all age x sex strata estimated for the six common alleles")

# ---- estimators labelled; carrier-only and few-adult cohorts degrade to 'not estimable' rather than failing
check(all(R$A$estimator %in% c("within-cohort", "reference-based")) && all(R$A$estimator_role[R$A$estimator == "within-cohort"] == "primary") && all(R$A$stratum_age[R$A$estimator == "reference-based"] == "lt18"),
      "every row labelled; within-cohort is primary; reference-based rows only for < 18 y")
b <- R$B[R$B$estimator == "within-cohort" & R$B$locus == "16p11.2_BP4_BP5" & R$B$allele == "DEL", ]
check(all(is.na(b$beta_sd)) && all(grepl("^not estimable", b$status)), "cohort B: within-cohort rows are 'not estimable' (no non-carriers)")
check(all(grepl("fallback", R$B$estimator_role[R$B$estimator == "reference-based"])), "cohort B: reference-based rows marked as fallback")
k <- R$K[R$K$estimator == "within-cohort" & R$K$locus == "16p11.2_BP4_BP5" & R$K$allele == "DEL" & R$K$stratum_age == "ge18" & R$K$stratum_sex == "all", ]
check(all(is.na(k$beta_sd)) && all(grepl("not estimable|suppressed", k$status)) && any(grepl("no adult head-circumference reference", k$status[k$trait == "hc"])),
      "cohort K: sparse adults give 'not estimable' rows; adult head circumference says there is no adult reference")
k <- R$K[R$K$estimator == "within-cohort" & R$K$stratum_age == "lt18" & R$K$stratum_sex == "all" & R$K$locus == "22q11.21_A_D" & R$K$allele == "DEL", ]
check(all(is.finite(k$beta_sd)) && all(grepl("age \\+ sex", k$covariates)), "cohort K: < 18 y within-cohort contrast estimated for all three traits, adjusted for age and sex")

# ---- one visit per person
kid <- with(syn$cohorts$B$pheno, tapply(age_years < 18 & !is.na(ofc_cm), subject, any)); g <- syn$cohorts$B$G[names(kid), "16p11.2_BP4_BP5"]
n_exp <- sum(kid & g == 1); n_got <- as.numeric(R$B$n_carriers[R$B$estimator == "reference-based" & R$B$locus == "16p11.2_BP4_BP5" & R$B$allele == "DEL" & R$B$trait == "hc" & R$B$stratum_age == "lt18" & R$B$stratum_sex == "all"])
check(n_got <= n_exp && n_got >= n_exp - 3 && nrow(syn$cohorts$B$pheno) > 2 * n_exp, sprintf("cohort B with repeated visits: %d carriers counted once each (%d persons with a child measurement)", n_got, n_exp))

# ---- suppression
for (l in names(R)) { n <- suppressWarnings(as.numeric(c(R[[l]]$n_carriers, R[[l]]$n_noncarriers, AG[[l]]$n_carriers_lt18, AG[[l]]$n_carriers_ge18, AG[[l]]$n_carriers_in_both_strata)))
  check(!any(n > 0 & n < 5, na.rm = TRUE), sprintf("cohort %s: no count between 1 and 4 in any output", l)) }
s <- R$A[R$A$locus == "15q13.3_BP4_BP5" & R$A$allele == "DEL" & R$A$stratum_age == "all" & R$A$stratum_sex == "all", ]
check(all(s$n_carriers == "<5") && all(is.na(s$beta_sd)) && all(is.na(s$mean_age_carriers)) && all(grepl("^suppressed", s$status)), "cohort A: the 3-carrier deletion is suppressed (count shown as <5, no statistics)")
small <- R$K[R$K$n_carriers == "<5", ]; check(nrow(small) > 0 && all(is.na(small$beta_sd)) && all(is.na(small$p)), "cohort K: every <5 cell has no statistics")
check(any(R$K$n_carriers == ">=5"), "cohort K: complementary cells have their count masked (>=5)")
hdr <- unlist(lapply(list.files(work, pattern = "^(main_effects|age_group_test|scoring_info)_.*tsv$", recursive = TRUE, full.names = TRUE), function(p) names(rd(p))))
txt <- unlist(lapply(list.files(work, pattern = "^(main_effects|age_group_test|scoring_info|run_log)_", recursive = TRUE, full.names = TRUE), readLines))
check(!any(tolower(hdr) %in% c("id", "iid", "fid", "identifier", "subject", "sampleid", "sample_id")) && !any(grepl("FAKE_", txt)), "outputs carry no id column and no participant id")
# strict mode as needed for All of Us: minimum 20, complementary cells suppressed outright
strict <- run_main_effects(list(cohort_label = "Kstrict", min_carriers = 20, complementary = "suppress", pheno = list(id = "subject", sex = "gender", age = "age_years", hc_cm = "ofc_cm")),
                           geno = syn$cohorts$K$geno, pheno = syn$cohorts$K$pheno, write = FALSE)$main
n <- suppressWarnings(as.numeric(strict$n_carriers)); cs <- grepl("complementary", strict$status)
check(!any(n > 0 & n < 20, na.rm = TRUE) && any(cs) && all(is.na(strict$beta_sd[cs])) && !any(strict$n_carriers == ">=20" & is.finite(strict$beta_sd)), "strict mode (min 20, complementary = suppress): no count 1-19, complementary cells carry no statistics")

# ---- adult scaling: fixed UK Biobank exponents by default; "estimate" gives per-cohort exponents and near-identical effects
si <- rd(file.path(work, "out_A", "scoring_info_A.tsv")); ad <- si[si$age_group == "ge18", ]; lam <- function(t, s) ad$lambda[ad$trait == t & ad$sex == s]
check(lam("height", "male") == 0.65 && lam("height", "female") == 0.60 && lam("bmi", "male") == -0.75 && lam("bmi", "female") == -1.05 && all(ad$lambda[ad$trait == "hc"] == 1) &&
      all(grepl("UK Biobank", ad$lambda_source[ad$trait != "hc"])), "adults: fixed Box-Cox exponents (height 0.65 / 0.60, BMI -0.75 / -1.05, head circumference untransformed)")
pcA <- list(id = "subject", sex = "gender", age = "age_years", height_cm = "height_cm", weight_kg = "weight_kg", hc_cm = "ofc_cm")
est <- run_main_effects(list(cohort_label = "A", adult_boxcox = "estimate", pheno = pcA), geno = syn$cohorts$A$geno, pheno = syn$cohorts$A$pheno, write = FALSE)
e1 <- est$main[est$main$stratum_age == "ge18" & est$main$estimator == "within-cohort", ]; f1 <- R$A[R$A$stratum_age == "ge18" & R$A$estimator == "within-cohort", ]
ok <- is.finite(e1$beta_sd) & is.finite(f1$beta_sd); el <- est$scoring[est$scoring$age_group == "ge18" & est$scoring$trait == "bmi", ]
check(all(e1$locus == f1$locus) && max(abs(e1$beta_sd[ok] - f1$beta_sd[ok])) < 0.1 && all(grepl("estimated", el$lambda_source)) && all(abs(el$lambda) < 0.6),
      sprintf("adult_boxcox = 'estimate': exponent estimated per cohort (simulated BMI is log-normal: %.2f, %.2f); adult effects within %.3f SD of the fixed-exponent run", el$lambda[1], el$lambda[2], max(abs(e1$beta_sd[ok] - f1$beta_sd[ok]))))
set.seed(5); yy <- rnorm(400, 100, 15)^(1 / 4)   # needs an exponent near 4, beyond the grid
wmsg <- tryCatch({ score_adult(yy * 10, runif(400, 20, 60), rep(c("male", "female"), 200), rep(TRUE, 400), "bmi", "estimate"); "" }, warning = function(w) conditionMessage(w))
check(grepl("edge of the search grid", wmsg), "adult_boxcox = 'estimate' warns when the exponent is on the edge of the grid")

# ---- age-group tests (16p11.2 deletion BMI simulated +0.3 in children, +0.9 in adults)
ag <- AG$A[AG$A$locus == "16p11.2_BP4_BP5" & AG$A$allele == "DEL" & AG$A$trait == "bmi" & AG$A$stratum_sex == "all" & AG$A$estimator_lt18 == "within-cohort", ]
w <- ag[grepl("^Wald", ag$test), ]; i <- ag[grepl("^interaction", ag$test), ]
check(nrow(w) == 1 && abs(w$diff_lt18_minus_ge18 + 0.6) < 4 * w$se_diff + 0.06 && w$p_diff < 1e-4, sprintf("age-group Wald test finds the simulated difference (%.2f, SE %.2f, P %.1e; simulated -0.60)", w$diff_lt18_minus_ge18, w$se_diff, w$p_diff))
check(nrow(i) == 1 && abs(i$diff_lt18_minus_ge18 + 0.6) < 4 * i$se_diff + 0.06 && i$p_diff < 1e-4, sprintf("interaction term finds it too (%.2f, SE %.2f, P %.1e)", i$diff_lt18_minus_ge18, i$se_diff, i$p_diff))
nul <- AG$A[AG$A$trait == "height" & AG$A$stratum_sex == "all" & AG$A$estimator_lt18 == "within-cohort" & is.finite(AG$A$p_diff), ]
check(nrow(nul) >= 10 && min(nul$p_diff) > 1e-4, sprintf("height (no simulated age difference): smallest of %d age-difference P values is %.3f", nrow(nul), min(nul$p_diff)))
check(any(grepl("mixed estimators", AG$A$note)) && all(AG$A$note[AG$A$estimator_lt18 == "reference-based"] != ""), "age tests that mix estimators are flagged")

# ---- meta-analysis
mdir <- file.path(work, "meta"); inp <- file.path(work, paste0("out_", names(R)), sprintf("main_effects_%s.tsv", names(R))); ainp <- sub("main_effects_", "age_group_test_", inp)
run(file.path(MOD, "meta_main_effects.R"), "--in", inp, "--age_in", ainp, "--out", mdir)
M <- rd(file.path(mdir, "meta_main_effects.tsv")); MA <- rd(file.path(mdir, "meta_age_group_test.tsv"))
recovered(M[M$estimator != "reference-based", ], "meta-analysis", est = "within-cohort", bcol = "beta_fe", scol = "se_fe")
wi <- M[M$estimator == "within-cohort", ]
check(!any(grepl("(^|;)B(;|$)", wi$cohorts)) && all(!wi$mixed_estimators) && all(wi$cohorts_reference_based == 0), "meta: within-cohort rows never include the carrier-only cohort")
mx <- M[M$mixed_estimators, ]
check(nrow(mx) > 0 && all(grepl("^mixed", mx$estimator)) && all(mx$cohorts_reference_based > 0 & mx$cohorts_within_cohort > 0) && all(mx$stratum_age == "lt18"), "meta: pooled rows that combine both estimators are labelled 'mixed' and flagged")
one <- M[M$locus == "16p11.2_BP4_BP5" & M$allele == "DEL" & M$trait == "hc" & M$stratum_age == "lt18" & M$stratum_sex == "all", ]
check(setequal(one$estimator, c("within-cohort", "reference-based", MIXED)) && one$k[one$estimator == "within-cohort"] == 3 && one$k[one$estimator == MIXED] == 4, "meta: 16p11.2 deletion head circumference < 18 y has k = 3 within-cohort and k = 4 mixed")
if (requireNamespace("metafor", quietly = TRUE)) {
  x <- do.call(rbind, lapply(R, function(d) d[d$locus == "22q11.21_A_D" & d$allele == "DEL" & d$trait == "height" & d$stratum_age == "lt18" & d$stratum_sex == "all" & d$estimator == "within-cohort", ]))
  mf <- metafor::rma(yi = x$beta_sd, sei = x$se, method = "FE"); mm <- M[M$locus == "22q11.21_A_D" & M$allele == "DEL" & M$trait == "height" & M$stratum_age == "lt18" & M$stratum_sex == "all" & M$estimator == "within-cohort", ]
  check(abs(mm$beta_fe - as.numeric(mf$beta)) < 1e-5 && abs(mm$se_fe - mf$se) < 1e-5 && abs(mm$Q - mf$QE) < 1e-4 && abs(mm$I2 - mf$I2 / 100) < 1e-4, "meta: beta, SE, Q and I2 equal metafor::rma(method = 'FE')")
} else cat("SKIP  metafor not installed: comparison with metafor::rma not run\n")
ma <- MA[MA$locus == "16p11.2_BP4_BP5" & MA$allele == "DEL" & MA$trait == "bmi" & MA$stratum_sex == "all" & MA$estimator_lt18 == "within-cohort", ]
check(nrow(ma) == 2 && all(abs(ma$diff_lt18_minus_ge18 + 0.6) < 4 * ma$se_diff + 0.06) && all(ma$p_diff < 1e-6), sprintf("meta age-group comparison: pooled difference %.2f (Wald) and %.2f (pooled interaction terms); simulated -0.60", ma$diff_lt18_minus_ge18[1], ma$diff_lt18_minus_ge18[2]))

# ---- figures
pdir <- file.path(work, "figures"); run(file.path(MOD, "plot_fig1_style.R"), "--main", file.path(mdir, "meta_main_effects.tsv"), "--age", file.path(mdir, "meta_age_group_test.tsv"), "--out", pdir)
pf <- file.path(pdir, c("fig1_style_main_effects.png", "fig1_style_age_groups.png")); check(all(file.exists(pf)) && all(file.size(pf) > 20000), "both figures written")
run(file.path(MOD, "plot_fig1_style.R"), "--main", inp[names(R) == "K"], "--age", ainp[names(R) == "K"], "--out", file.path(work, "figures_K"))
check(file.exists(file.path(work, "figures_K", "fig1_style_main_effects.png")), "figure also drawn from a single cohort's table")

cat(sprintf("\n%s  (%d failed)   work folder: %s\n", if (nfail == 0) "ALL TESTS PASSED" else "TESTS FAILED", nfail, normalizePath(work)))
quit(status = if (nfail == 0) 0 else 1)
