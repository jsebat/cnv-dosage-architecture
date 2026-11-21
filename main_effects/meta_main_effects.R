# meta_main_effects.R  (2026-11-20)
# What: fixed-effects inverse-variance meta-analysis of the per-cohort recurrent-CNV main effects written by
#   main_effects.R, by locus x allele x trait x age stratum x sex stratum, with Cochran Q and I2 (the same numbers as
#   metafor::rma(yi, sei = , method = "FE"); checked in tests/), and the < 18 vs >= 18 comparison on the pooled estimates.
#   The two estimators are never pooled silently: "within-cohort" and "reference-based" rows are meta-analysed
#   separately, and a third row, estimator = "mixed", pools one estimate per cohort (within-cohort when the cohort has
#   one, otherwise reference-based) with mixed_estimators = TRUE. It is written only when it really mixes the two.
# Inputs: --in  one or more main_effects_<cohort>.tsv (files or folders)
#         --age_in (optional) age_group_test_<cohort>.tsv files, to pool the per-cohort interaction terms
#         --out output folder
# Outputs: meta_main_effects.tsv, meta_age_group_test.tsv
# null model: none; pooled effect = Wald z-test of the inverse-variance weighted mean; heterogeneity = Cochran Q
#   (chi-square, k - 1 df); age-group difference = Wald z-test of (pooled < 18) - (pooled >= 18).
# Usage: Rscript meta_main_effects.R --in siteA/main_effects_A.tsv siteB/main_effects_B.tsv --out meta_out

ivw <- function(b, s) {
  ok <- is.finite(b) & is.finite(s) & s > 0; b <- b[ok]; s <- s[ok]; k <- length(b)
  if (!k) return(NULL)
  w <- 1 / s^2; bf <- sum(w * b) / sum(w); sf <- sqrt(1 / sum(w)); Q <- if (k > 1) sum(w * (b - bf)^2) else 0
  data.frame(k = k, beta_fe = bf, se_fe = sf, z = bf / sf, p_fe = 2 * pnorm(-abs(bf / sf)), Q = if (k > 1) Q else NA, Q_df = k - 1,
             p_Q = if (k > 1) pchisq(Q, k - 1, lower.tail = FALSE) else NA, I2 = if (k > 1) (if (Q > 0) max(0, (Q - (k - 1)) / Q) else 0) else NA)
}
KEY <- c("locus", "allele", "trait", "stratum_age", "stratum_sex")
MIXED <- "mixed (within-cohort + reference-based)"

read_many <- function(paths, pattern) {
  f <- unlist(lapply(paths, function(p) if (dir.exists(p)) list.files(p, pattern = pattern, full.names = TRUE) else p))
  if (!length(f)) stop("no input files")
  do.call(rbind, lapply(f, function(x) read.delim(x, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character", na.strings = "NA")))
}

meta_main <- function(D) {
  D$beta_sd <- suppressWarnings(as.numeric(D$beta_sd)); D$se <- suppressWarnings(as.numeric(D$se))
  D <- D[is.finite(D$beta_sd) & is.finite(D$se) & D$se > 0, , drop = FALSE]
  dup <- duplicated(D[, c("cohort", KEY, "estimator")]); if (any(dup)) stop("the same cohort label appears twice for one cell: ", paste(unique(D$cohort[dup]), collapse = ", "))
  nnum <- suppressWarnings(as.numeric(D$n_carriers)); D$n_lb <- ifelse(is.na(nnum), suppressWarnings(as.numeric(sub("^>=", "", D$n_carriers))), nnum); D$n_masked <- is.na(nnum)
  one <- function(g, est, mixed) { r <- ivw(g$beta_sd, g$se); if (is.null(r)) return(NULL)
    cbind(g[1, KEY], estimator = est, mixed_estimators = mixed, r, cohorts = paste(g$cohort, collapse = ";"),
          cohorts_within_cohort = sum(g$estimator == "within-cohort"), cohorts_reference_based = sum(g$estimator == "reference-based"),
          n_carriers_total = sum(g$n_lb), n_carriers_is_lower_bound = any(g$n_masked), scoring = paste(unique(g$scoring), collapse = " | "), row.names = NULL) }
  out <- list()
  for (g in split(D, D[, KEY], drop = TRUE)) {
    w <- g[g$estimator == "within-cohort", ]; r <- g[g$estimator == "reference-based", ]
    if (nrow(w)) out[[length(out) + 1]] <- one(w, "within-cohort", FALSE)
    if (nrow(r)) out[[length(out) + 1]] <- one(r, "reference-based", FALSE)
    extra <- r[!r$cohort %in% w$cohort, ]
    if (nrow(w) && nrow(extra)) out[[length(out) + 1]] <- one(rbind(w, extra), MIXED, TRUE)
  }
  M <- do.call(rbind, out); M[order(M$trait, M$locus, M$allele, M$stratum_age, M$stratum_sex, M$estimator), ]
}

meta_age <- function(M) {
  ch <- M[M$stratum_age == "lt18", ]; ad <- M[M$stratum_age == "ge18" & M$estimator == "within-cohort", ]
  m <- merge(ch, ad, by = c("locus", "allele", "trait", "stratum_sex"), suffixes = c("_lt18", "_ge18"))
  if (!nrow(m)) return(NULL)
  d <- m$beta_fe_lt18 - m$beta_fe_ge18; se <- sqrt(m$se_fe_lt18^2 + m$se_fe_ge18^2)
  data.frame(m[, c("locus", "allele", "trait", "stratum_sex")], test = "Wald (difference of the two pooled stratum estimates)", estimator_lt18 = m$estimator_lt18, estimator_ge18 = m$estimator_ge18,
             mixed_estimators = m$estimator_lt18 != "within-cohort", k_lt18 = m$k_lt18, k_ge18 = m$k_ge18, n_carriers_lt18 = m$n_carriers_total_lt18, n_carriers_ge18 = m$n_carriers_total_ge18,
             beta_lt18 = m$beta_fe_lt18, se_lt18 = m$se_fe_lt18, beta_ge18 = m$beta_fe_ge18, se_ge18 = m$se_fe_ge18,
             diff_lt18_minus_ge18 = d, se_diff = se, z = d / se, p_diff = 2 * pnorm(-abs(d / se)), Q = NA, p_Q = NA, I2 = NA, cohorts = paste(m$cohorts_lt18, m$cohorts_ge18, sep = " / "), stringsAsFactors = FALSE)
}

meta_interaction <- function(A) {   # pools the per-cohort interaction terms (each from one model within one cohort)
  A <- A[grepl("^interaction", A$test), ]; A$d <- suppressWarnings(as.numeric(A$diff_lt18_minus_ge18)); A$s <- suppressWarnings(as.numeric(A$se_diff))
  A <- A[is.finite(A$d) & is.finite(A$s) & A$s > 0, ]; if (!nrow(A)) return(NULL)
  do.call(rbind, lapply(split(A, A[, c("locus", "allele", "trait", "stratum_sex")], drop = TRUE), function(g) { r <- ivw(g$d, g$s)
    data.frame(g[1, c("locus", "allele", "trait", "stratum_sex")], test = "interaction term, pooled across cohorts", estimator_lt18 = "within-cohort", estimator_ge18 = "within-cohort",
               mixed_estimators = FALSE, k_lt18 = r$k, k_ge18 = r$k, n_carriers_lt18 = NA, n_carriers_ge18 = NA, beta_lt18 = NA, se_lt18 = NA, beta_ge18 = NA, se_ge18 = NA,
               diff_lt18_minus_ge18 = r$beta_fe, se_diff = r$se_fe, z = r$z, p_diff = r$p_fe, Q = r$Q, p_Q = r$p_Q, I2 = r$I2, cohorts = paste(g$cohort, collapse = ";"), row.names = NULL, stringsAsFactors = FALSE) }))
}

run_meta <- function(inputs, out_dir, age_inputs = NULL) {
  M <- meta_main(read_many(inputs, "^main_effects_.*\\.tsv$")); A <- meta_age(M)
  if (!is.null(age_inputs)) A <- rbind(A, meta_interaction(read_many(age_inputs, "^age_group_test_.*\\.tsv$")))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(M, file.path(out_dir, "meta_main_effects.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  if (!is.null(A)) write.table(A, file.path(out_dir, "meta_age_group_test.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  message("wrote ", nrow(M), " pooled cells to ", normalizePath(out_dir)); invisible(list(main = M, age = A))
}

if (!exists("MAIN_EFFECTS_SOURCE_ONLY") && !interactive()) {
  a <- commandArgs(TRUE); grab <- function(k) { i <- match(k, a); if (is.na(i)) return(NULL); j <- i + 1; e <- j
    while (e <= length(a) && !grepl("^--", a[e])) e <- e + 1; if (e == j) NULL else a[j:(e - 1)] }
  if (is.null(grab("--in")) || is.null(grab("--out"))) stop("usage: Rscript meta_main_effects.R --in main_effects_A.tsv main_effects_B.tsv [--age_in age_group_test_A.tsv ...] --out DIR")
  run_meta(grab("--in"), grab("--out"), grab("--age_in"))
}
