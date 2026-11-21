# main_effects.R  (2026-11-20)
# What: main effect of every recurrent CNV (deletion and duplication at each locus) on height, BMI and head
#   circumference within ONE cohort, for all ages, < 18 y and >= 18 y, each also by sex, plus the CNV x age-group test.
#   Federated use: run at the site that holds the data; only the aggregate tables written here leave the site.
# Inputs: (1) genotype matrix from recurrentCNVgenotyping/GenotypeRecurrentCNVs_tutorial.ipynb (one row per sample,
#   one column per locus, 1 = DEL, 2 = no CNV, 3 = DUP, empty = failed QC); (2) phenotype table (id, sex, age at
#   measurement, height cm, weight kg or BMI, head circumference cm; optional PC1-PC10, site, repeated visits);
#   (3) a JSON config naming the files and columns (config_example.json). Needs 00_scoring.R + growth_refs_lms.csv
#   in the same folder. R packages: MASS, jsonlite (data.table is used for reading if installed).
# Outputs (aggregates only, in out_dir): main_effects_<cohort>.tsv, age_group_test_<cohort>.tsv,
#   scoring_info_<cohort>.tsv, run_log_<cohort>.json. Cells under the minimum carrier count are suppressed.
# Scoring: children (< 18 y) = age- and sex-normed Z from growth references (WHO height and BMI, UK 1990 head
#   circumference), |Z| > 5 set to missing. Adults = cohort-relative Z as in the companion BMI paper: within sex, Box-Cox
#   with a fixed exponent (estimated once in UK Biobank non-carriers; head circumference untransformed), then residual
#   of lm(y ~ age + age^2) / residual SD, both from the cohort's non-carrier adults (adult_boxcox in the config).
# Estimators (every output row is labelled):
#   "within-cohort"   (primary) carriers vs the cohort's own non-carriers:
#                     lm(Z ~ factor(CNV genotype, ref = no CNV) + sex + age + PC1..PC10 [+ site])
#   "reference-based" (sensitivity; fallback for carrier-only clinics; children only) mean Z of carriers against the
#                     growth reference, one-sample.
# null model: none; within-cohort = OLS t-test of the genotype coefficient; reference-based = one-sample t-test of
#   the carrier mean Z against 0; age-group difference = Wald z-test of (effect < 18) - (effect >= 18) and the OLS
#   t-test of the genotype x age-group interaction term.
# Usage: Rscript main_effects.R --config my_config.json [--geno F] [--pheno F] [--out DIR] [--cohort LABEL] [--min_carriers N]
suppressMessages({ library(MASS); library(jsonlite) })

`%||%` <- function(a, b) if (is.null(a)) b else a
.script_dir <- function() {
  if (exists("MAIN_EFFECTS_DIR")) return(MAIN_EFFECTS_DIR)   # set it before source()-ing this file from another script
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(gsub("~+~", " ", f[1], fixed = TRUE))) else getwd()
}
HERE <- .script_dir()
source(file.path(HERE, "00_scoring.R"))   # child LMS scorers (reads growth_refs_lms.csv from HERE)

AGES <- c("all", "lt18", "ge18"); SEXES <- c("all", "male", "female"); TRAITS <- c("height", "bmi", "hc"); ALLELES <- c(DEL = 1L, DUP = 3L)
DEFAULTS <- list(
  cohort_label = "cohort", genotype_file = NULL, phenotype_file = NULL, out_dir = "main_effects_out",
  geno = list(id = "identifier", loci = NULL, include_oneway = FALSE),
  pheno = list(id = "id", sex = "sex", age = "age", age_unit = "years", height_cm = "height_cm", weight_kg = "weight_kg",
               bmi = "bmi", hc_cm = "hc_cm", pcs = "auto", site = NULL, cohort = NULL,
               male = c("male", "Male", "MALE", "M", "m", "1"), female = c("female", "Female", "FEMALE", "F", "f", "2")),
  min_carriers = 5, min_controls = 20, min_adult_ref = 50, adult_age = 18, z_max_child = 5,
  visit_rule = "latest", adult_boxcox = "fixed_ukbb", adjust_age_sex = TRUE, reference_based = TRUE, complementary = "mask_count",
  plausible = list(height_cm = c(30, 250), weight_kg = c(1, 400), bmi = c(5, 100), hc_cm = c(20, 80)))
CHILD_SCORING <- c(height = "WHO 2006/2007 growth reference LMS Z", bmi = "WHO 2006/2007 growth reference LMS Z", hc = "UK 1990 growth reference LMS Z")

read_table <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(as.data.frame(readRDS(path)))
  if (requireNamespace("data.table", quietly = TRUE)) return(as.data.frame(data.table::fread(path, data.table = FALSE, check.names = FALSE, na.strings = c("", "NA", "NaN", "nan", "None"))))
  sep <- if (grepl("\\.csv(\\.gz)?$", path, ignore.case = TRUE)) "," else "\t"
  read.table(path, sep = sep, header = TRUE, check.names = FALSE, quote = "\"", comment.char = "", na.strings = c("", "NA", "NaN", "nan", "None"), stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------------------------- inputs
load_geno <- function(g, cfg) {
  idc <- cfg$geno$id; if (!idc %in% names(g)) stop("genotype matrix: id column '", idc, "' not found")
  g <- g[!is.na(g[[idc]]), , drop = FALSE]; ids <- as.character(g[[idc]])
  if (anyDuplicated(ids)) stop("genotype matrix: duplicated sample ids (one row per sample is required)")
  loci <- cfg$geno$loci
  skipped <- character(0)
  if (is.null(loci)) {   # auto-detect: cytoband-like column names whose values are all 1 / 2 / 3 / missing
    cand <- grep("^([0-9]{1,2}|X|Y)[pq][0-9]", names(g), value = TRUE)
    if (!isTRUE(cfg$geno$include_oneway)) cand <- cand[!grepl("_oneway$", cand)]
    ok <- vapply(cand, function(l) { v <- suppressWarnings(as.numeric(g[[l]])); all(is.na(g[[l]]) | v %in% c(1, 2, 3)) }, logical(1))
    loci <- cand[ok]; skipped <- cand[!ok]
  } else { miss <- setdiff(loci, names(g)); if (length(miss)) stop("genotype matrix: loci not found: ", paste(miss, collapse = ", ")) }
  if (!length(loci)) stop("genotype matrix: no locus columns found (set geno.loci in the config)")
  G <- vapply(loci, function(l) as.integer(suppressWarnings(as.numeric(g[[l]]))), integer(nrow(g)))
  G <- matrix(G, nrow = nrow(g), dimnames = list(ids, loci))
  if (any(!is.na(G) & !G %in% 1:3)) stop("genotype matrix: values other than 1, 2, 3 in the locus columns")
  attr(G, "skipped") <- skipped; G
}

load_pheno <- function(p, cfg) {
  pc <- cfg$pheno; need <- c(pc$id, pc$sex, pc$age); miss <- setdiff(need, names(p))
  if (length(miss)) stop("phenotype table: required column(s) not found: ", paste(miss, collapse = ", "))
  col <- function(nm) if (!is.null(nm) && nm %in% names(p)) suppressWarnings(as.numeric(p[[nm]])) else rep(NA_real_, nrow(p))
  sx <- as.character(p[[pc$sex]]); sex <- ifelse(sx %in% pc$male, "male", ifelse(sx %in% pc$female, "female", NA))
  age <- suppressWarnings(as.numeric(p[[pc$age]])) / switch(pc$age_unit, years = 1, months = 12, days = 365.25, stop("pheno.age_unit must be years, months or days"))
  d <- data.frame(id = as.character(p[[pc$id]]), sex = sex, age = age, height = col(pc$height_cm), weight = col(pc$weight_kg), bmi = col(pc$bmi), hc = col(pc$hc_cm), stringsAsFactors = FALSE)
  for (v in c("height", "weight", "hc")) { m <- median(d[[v]], na.rm = TRUE)
    lo <- c(height = 60, weight = 3, hc = 25)[[v]]; hi <- c(height = 220, weight = 250, hc = 70)[[v]]
    if (is.finite(m) && (m < lo * 0.5 || m > hi)) stop(sprintf("phenotype table: median %s = %.3g; expected cm (height, head circumference) and kg (weight)", v, m)) }
  rng <- cfg$plausible; n_impl <- c(height = 0, weight = 0, bmi = 0, hc = 0)
  for (v in names(n_impl)) { r <- rng[[c(height = "height_cm", weight = "weight_kg", bmi = "bmi", hc = "hc_cm")[[v]]]]
    bad <- !is.na(d[[v]]) & (d[[v]] < r[1] | d[[v]] > r[2]); n_impl[v] <- sum(bad); d[[v]][bad] <- NA }
  calc <- d$weight / (d$height / 100)^2; d$bmi <- ifelse(is.na(d$bmi), calc, d$bmi)
  bad <- !is.na(d$bmi) & (d$bmi < rng$bmi[1] | d$bmi > rng$bmi[2]); n_impl["bmi"] <- n_impl["bmi"] + sum(bad); d$bmi[bad] <- NA
  pcs <- pc$pcs; if (identical(pcs, "auto")) pcs <- head(grep("^PC[0-9]+$", names(p), value = TRUE), 10)
  pcs <- intersect(pcs %||% character(0), names(p))
  for (k in pcs) d[[k]] <- suppressWarnings(as.numeric(p[[k]]))
  d$site <- if (!is.null(pc$site) && pc$site %in% names(p)) as.character(p[[pc$site]]) else NA_character_
  d$cohort <- if (!is.null(pc$cohort) && pc$cohort %in% names(p)) as.character(p[[pc$cohort]]) else cfg$cohort_label
  attr(d, "pcs") <- pcs; attr(d, "n_implausible") <- n_impl; attr(d, "has_site") <- !all(is.na(d$site)); d
}

# one visit per person: among rows flagged keep, the latest (or earliest) age; returns row indices, one per id
pick_visit <- function(id, age, keep, rule = "latest") {
  i <- which(keep & !is.na(age)); if (!length(i)) return(integer(0))
  o <- i[order(id[i], if (rule == "latest") -age[i] else age[i])]; o[!duplicated(id[o])]
}

# adults: the convention of the companion BMI paper (cnv-growth-bmi 00_scoring.R, score_adult_bc(), called here so that
# the two repositories agree): within sex, Box-Cox transform, then lm(transformed y ~ age + age^2) in the cohort's
# non-carrier adults, Z = residual / non-carrier residual SD (|Z| > 10 set to missing by that function).
# adult_boxcox = "fixed_ukbb" (default): exponent fixed at the value estimated once in UK Biobank non-carriers
#   (ADULT_LAMBDA_UKBB); head circumference is not transformed (exponent 1; UK Biobank has no head circumference).
# adult_boxcox = "estimate": exponent estimated in this cohort's non-carrier adults (profile likelihood, grid -2 to 2);
#   only sensible for biobank-size cohorts; a warning is given when the estimate sits on the edge of the grid.
ADULT_LAMBDA_UKBB <- list(height = c(male = 0.65, female = 0.60), bmi = c(male = -0.75, female = -1.05), hc = c(male = 1, female = 1))
score_adult <- function(y, age, sex, is_ref, trait, mode = "fixed_ukbb", min_ref = 50, min_age = 18, label = "") {
  z <- rep(NA_real_, length(y)); sx <- c("male", "female"); ok <- !is.na(y) & !is.na(age) & !is.na(sex) & y > 0 & age >= min_age
  info <- data.frame(sex = sx, n_ref = vapply(sx, function(s) sum(ok & is_ref & sex == s), numeric(1)), n_scored = 0, lambda = NA_real_, lambda_source = NA_character_, stringsAsFactors = FALSE)
  for (s in sx) { i <- which(ok & sex == s); if (!length(i) || info$n_ref[info$sex == s] < min_ref) next
    fixed <- mode == "fixed_ukbb" || trait == "hc"
    zz <- score_adult_bc(y[i], age[i], sex[i], rep("cohort", length(i)), is_ref[i], min_age = min_age, min_ref = min_ref,
                         lambda_from = "(none: estimate in this cohort)", lambda_fixed = if (fixed) ADULT_LAMBDA_UKBB[[trait]][[s]] else NULL)
    a <- attr(zz, "info"); if (is.null(a)) next
    z[i] <- as.numeric(zz); info$n_scored[info$sex == s] <- sum(!is.na(zz)); info$lambda[info$sex == s] <- a$lambda[1]
    info$lambda_source[info$sex == s] <- if (trait == "hc") "none (not transformed)" else if (fixed) "fixed, UK Biobank non-carriers" else "estimated in this cohort"
    if (!fixed && abs(abs(a$lambda[1]) - 2) < 1e-8) warning(sprintf("%s %s %s: the estimated Box-Cox exponent (%.2f) is on the edge of the search grid; use adult_boxcox = 'fixed_ukbb'", label, trait, s, a$lambda[1]), call. = FALSE)
  }
  attr(z, "info") <- info; z
}
adult_scoring_label <- function(trait, mode) sprintf("cohort-relative Z: within sex, %s, then residual of lm(y ~ age + age^2) / residual SD, both from the cohort's non-carrier adults",
  if (trait == "hc") "no transformation" else if (mode == "fixed_ukbb") sprintf("Box-Cox with fixed exponent (UK Biobank non-carriers: men %s, women %s)", ADULT_LAMBDA_UKBB[[trait]][["male"]], ADULT_LAMBDA_UKBB[[trait]][["female"]]) else "Box-Cox with the exponent estimated in this cohort")

# OLS with rank handling; returns coefficient, SE and residual df (checked against lm() in tests/)
ols <- function(X, y) {
  q <- qr(X); r <- q$rank; piv <- q$pivot[seq_len(r)]; b <- qr.coef(q, y); df <- length(y) - r
  res <- y - X[, piv, drop = FALSE] %*% b[piv]; s2 <- sum(res^2) / df
  se <- rep(NA_real_, ncol(X)); se[piv] <- sqrt(diag(chol2inv(qr.R(q)[seq_len(r), seq_len(r), drop = FALSE])) * s2)
  names(se) <- names(b) <- colnames(X); list(b = b, se = se, df = df)
}

# small-cell closure on the 3 x 3 grid (age stratum x sex) of carrier counts for one locus x allele x trait.
# primary: 0 < n < min. complementary: any line (two parts and their total) with exactly one hidden cell gets a second one.
close33 <- function(n, minc) {
  P <- n > 0 & n < minc; S <- P
  lines <- c(lapply(1:3, function(a) cbind(a, c(2, 3, 1))), lapply(1:3, function(s) cbind(c(2, 3, 1), s)))   # parts first, total last
  repeat { changed <- FALSE
    for (ln in lines) { s <- S[ln]; if (sum(s) != 1) next
      cand <- which(!s & n[ln] > 0); if (!length(cand)) next
      S[ln[cand[which.min(n[ln][cand])], , drop = FALSE]] <- TRUE; changed <- TRUE }
    if (!changed) break }
  list(primary = P, complementary = S & !P)
}

fmt_n <- function(n, minc) ifelse(is.na(n), NA_character_, ifelse(n > 0 & n < minc, paste0("<", minc), as.character(n)))

# ---------------------------------------------------------------------------------------------- one cohort
run_cohort <- function(G, ph, cfg, label) {
  minc <- cfg$min_carriers; pcs <- attr(ph, "pcs"); has_site <- attr(ph, "has_site"); log <- list(cohort = label)
  ph <- ph[ph$id %in% rownames(G), , drop = FALSE]
  log$n_pheno_rows <- nrow(ph); bad <- is.na(ph$sex) | is.na(ph$age); log$n_rows_missing_sex_or_age <- sum(bad); ph <- ph[!bad, , drop = FALSE]
  sx <- tapply(ph$sex, ph$id, function(v) length(unique(v))); incons <- names(sx)[sx > 1]; log$n_persons_inconsistent_sex <- length(incons)
  ph <- ph[!ph$id %in% incons, , drop = FALSE]
  if (!nrow(ph)) stop("cohort ", label, ": no phenotype rows match the genotype matrix ids")
  first <- !duplicated(ph$id); P <- ph[first, c("id", "sex", pcs, "site"), drop = FALSE]
  cc <- (if (length(pcs)) complete.cases(P[, pcs, drop = FALSE]) else rep(TRUE, nrow(P))) & (!has_site | !is.na(P$site)); log$n_persons_missing_covariates <- sum(!cc); P <- P[cc, , drop = FALSE]
  ph <- ph[ph$id %in% P$id, , drop = FALSE]; Gc <- G[P$id, , drop = FALSE]; loci <- colnames(Gc)
  P$control <- rowSums(is.na(Gc) | Gc != 2L) == 0   # non-carrier = no CNV at any analysed locus
  log$n_persons <- nrow(P); log$n_noncarriers <- sum(P$control); log$n_loci <- length(loci)
  C0 <- if (length(pcs)) as.matrix(P[, pcs, drop = FALSE]) else NULL
  if (has_site && length(unique(P$site)) > 1) { S <- model.matrix(~ factor(P$site))[, -1, drop = FALSE]; colnames(S) <- paste0("site", seq_len(ncol(S))); C0 <- cbind(C0, S) }
  male <- as.numeric(P$sex == "male"); pos <- match(ph$id, P$id); child <- ph$age < cfg$adult_age
  rows <- list(); arow <- list(); sinfo <- list(); cnt <- list(); k <- 0; ka <- 0
  add <- function(...) { k <<- k + 1; rows[[k]] <<- list(...) }
  for (tr in TRAITS) {
    y <- ph[[tr]]; provided <- any(!is.na(y))
    zc <- rep(NA_real_, nrow(ph))
    if (provided) { zc[child] <- switch(tr, height = score_height(y[child], ph$age[child], ph$sex[child]), bmi = score_bmi_child(y[child], ph$age[child], ph$sex[child]),
                                        hc = score_hc(y[child], ph$age[child], ph$sex[child]))
      log[[paste0("n_child_z_removed_", tr)]] <- sum(!is.na(zc) & abs(zc) > cfg$z_max_child); zc[!is.na(zc) & abs(zc) > cfg$z_max_child] <- NA }
    ic <- pick_visit(ph$id, ph$age, child & !is.na(zc), cfg$visit_rule); ia <- pick_visit(ph$id, ph$age, !child & !is.na(y), cfg$visit_rule)
    z_lt <- age_lt <- z_ge <- age_ge <- y_ge <- rep(NA_real_, nrow(P))
    z_lt[pos[ic]] <- zc[ic]; age_lt[pos[ic]] <- ph$age[ic]; y_ge[pos[ia]] <- y[ia]; age_ge[pos[ia]] <- ph$age[ia]
    za <- score_adult(y_ge, age_ge, P$sex, P$control, tr, cfg$adult_boxcox, cfg$min_adult_ref, cfg$adult_age, label); ai <- attr(za, "info"); z_ge <- as.numeric(za)
    ADULT_SCORING <- adult_scoring_label(tr, cfg$adult_boxcox)
    for (s in c("male", "female")) {
      w <- P$control & P$sex == s & !is.na(z_lt)
      sinfo[[length(sinfo) + 1]] <- data.frame(cohort = label, trait = tr, age_group = "lt18", sex = s, scoring = CHILD_SCORING[[tr]], n_noncarriers_scored = fmt_n(sum(w), minc), lambda = NA, lambda_source = NA,
        mean_z_noncarriers = if (sum(w) >= minc) mean(z_lt[w]) else NA, sd_z_noncarriers = if (sum(w) >= max(minc, 2)) sd(z_lt[w]) else NA, note = "", stringsAsFactors = FALSE)
      a1 <- ai[ai$sex == s, ]; scored <- !is.na(a1$lambda)
      sinfo[[length(sinfo) + 1]] <- data.frame(cohort = label, trait = tr, age_group = "ge18", sex = s, scoring = if (scored) ADULT_SCORING else "none", n_noncarriers_scored = fmt_n(a1$n_ref, minc), lambda = a1$lambda, lambda_source = a1$lambda_source,
        mean_z_noncarriers = if (scored) 0 else NA, sd_z_noncarriers = if (scored) 1 else NA,
        note = if (scored) "" else sprintf("adults not scored: fewer than %d non-carrier adults of this sex with the trait", cfg$min_adult_ref), stringsAsFactors = FALSE)
    }
    adult_first <- cfg$visit_rule == "latest"; use_ge <- if (adult_first) !is.na(z_ge) else is.na(z_lt) & !is.na(z_ge)
    Z <- list(lt18 = z_lt, ge18 = z_ge, all = ifelse(use_ge, z_ge, z_lt)); AGE <- list(lt18 = age_lt, ge18 = age_ge, all = ifelse(use_ge, age_ge, age_lt))
    grp_child <- as.numeric(!use_ge)   # in the all-ages set: 1 = the person's visit is < 18 y
    no_adult <- if (tr == "hc") sprintf("not estimable: there is no adult head-circumference reference and this cohort has fewer than %d non-carrier adults per sex", cfg$min_adult_ref) else
      sprintf("not estimable: adults are scored relative to the cohort and it has fewer than %d non-carrier adults per sex", cfg$min_adult_ref)
    for (a in AGES) for (s in SEXES) {
      insex <- if (s == "all") rep(TRUE, nrow(P)) else P$sex == s
      z <- Z[[a]]; age <- AGE[[a]]; base <- insex & !is.na(z)
      raw_only <- a == "ge18" && !any(base)   # adults measured but not scorable
      if (raw_only) { base <- insex & !is.na(y_ge); age <- age_ge }
      scoring <- if (!provided) "none" else if (raw_only) "none" else switch(a, lt18 = CHILD_SCORING[[tr]], ge18 = ADULT_SCORING,
        all = paste(c(if (any(base & grp_child == 1)) paste("< 18 y:", CHILD_SCORING[[tr]]), if (any(base & grp_child == 0)) paste(">= 18 y:", ADULT_SCORING)), collapse = "; "))
      Cx <- C0
      if (isTRUE(cfg$adjust_age_sex)) { Cx <- cbind(Cx, age = age); if (s == "all") Cx <- cbind(Cx, male = male) }
      if (a == "all") Cx <- cbind(Cx, child = grp_child)
      covtxt <- paste(c(if (isTRUE(cfg$adjust_age_sex)) c("age", if (s == "all") "sex"), pcs, if (!is.null(C0) && any(grepl("^site", colnames(C0)))) "site", if (a == "all") "age group"), collapse = " + ")
      for (L in loci) {
        g <- Gc[, L]; u <- base & !is.na(g) & (g != 2L | P$control); n0 <- sum(u & g == 2L)
        nA <- c(DEL = sum(u & g == 1L), DUP = sum(u & g == 3L)); mage <- c(DEL = mean(age[u & g == 1L]), DUP = mean(age[u & g == 3L]))
        cnt[[paste(tr, L, a, s, sep = "\r")]] <- nA
        st <- NULL; fit <- NULL
        if (!provided) st <- "not estimable: trait not provided" else if (raw_only && any(nA > 0)) st <- no_adult
        else if (n0 < cfg$min_controls) st <- sprintf("not estimable: fewer than %d non-carriers in this cohort and stratum", cfg$min_controls)
        else if (max(nA) >= minc) {
          X <- cbind(`(Intercept)` = 1, DEL = as.numeric(g == 1L), DUP = as.numeric(g == 3L), Cx)[u, , drop = FALSE]
          X <- X[, c(TRUE, nA > 0, rep(TRUE, ncol(X) - 3)), drop = FALSE]; fit <- ols(X, z[u])
        }
        for (al in names(ALLELES)) {
          b <- se <- NA_real_; df <- NA; s1 <- st %||% (if (nA[[al]] == 0) "no carriers" else if (nA[[al]] < minc) "suppressed" else "ok")
          if (is.null(st) && !is.null(fit) && nA[[al]] >= minc) { b <- fit$b[[al]]; se <- fit$se[[al]]; df <- fit$df; if (is.na(b)) s1 <- "not estimable: carrier term collinear with covariates" }
          add(cohort = label, locus = L, allele = al, trait = tr, stratum_age = a, stratum_sex = s, estimator = "within-cohort", estimator_role = "primary",
              n_carriers = nA[[al]], n_noncarriers = n0, beta_sd = b, se = se, df = df, mean_age_carriers = mage[[al]], scoring = scoring, covariates = if (is.na(b)) "" else covtxt, status = s1)
          if (isTRUE(cfg$reference_based) && a == "lt18") {
            v <- z[u & g == ALLELES[[al]]]; okr <- provided && length(v) >= max(minc, 2)
            add(cohort = label, locus = L, allele = al, trait = tr, stratum_age = a, stratum_sex = s, estimator = "reference-based", estimator_role = if (n0 < cfg$min_controls) "fallback (no non-carriers)" else "sensitivity",
                n_carriers = nA[[al]], n_noncarriers = NA, beta_sd = if (okr) mean(v) else NA, se = if (okr) sd(v) / sqrt(length(v)) else NA, df = if (okr) length(v) - 1 else NA,
                mean_age_carriers = mage[[al]], scoring = scoring, covariates = "", status = if (!provided) "not estimable: trait not provided" else if (!length(v)) "no carriers" else if (!okr) "suppressed" else "ok")
          }
        }
        # CNV x age group in one model (each person once, the all-ages visit): interaction = effect(< 18) - effect(>= 18)
        if (a == "all" && provided) {
          nb <- sum(insex & !is.na(g) & (g != 2L | P$control) & !is.na(z_lt) & !is.na(z_ge))   # persons contributing to both age strata (Wald test)
          for (al in names(ALLELES)) {
            ca <- u & g == ALLELES[[al]]; nch <- sum(ca & grp_child == 1); nad <- sum(ca & grp_child == 0)
            c0 <- u & g == 2L; okI <- nch >= minc && nad >= minc && sum(c0 & grp_child == 1) >= cfg$min_controls && sum(c0 & grp_child == 0) >= cfg$min_controls
            d <- sd_ <- dfI <- NA_real_
            if (okI) { uu <- u & (g == 2L | g == ALLELES[[al]]); cv <- as.numeric(g == ALLELES[[al]])
              X <- cbind(`(Intercept)` = 1, carrier = cv, inter = cv * grp_child, Cx)[uu, , drop = FALSE]; f <- ols(X, z[uu]); d <- f$b[["inter"]]; sd_ <- f$se[["inter"]]; dfI <- f$df }
            ka <- ka + 1; arow[[ka]] <- list(cohort = label, locus = L, allele = al, trait = tr, stratum_sex = s, test = "interaction term (one model, one visit per person)",
              estimator_lt18 = "within-cohort", estimator_ge18 = "within-cohort", beta_lt18 = NA_real_, se_lt18 = NA_real_, beta_ge18 = NA_real_, se_ge18 = NA_real_,
              diff_lt18_minus_ge18 = d, se_diff = sd_, df = dfI, n_both = sum(ca & !is.na(z_lt) & !is.na(z_ge)), min_n = min(nch, nad),
              status = if (okI) "ok" else "not estimable: too few carriers or non-carriers in one age group")
          }
        }
      }
    }
  }
  R <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
  # ---- small cells: primary and complementary flags per locus x allele x trait
  R$prim <- FALSE; R$comp <- FALSE; key <- paste(R$trait, R$locus, R$stratum_age, R$stratum_sex, sep = "\r")
  for (tr in TRAITS) for (L in loci) for (al in names(ALLELES)) {
    n <- matrix(0, 3, 3, dimnames = list(AGES, SEXES)); for (a in AGES) for (s in SEXES) n[a, s] <- cnt[[paste(tr, L, a, s, sep = "\r")]][[al]]
    if (!any(n > 0 & n < minc)) next
    cl <- close33(n, minc)
    for (a in AGES) for (s in SEXES) if (cl$primary[a, s] || cl$complementary[a, s]) { w <- key == paste(tr, L, a, s, sep = "\r") & R$allele == al
      R$prim[w] <- cl$primary[a, s]; R$comp[w] <- cl$complementary[a, s] }
  }
  R$mean_age_carriers[is.nan(R$mean_age_carriers)] <- NA
  R$t <- R$beta_sd / R$se; R$p <- 2 * pt(-abs(R$t), R$df)
  # ---- age-group Wald tests from the stratum estimates (uses the unsuppressed values, then applies the same flags)
  strict <- identical(cfg$complementary, "suppress"); hide <- R$prim | (strict & R$comp)
  W <- list()
  for (e1 in c("within-cohort", if (isTRUE(cfg$reference_based)) "reference-based")) {
    a1 <- R[R$stratum_age == "lt18" & R$estimator == e1, ]; a2 <- R[R$stratum_age == "ge18" & R$estimator == "within-cohort", ]
    m <- merge(a1, a2, by = c("cohort", "locus", "allele", "trait", "stratum_sex"), suffixes = c("_lt18", "_ge18"))
    d <- m$beta_sd_lt18 - m$beta_sd_ge18; sed <- sqrt(m$se_lt18^2 + m$se_ge18^2); ok <- is.finite(d) & is.finite(sed)
    W[[e1]] <- data.frame(m[, c("cohort", "locus", "allele", "trait", "stratum_sex")], test = "Wald (difference of the two stratum estimates)", estimator_lt18 = e1, estimator_ge18 = "within-cohort",
      beta_lt18 = m$beta_sd_lt18, se_lt18 = m$se_lt18, beta_ge18 = m$beta_sd_ge18, se_ge18 = m$se_ge18, diff_lt18_minus_ge18 = d, se_diff = sed, df = NA_real_,
      n_carriers_lt18 = m$n_carriers_lt18, n_carriers_ge18 = m$n_carriers_ge18, hide = (m$prim_lt18 | m$prim_ge18 | (strict & (m$comp_lt18 | m$comp_ge18))),
      mask_lt18 = m$comp_lt18, mask_ge18 = m$comp_ge18,
      status = ifelse(ok, "ok", "not estimable: one or both stratum estimates missing"), stringsAsFactors = FALSE)
  }
  W <- do.call(rbind, W)
  if (length(arow)) {
    I <- do.call(rbind, lapply(arow, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
    lk <- function(a) { r <- R[R$stratum_age == a & R$estimator == "within-cohort", ]; r[match(paste(I$locus, I$allele, I$trait, I$stratum_sex), paste(r$locus, r$allele, r$trait, r$stratum_sex)), ] }
    r1 <- lk("lt18"); r2 <- lk("ge18")
    I$n_carriers_lt18 <- r1$n_carriers; I$n_carriers_ge18 <- r2$n_carriers; I$hide <- r1$prim | r2$prim | (strict & (r1$comp | r2$comp)) | (I$min_n > 0 & I$min_n < minc)
    I$mask_lt18 <- r1$comp; I$mask_ge18 <- r2$comp
    nb <- I$n_both; names(nb) <- paste(I$locus, I$allele, I$trait, I$stratum_sex); W$n_both <- nb[paste(W$locus, W$allele, W$trait, W$stratum_sex)]
    I$min_n <- NULL; W <- rbind(W, I[, names(W)])
  } else W$n_both <- NA
  W$z_or_t <- W$diff_lt18_minus_ge18 / W$se_diff; W$p_diff <- ifelse(is.na(W$df), 2 * pnorm(-abs(W$z_or_t)), 2 * pt(-abs(W$z_or_t), W$df))
  numA <- c("beta_lt18", "se_lt18", "beta_ge18", "se_ge18", "diff_lt18_minus_ge18", "se_diff", "df", "z_or_t", "p_diff")
  W[W$hide, numA] <- NA; W$status[W$hide] <- "suppressed"
  mk <- function(n, mask) ifelse(mask & n >= minc, paste0(">=", minc), fmt_n(n, minc))
  W$n_carriers_lt18 <- mk(W$n_carriers_lt18, W$mask_lt18); W$n_carriers_ge18 <- mk(W$n_carriers_ge18, W$mask_ge18)
  W$n_carriers_in_both_strata <- fmt_n(W$n_both, minc)
  W$note <- ifelse(W$estimator_lt18 != W$estimator_ge18, "mixed estimators: reference-based < 18 vs within-cohort >= 18", "")
  W <- W[, c("cohort", "locus", "allele", "trait", "stratum_sex", "test", "estimator_lt18", "estimator_ge18", "n_carriers_lt18", "n_carriers_ge18", "n_carriers_in_both_strata",
             "beta_lt18", "se_lt18", "beta_ge18", "se_ge18", "diff_lt18_minus_ge18", "se_diff", "z_or_t", "df", "p_diff", "status", "note")]
  # ---- apply suppression to the main table
  num <- c("beta_sd", "se", "t", "df", "p", "mean_age_carriers")
  R[hide, num] <- NA; R$covariates[hide] <- ""
  R$status[R$prim] <- sprintf("suppressed: fewer than %d carriers", minc)
  R$status[R$comp & !R$prim] <- if (strict) "suppressed: complementary (protects a small cell in another stratum)" else ifelse(R$status[R$comp & !R$prim] == "ok", "ok (carrier count masked to protect a small cell in another stratum)", R$status[R$comp & !R$prim])
  R$status[R$status == "suppressed"] <- sprintf("suppressed: fewer than %d carriers", minc)
  R$n_carriers <- ifelse(R$comp & !R$prim, paste0(">=", minc), fmt_n(R$n_carriers, minc)); R$n_noncarriers <- fmt_n(R$n_noncarriers, minc)
  R$min_carriers <- minc
  R <- R[, c("cohort", "locus", "allele", "trait", "stratum_age", "stratum_sex", "estimator", "estimator_role", "n_carriers", "n_noncarriers", "beta_sd", "se", "t", "df", "p",
             "mean_age_carriers", "scoring", "covariates", "status", "min_carriers")]
  list(main = R, age = W, scoring = do.call(rbind, sinfo), log = log)
}

# ---------------------------------------------------------------------------------------------- driver
run_main_effects <- function(cfg = list(), geno = NULL, pheno = NULL, write = TRUE) {
  cfg <- modifyList(DEFAULTS, cfg)
  if (!cfg$visit_rule %in% c("latest", "earliest")) stop("visit_rule must be 'latest' or 'earliest'")
  if (!cfg$adult_boxcox %in% c("fixed_ukbb", "estimate")) stop("adult_boxcox must be 'fixed_ukbb' or 'estimate'")
  if (!cfg$complementary %in% c("mask_count", "suppress")) stop("complementary must be 'mask_count' or 'suppress'")
  geno <- geno %||% read_table(cfg$genotype_file); pheno <- pheno %||% read_table(cfg$phenotype_file)
  G <- load_geno(geno, cfg); ph <- load_pheno(pheno, cfg)
  out <- list()
  for (co in unique(ph$cohort)) { sub <- ph[ph$cohort == co, , drop = FALSE]; attr(sub, "pcs") <- attr(ph, "pcs"); attr(sub, "has_site") <- attr(ph, "has_site") && !all(is.na(sub$site))
    message("cohort ", co, ": ", length(unique(sub$id)), " persons with phenotypes, ", ncol(G), " loci"); out[[co]] <- run_cohort(G, sub, cfg, co) }
  res <- list(main = do.call(rbind, lapply(out, `[[`, "main")), age = do.call(rbind, lapply(out, `[[`, "age")), scoring = do.call(rbind, lapply(out, `[[`, "scoring")))
  res$log <- list(module = "main_effects.R 2026-11-20", R = R.version.string, date = as.character(Sys.Date()),
                  settings = cfg[c("min_carriers", "min_controls", "min_adult_ref", "adult_age", "z_max_child", "visit_rule", "adult_boxcox", "adjust_age_sex", "reference_based", "complementary", "plausible")],
                  covariate_pcs = attr(ph, "pcs"), loci = colnames(G), locus_like_columns_skipped = attr(G, "skipped"),
                  n_values_outside_plausible_range = as.list(attr(ph, "n_implausible")), cohorts = unname(lapply(out, `[[`, "log")))
  if (write) {
    dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE); lab <- gsub("[^A-Za-z0-9_.-]", "_", cfg$cohort_label)
    wr <- function(d, f) write.table(d, file.path(cfg$out_dir, f), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    wr(res$main, sprintf("main_effects_%s.tsv", lab)); wr(res$age, sprintf("age_group_test_%s.tsv", lab)); wr(res$scoring, sprintf("scoring_info_%s.tsv", lab))
    writeLines(toJSON(res$log, auto_unbox = TRUE, pretty = TRUE, null = "null"), file.path(cfg$out_dir, sprintf("run_log_%s.json", lab)))
    message("wrote aggregate tables to ", normalizePath(cfg$out_dir))
  }
  invisible(res)
}

if (!exists("MAIN_EFFECTS_SOURCE_ONLY") && !interactive()) {
  a <- commandArgs(TRUE); opt <- function(k) { i <- match(k, a); if (is.na(i)) NULL else a[i + 1] }
  if (is.null(opt("--config")) && is.null(opt("--geno"))) stop("usage: Rscript main_effects.R --config my_config.json [--geno F] [--pheno F] [--out DIR] [--cohort LABEL] [--min_carriers N]")
  cfg <- if (!is.null(opt("--config"))) fromJSON(opt("--config"), simplifyVector = TRUE) else list()
  cfg <- cfg[!grepl("^_", names(cfg))]
  if (!is.null(opt("--geno"))) cfg$genotype_file <- opt("--geno")
  if (!is.null(opt("--pheno"))) cfg$phenotype_file <- opt("--pheno")
  if (!is.null(opt("--out"))) cfg$out_dir <- opt("--out")
  if (!is.null(opt("--cohort"))) cfg$cohort_label <- opt("--cohort")
  if (!is.null(opt("--min_carriers"))) cfg$min_carriers <- as.numeric(opt("--min_carriers"))
  run_main_effects(cfg)
}
