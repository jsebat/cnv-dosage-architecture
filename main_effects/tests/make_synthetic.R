# make_synthetic.R  (2026-11-20)
# What: SYNTHETIC test cohorts for the main_effects module, with known CNV effects. Nothing here is real: ids are
#   FAKE_<cohort>_<number>, genotypes and measurements are simulated. Four cohorts:
#     A, C  population-like cohorts with children and adults, non-carriers, PC1-PC10
#     K     paediatric clinic with its own non-carriers (same age distribution), few adults, repeated visits
#     B     carrier-only paediatric clinic with repeated visits (reference-based fallback)
# Inputs: none (growth_refs_lms.csv of the module is used to turn simulated Z scores into cm / kg).
# Outputs (in the folder given as argument; default a temporary folder, never the repository): for each cohort a
#   genotype matrix in the format of GenotypeRecurrentCNVs_tutorial.ipynb (genomatrix_<cohort>.csv), a phenotype table
#   (pheno_<cohort>.tsv) and a config (config_<cohort>.json); plus truth.tsv with the simulated effects.
# null model: none; this script only simulates data.
# Usage: Rscript tests/make_synthetic.R [out_dir]      (also sourced by tests/test_main_effects.R)

# simulated effects in Z units (SD); bmi_adult differs from bmi_child at 16p11.2 DEL to test the age-group comparison
TRUTH <- data.frame(
  locus  = c("16p11.2_BP4_BP5", "16p11.2_BP4_BP5", "22q11.21_A_D", "22q11.21_A_D", "1q21.1_BP3_BP4", "1q21.1_BP3_BP4", "15q13.3_BP4_BP5", "15q13.3_BP4_BP5"),
  allele = c("DEL", "DUP", "DEL", "DUP", "DEL", "DUP", "DEL", "DUP"),
  height = c(-0.4, -0.3, -0.8, 0, -0.3, 0, 0, 0), bmi_child = c(0.3, -0.4, 0, 0, 0, 0, 0, 0), bmi_adult = c(0.9, -0.4, 0, 0, 0, 0, 0, 0),
  hc = c(0.6, -0.5, -0.6, 0, -0.8, 0.6, 0, 0), stringsAsFactors = FALSE)
LOCI <- unique(TRUTH$locus)

.inv_lms <- function(z, age, sex, ref, item, LMS) {   # Z -> measurement with the same table and interpolation as the scorer
  y <- rep(NA_real_, length(z))
  for (s in c("male", "female")) { t <- LMS[LMS$ref == ref & LMS$item == item & LMS$sex == s, ]; t <- t[!is.na(t$mu) & !is.na(t$sigma) & !is.na(t$nu), ]; t <- t[order(t$age), ]
    i <- which(sex == s); if (!length(i)) next
    L <- approx(t$age, t$nu, age[i])$y; M <- approx(t$age, t$mu, age[i])$y; S <- approx(t$age, t$sigma, age[i])$y
    y[i] <- ifelse(abs(L) < 1e-8, M * exp(S * z[i]), M * (1 + L * S * z[i])^(1 / L)) }
  y
}

# n_car: named vector "locus:allele" -> carriers; visits: maximum visits per person (1 = cross-sectional)
make_cohort <- function(label, n_ctrl, n_car, frac_adult, seed, LMS, visits = 1, with_pcs = TRUE) {
  set.seed(seed); n <- n_ctrl + sum(n_car)
  id <- sprintf("FAKE_%s_%06d", label, seq_len(n)); sex <- sample(c("male", "female"), n, TRUE)
  G <- matrix(2L, n, length(LOCI), dimnames = list(id, LOCI)); at <- n_ctrl
  for (k in names(n_car)) { p <- strsplit(k, ":")[[1]]; if (n_car[[k]] > 0) G[at + seq_len(n_car[[k]]), p[1]] <- c(DEL = 1L, DUP = 3L)[[p[2]]]; at <- at + n_car[[k]] }
  eff <- function(col) { e <- rep(0, n); for (r in seq_len(nrow(TRUTH))) { w <- G[, TRUTH$locus[r]] == c(DEL = 1L, DUP = 3L)[[TRUTH$allele[r]]]; e[w] <- e[w] + TRUTH[[col]][r] }; e }
  adult <- runif(n) < frac_adult; e_h <- rnorm(n); e_b <- rnorm(n); e_c <- rnorm(n)
  rows <- list()
  for (v in seq_len(visits)) {
    keep <- if (v == 1) rep(TRUE, n) else runif(n) < 0.6
    age <- ifelse(adult, runif(n, 18.5, 70), runif(n, 1 + (v - 1) * 2, min(17.9, 9 + (v - 1) * 2.2)))   # later visits are older, all children stay < 18
    nz <- if (visits > 1) 0.15 else 0   # small visit-to-visit noise
    zh <- (e_h + rnorm(n, 0, nz)) / sqrt(1 + nz^2) + eff("height"); zc <- (e_c + rnorm(n, 0, nz)) / sqrt(1 + nz^2) + eff("hc")
    zb <- (e_b + rnorm(n, 0, nz)) / sqrt(1 + nz^2) + ifelse(adult, eff("bmi_adult"), eff("bmi_child"))
    young <- age < 5
    h <- ifelse(adult, ifelse(sex == "male", 176, 163) - 0.03 * (age - 45) + 6.5 * zh,
                ifelse(young, .inv_lms(zh, age, sex, "who", "height", LMS), .inv_lms(zh, age, sex, "who2007", "height", LMS)))
    b <- ifelse(adult, exp(log(26) + 0.004 * (age - 45) + 0.16 * zb),
                ifelse(young, .inv_lms(zb, age, sex, "who", "bmi", LMS), .inv_lms(zb, age, sex, "who2007", "bmi", LMS)))
    hc <- ifelse(adult, ifelse(sex == "male", 57, 55.3) + 1.6 * zc, .inv_lms(zc, age, sex, "uk1990", "hc", LMS))
    d <- data.frame(subject = id, gender = ifelse(sex == "male", "M", "F"), age_years = round(age, 3), height_cm = round(h, 2), weight_kg = round(b * (h / 100)^2, 3), ofc_cm = round(hc, 2), stringsAsFactors = FALSE)
    rows[[v]] <- d[keep & (v == 1 | !adult), ]
  }
  ph <- do.call(rbind, rows)
  bad <- sample(which(ph$age_years < 17 & ph$age_years > 2), 6); ph$height_cm[bad] <- ph$height_cm[bad] * 1.9   # entry errors: must be removed by |Z| > 5 / plausible range
  if (with_pcs) { set.seed(seed + 1); pc <- matrix(rnorm(n * 10), n, 10, dimnames = list(NULL, paste0("PC", 1:10))); ph <- cbind(ph, pc[match(ph$subject, id), ]) }
  ph <- ph[sample(nrow(ph)), ]
  # genotype matrix as written by the genotyping notebook: index column, identifier, manifest columns, loci, "_oneway" and helper columns
  gm <- data.frame(identifier = id, sampleID = id, QC = "pass", G, check.names = FALSE, stringsAsFactors = FALSE)
  gm[["16p11.2_BP4_BP5_oneway"]] <- G[, "16p11.2_BP4_BP5"]; gm[["16p11.2"]] <- paste0("222", G[, "16p11.2_BP4_BP5"]); gm[["22q11.21_AD"]] <- strrep(G[, "22q11.21_A_D"], 3)
  list(geno = gm, pheno = ph, G = G)
}

write_cohort <- function(x, label, dir, min_carriers = 5) {
  write.csv(x$geno, file.path(dir, sprintf("genomatrix_%s.csv", label)))   # row names = the pandas index column
  write.table(x$pheno, file.path(dir, sprintf("pheno_%s.tsv", label)), sep = "\t", quote = FALSE, row.names = FALSE)
  cfg <- list(cohort_label = label, genotype_file = file.path(dir, sprintf("genomatrix_%s.csv", label)), phenotype_file = file.path(dir, sprintf("pheno_%s.tsv", label)),
              out_dir = file.path(dir, paste0("out_", label)), min_carriers = min_carriers,
              pheno = list(id = "subject", sex = "gender", age = "age_years", height_cm = "height_cm", weight_kg = "weight_kg", hc_cm = "ofc_cm"))
  f <- file.path(dir, sprintf("config_%s.json", label)); writeLines(jsonlite::toJSON(cfg, auto_unbox = TRUE, pretty = TRUE), f); f
}

make_all <- function(dir, lms_file) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE); LMS <- read.csv(lms_file)
  big <- c("16p11.2_BP4_BP5:DEL" = 400, "16p11.2_BP4_BP5:DUP" = 400, "22q11.21_A_D:DEL" = 400, "22q11.21_A_D:DUP" = 300, "1q21.1_BP3_BP4:DEL" = 300, "1q21.1_BP3_BP4:DUP" = 300, "15q13.3_BP4_BP5:DEL" = 3)
  co <- list(A = make_cohort("A", 20000, big, 0.55, 11, LMS), C = make_cohort("C", 12000, round(big * 0.7), 0.5, 23, LMS),
             K = make_cohort("K", 1500, c("16p11.2_BP4_BP5:DEL" = 300, "16p11.2_BP4_BP5:DUP" = 150, "22q11.21_A_D:DEL" = 300, "1q21.1_BP3_BP4:DUP" = 4), 0.03, 37, LMS, visits = 4, with_pcs = FALSE),
             B = make_cohort("B", 0, c("16p11.2_BP4_BP5:DEL" = 250, "22q11.21_A_D:DEL" = 250, "1q21.1_BP3_BP4:DUP" = 3), 0.06, 41, LMS, visits = 4, with_pcs = FALSE))
  cfgs <- vapply(names(co), function(l) write_cohort(co[[l]], l, dir), character(1))
  write.table(TRUTH, file.path(dir, "truth.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  list(cohorts = co, configs = cfgs, truth = TRUTH)
}

if (!exists("MAIN_EFFECTS_SOURCE_ONLY") && !interactive()) {
  a <- commandArgs(FALSE); f <- gsub("~+~", " ", sub("^--file=", "", a[grep("^--file=", a)][1]), fixed = TRUE)
  out <- commandArgs(TRUE)[1]; if (is.na(out)) out <- file.path(tempdir(), "main_effects_synthetic")
  make_all(out, file.path(dirname(dirname(normalizePath(f))), "growth_refs_lms.csv")); message("synthetic cohorts written to ", normalizePath(out))
}
