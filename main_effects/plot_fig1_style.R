# plot_fig1_style.R  (2026-11-20)
# What: two figures from the aggregate tables of this module, in the colours of fig1_main_effects (deletion red,
#   duplication blue).
#   (1) fig1_style_main_effects.png: deletion and duplication effect per locus with 95% CI for height, BMI and head
#       circumference. Head circumference is marked with an asterisk and a footnote ("paediatric clinical cohorts")
#       and is drawn from the < 18 y stratum by default; height and BMI from the all-ages stratum.
#   (2) fig1_style_age_groups.png: effect < 18 y and >= 18 y per locus as paired points joined by a line, with the P
#       value of the difference.
#   Filled points = within-cohort estimates; open points = reference-based or mixed estimates (flagged in the caption).
# Inputs: --main meta_main_effects.tsv (or one cohort's main_effects_<cohort>.tsv)
#         --age  meta_age_group_test.tsv (or age_group_test_<cohort>.tsv); optional
#         --out  output folder; optional --stratum_age all --hc_stratum lt18 --sex all --width 11
# Outputs: the two PNG files and fig1_style_plotted_values.tsv (the numbers drawn)
# null model: none; this script only draws estimates (95% CI = estimate +- 1.96 SE) and prints the P values computed
#   by main_effects.R / meta_main_effects.R (Wald test of the age-group difference).
suppressMessages(library(ggplot2))
COLS <- c(DEL = "#d62728", DUP = "#1f4e9c"); TLAB <- c(height = "Height", bmi = "BMI", hc = "Head circumference*")
EST_ORDER <- c("within-cohort", "mixed (within-cohort + reference-based)", "reference-based")
FOOT <- "* Head circumference: paediatric clinical cohorts (age < 18 y; Z scores from the UK 1990 growth reference)."

read_est <- function(f) {
  d <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE, na.strings = "NA")
  if ("beta_fe" %in% names(d)) { d$beta <- d$beta_fe; d$se <- d$se_fe } else d$beta <- d$beta_sd
  d$beta <- suppressWarnings(as.numeric(d$beta)); d$se <- suppressWarnings(as.numeric(d$se)); d[is.finite(d$beta) & is.finite(d$se), ]
}
pretty_locus <- function(x) { s <- gsub("_", " ", x); s <- gsub("(BP[0-9]) (BP[0-9])", "\\1-\\2", s); gsub(" ([A-H]) ([A-H])$", " \\1-\\2", s) }
locus_order <- function(x) { u <- unique(x); chr <- suppressWarnings(as.numeric(sub("^([0-9]+).*", "\\1", u))); arm <- sub("^[0-9XY]+([pq]).*", "\\1", u); u[order(chr, arm, u)] }
best <- function(d, by) {   # one estimator per cell: within-cohort, else mixed, else reference-based
  d$rank <- match(d$estimator, EST_ORDER); d <- d[order(d$rank), ]; d[!duplicated(d[, by]), ] }
theme_f1 <- function() theme_classic(base_size = 10) + theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.4), strip.background = element_blank(),
  strip.text = element_text(size = 10.5), legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 8), panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3))

plot_main <- function(D, stratum_age = "all", hc_stratum = "lt18", sex = "all") {
  d <- D[D$stratum_sex == sex & ((D$trait != "hc" & D$stratum_age == stratum_age) | (D$trait == "hc" & D$stratum_age == hc_stratum)), ]
  d <- best(d, c("locus", "allele", "trait")); if (!nrow(d)) stop("nothing to plot")
  lv <- rev(locus_order(d$locus)); d$locus_f <- factor(d$locus, levels = lv, labels = pretty_locus(lv))
  d$trait_f <- factor(TLAB[d$trait], levels = TLAB); d$within <- ifelse(d$estimator == "within-cohort", "within-cohort contrast", "reference-based or mixed")
  d$within <- factor(d$within, levels = c("within-cohort contrast", "reference-based or mixed"))
  pd <- position_dodge(width = 0.55)
  p <- ggplot(d, aes(beta, locus_f, colour = allele, shape = within)) + geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_linerange(aes(xmin = beta - 1.96 * se, xmax = beta + 1.96 * se), linewidth = 0.4, position = pd) + geom_point(size = 1.9, position = pd, fill = "white") +
    facet_wrap(~trait_f, nrow = 1, scales = "free_x") + scale_colour_manual(values = COLS, name = NULL, labels = c(DEL = "Deletion", DUP = "Duplication")) +
    scale_shape_manual(values = c("within-cohort contrast" = 16, "reference-based or mixed" = 21), name = NULL, drop = TRUE) +
    labs(x = "Effect size (s.d.) with 95% CI", y = NULL, caption = paste0(FOOT, "\nHeight and BMI: ", c(all = "all ages", lt18 = "age < 18 y", ge18 = "age >= 18 y")[[stratum_age]],
         ". Open symbols: carriers compared with the growth reference (no within-cohort non-carriers).")) + theme_f1()
  list(plot = p, data = d, n_loci = length(lv))
}

plot_age <- function(D, A, sex = "all") {
  d <- D[D$stratum_sex == sex & D$stratum_age %in% c("lt18", "ge18"), ]; d <- best(d, c("locus", "allele", "trait", "stratum_age"))
  both <- aggregate(stratum_age ~ locus + allele + trait, d, function(v) length(unique(v))); both <- both[both$stratum_age == 2, c("locus", "allele", "trait")]
  d <- merge(d, both); if (!nrow(d)) return(NULL)
  lv <- rev(locus_order(d$locus)); d$locus_f <- factor(d$locus, levels = lv, labels = pretty_locus(lv)); d$trait_f <- factor(TLAB[d$trait], levels = TLAB)
  d$age_f <- factor(ifelse(d$stratum_age == "lt18", "< 18 y", ">= 18 y"), levels = c("< 18 y", ">= 18 y"))
  d$allele_f <- factor(c(DEL = "Deletion", DUP = "Duplication")[d$allele], levels = c("Deletion", "Duplication"))
  w <- reshape(d[, c("locus_f", "allele", "allele_f", "trait_f", "stratum_age", "beta")], idvar = c("locus_f", "allele", "allele_f", "trait_f"), timevar = "stratum_age", direction = "wide")
  lab <- NULL
  if (!is.null(A)) { a <- A[A$stratum_sex == sex & grepl("^Wald", A$test), ]; a$p_diff <- suppressWarnings(as.numeric(a$p_diff)); a <- a[is.finite(a$p_diff), ]
    a$rank <- match(a$estimator_lt18, EST_ORDER); a <- a[order(a$rank), ]; a <- a[!duplicated(a[, c("locus", "allele", "trait")]), ]
    lab <- merge(unique(d[, c("locus", "locus_f", "allele", "allele_f", "trait", "trait_f")]), a[, c("locus", "allele", "trait", "p_diff", "estimator_lt18")])
    if (nrow(lab)) { lab$txt <- paste0("P = ", formatC(lab$p_diff, format = "g", digits = 2), ifelse(lab$estimator_lt18 != "within-cohort", "†", ""))
      xr <- aggregate(cbind(hi = beta + 1.96 * se) ~ trait_f, d, max); lab <- merge(lab, xr) } else lab <- NULL }
  p <- ggplot(d, aes(beta, locus_f, colour = allele)) + geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_segment(data = w, aes(x = beta.lt18, xend = beta.ge18, y = locus_f, yend = locus_f), linewidth = 0.5, alpha = 0.6) +
    geom_linerange(aes(xmin = beta - 1.96 * se, xmax = beta + 1.96 * se), linewidth = 0.3, alpha = 0.45) +
    geom_point(aes(shape = age_f), size = 2.1, fill = "white") + facet_grid(allele_f ~ trait_f, scales = "free_x") +
    scale_colour_manual(values = COLS, guide = "none") + scale_shape_manual(values = c("< 18 y" = 21, ">= 18 y" = 16), name = "Age at measurement") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.32))) +
    labs(x = "Effect size (s.d.) with 95% CI", y = NULL, caption = paste0("P: Wald test of the difference between the < 18 y and >= 18 y effects. ", FOOT,
         "\n† < 18 y effect from a reference-based or mixed estimate. Adult head circumference is cohort-relative (no adult reference).")) + theme_f1()
  if (!is.null(lab)) p <- p + geom_text(data = lab, aes(x = hi, y = locus_f, label = txt), inherit.aes = FALSE, hjust = -0.12, size = 2.5, colour = "grey25")
  list(plot = p, data = d, n_loci = length(lv))
}

run_plots <- function(main, out_dir, age = NULL, stratum_age = "all", hc_stratum = "lt18", sex = "all", width = 11) {
  D <- read_est(main); A <- if (!is.null(age)) read.delim(age, check.names = FALSE, stringsAsFactors = FALSE, na.strings = "NA") else NULL
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE); keep <- c("locus", "allele", "trait", "stratum_age", "stratum_sex", "estimator", "beta", "se")
  m <- plot_main(D, stratum_age, hc_stratum, sex); h <- max(3, 1.6 + 0.2 * m$n_loci)
  ggsave(file.path(out_dir, "fig1_style_main_effects.png"), m$plot, width = width, height = h, dpi = 200, limitsize = FALSE); vals <- cbind(figure = "main_effects", m$data[, keep])
  g <- plot_age(D, A, sex)
  if (!is.null(g)) { ggsave(file.path(out_dir, "fig1_style_age_groups.png"), g$plot, width = width, height = max(4, 2.2 + 0.36 * g$n_loci), dpi = 200, limitsize = FALSE)
    vals <- rbind(vals, cbind(figure = "age_groups", g$data[, keep])) } else message("no locus has both a < 18 y and a >= 18 y estimate: age-group figure not drawn")
  write.table(vals, file.path(out_dir, "fig1_style_plotted_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  message("wrote figures to ", normalizePath(out_dir)); invisible(list(main = m, age = g))
}

if (!exists("MAIN_EFFECTS_SOURCE_ONLY") && !interactive()) {
  a <- commandArgs(TRUE); opt <- function(k, d = NULL) { i <- match(k, a); if (is.na(i)) d else a[i + 1] }
  if (is.null(opt("--main")) || is.null(opt("--out"))) stop("usage: Rscript plot_fig1_style.R --main meta_main_effects.tsv [--age meta_age_group_test.tsv] --out DIR")
  run_plots(opt("--main"), opt("--out"), opt("--age"), opt("--stratum_age", "all"), opt("--hc_stratum", "lt18"), opt("--sex", "all"), as.numeric(opt("--width", "11")))
}
