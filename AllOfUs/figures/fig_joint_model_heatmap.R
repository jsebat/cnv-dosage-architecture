# fig_joint_model_heatmap.R -- one-panel summary of the joint models of the height and BMI CNV dose response (Jon 2026-11-15: "a heat map";
#   "see how it's done in the Shanta paper with asterisks and boxes, and the grid is very small" -> style of Engchuan/conga_v2/code/fig2_mirror_matrix.R:
#   tile = z capped at +-4, * P < 0.05, ** BH q < 0.05 over the tiles of the figure, black box = P < 0.05 and the term holds in EVERY single-locus removal,
#   each removal against its own matched null; small format, 7-pt text).
#   Rows = joint-model terms; columns = level x robustness drop:
#     recurrent CNVs (34 Sacks intervals; bivariate REML joint model, growth_joint_bivariate.R, growth_joint_bivariate_strata.R, growth_mirror_decomp.R):
#       all (single dosage-sensitive gene in the model), without single-gene DUPs (union-gene DUP arms dropped), without single-gene CNVs, without the 13 asymmetric loci;
#     UK Biobank genome-wide (UKBB session; results_aggregate/robust_strata_2026-11-15 S0-S3, union_strata_2026-11-15/joint_noUnionAny):
#       all (union presence as a covariate), without single-gene DUP carriers, without single-gene CNV carriers, without asymmetric-locus carriers, without both.
#   Single-locus removals: 34 intervals (recurrent CNVs), 24 recurrent units (UK Biobank). Mirror rows: z of the deletion-duplication anticorrelation
#   (likelihood ratio; recurrent CNVs, all only; sign set so a mirror is positive). Aggregates only. Output: fig_joint_model_heatmap.{pdf,png}, joint_model_heatmap.csv
suppressMessages({ library(data.table); library(ggplot2) })
# AllOfUs copy (2026-11-16). Locus columns: lab Dropbox (SEBATLAB_DROPBOX). UK Biobank columns: the UKBB aggregates by default, or, with
#   GENOMEWIDE_META=<meta_sacks2026.py output>, the UK Biobank + All of Us meta-analysis (<meta>/joint/<carrier subset>/<trait>/grid_summary.tsv).
#   FIG_OUT = output folder (default: the locus results folder V). Rscript fig_joint_model_heatmap.R
DBX <- Sys.getenv("SEBATLAB_DROPBOX", file.path(Sys.getenv("USERPROFILE", Sys.getenv("HOME")), "sebatlab Dropbox"))
JC <- file.path(DBX, "Jon/Claude"); V <- file.path(JC, "two_axis_CNV_dose_response/v2_2026-11-08/results/eight_traits/growth_conga")
RA <- file.path(JC, "UKBB_anthro_CNV/results_aggregate"); META <- Sys.getenv("GENOMEWIDE_META", ""); FO <- Sys.getenv("FIG_OUT", V)
SUBSET <- c(G_all = "all", G_noUD = "noUnionDUP", G_noU = "noUnionAny", G_noA = "noAsym13", G_noUDA = "noUnionDUP_noAsym13")
# Jon 2026-11-15: rows grouped by the axis each term captures, axis named on the left
TLAB <- c(mirror_tot = "Total", mirror_tw = "Explained by TWAS genes", tw_DEL = "   deletions", tw_DUP = "   duplications", mirror_res = "Not explained by TWAS",
          hit_DEL = "Deletion", hit_DUP = "Duplication", con_DEL = "Deletion", con_DUP = "Duplication", nall_DEL = "Deletion", nall_DUP = "Duplication")
AXIS <- c(mirror_tot = "Mirror", mirror_tw = "Mirror", tw_DEL = "Mirror", tw_DUP = "Mirror", mirror_res = "Mirror",
          hit_DEL = "Single\nhit", hit_DUP = "Single\nhit", con_DEL = "Constraint", con_DUP = "Constraint", nall_DEL = "Mean\nshift", nall_DUP = "Mean\nshift")
AXL <- c("Mirror", "Constraint", "Single\nhit", "Mean\nshift")   # Jon 2026-11-15: axis names and order
CL <- c("L_all", "L_noUD", "L_noU", "L_noA", "G_all", "G_noUD", "G_noU", "G_noA", "G_noUDA")
CLAB <- c(L_all = "All CNVs", L_noUD = "No single-gene DUPs", L_noU = "No single-gene CNVs", L_noA = "No asymmetric loci",
          G_all = "All CNVs", G_noUD = "No single-gene DUPs", G_noU = "No single-gene CNVs", G_noA = "No asymmetric loci", G_noUDA = "Both removed")
out <- list()
add <- function(trait, col, term, z, P, all_held) out[[length(out) + 1]] <<- data.table(trait, col, term, z, P, all_held)
allh <- function(s) { p <- as.numeric(strsplit(s, "/")[[1]]); p[2] > 0 && p[1] == p[2] }
## recurrent CNVs
lmap <- c(mirror = "mirror_tw", tw_DEL = "tw_DEL", tw_DUP = "tw_DUP", hit_DEL = "hit_DEL", hit_DUP = "hit_DUP", con_DEL = "con_DEL", con_DUP = "con_DUP",
          nall_DEL = "nall_DEL", nall_DUP = "nall_DUP")
F0 <- fread(file.path(V, "joint_bivariate_2026-11-15/joint_bivariate_summary.csv"))[grepl("^unsigned", twas) & readout %in% names(lmap)]
# Jon 2026-11-16: the mean shift (total genes) is not rotated by the shift null, and under that null it absorbs the effect of the scrambled
#   annotation terms (UKBB BMI DEL: observed t +4.7, null mean t +7.35 -> z_emp -4.3, the opposite sign of the raw effect, which is upper right as in
#   the burden bins). Mean-shift tiles therefore show the analytic z (t against 0); the annotation terms keep the shift-null z.
nall_z <- function(ro, za, ze) if (grepl("^nall", ro)) za else ze
nall_P <- function(ro, za, pe) if (grepl("^nall", ro)) 2 * pnorm(-abs(za)) else pe
nall_h <- function(ro, h) if (grepl("^nall", ro)) FALSE else h     # per-removal analytic z not saved at locus level (no mean-shift tile has P < 0.05 there)
for (i in seq_len(nrow(F0))) add(F0$trait[i], "L_all", lmap[[F0$readout[i]]], nall_z(F0$readout[i], F0$z_analytic[i], F0$z_emp[i]),
                                 nall_P(F0$readout[i], F0$z_analytic[i], F0$P_emp[i]), nall_h(F0$readout[i], allh(F0$int_held[i])))
ST <- fread(file.path(V, "joint_bivariate_strata_2026-11-15/joint_bivariate_summary.csv"))[readout %in% names(lmap)]
smap <- c(noUnionDUP = "L_noUD", noUnionAny = "L_noU", noAsym13 = "L_noA")
for (i in seq_len(nrow(ST))) add(ST$trait[i], smap[[ST$stratum[i]]], lmap[[ST$readout[i]]], nall_z(ST$readout[i], ST$z_analytic[i], ST$z_emp[i]),
                                 nall_P(ST$readout[i], ST$z_analytic[i], ST$P_emp[i]), nall_h(ST$readout[i], allh(ST$int_held[i])))
DC <- fread(file.path(V, "joint_bivariate_2026-11-15/mirror_decomposition_union_presence.csv"))[twas == "twU"]
for (i in seq_len(nrow(DC))) {
  add(DC$trait[i], "L_all", "mirror_tot", -DC$z_total[i], 2 * pnorm(-abs(DC$z_total[i])), allh(DC$z_total_int_held[i]))
  add(DC$trait[i], "L_all", "mirror_res", -DC$z_resid[i], 2 * pnorm(-abs(DC$z_resid[i])), allh(DC$z_resid_int_held[i])) }
DS <- fread(file.path(V, "joint_bivariate_strata_2026-11-15/mirror_decomposition_strata_union_presence.csv"))   # growth_mirror_decomp_strata.R
for (i in seq_len(nrow(DS))) { cn <- smap[[DS$stratum[i]]]
  add(DS$trait[i], cn, "mirror_tot", -DS$z_total[i], 2 * pnorm(-abs(DS$z_total[i])), allh(DS$z_total_int_held[i]))
  add(DS$trait[i], cn, "mirror_res", -DS$z_resid[i], 2 * pnorm(-abs(DS$z_resid[i])), allh(DS$z_resid_int_held[i])) }
## UK Biobank
gmap <- function(ro) { m <- sub("^g:", "", ro); tm <- sub("\\|.*", "", m); ar <- sub(".*\\|", "", m)
  if (tm == "n_all") return(paste0("nall_", ar)); if (tm == "none_H") return(if (ar == "mirror") NA else paste0("con_", ar))
  if (tm %in% c("pres:n95", "pres:n95HI")) return(if (ar == "mirror") NA else paste0("hit_", ar)); if (ar == "mirror") "mirror_tw" else paste0("tw_", ar) }
gfiles <- c(G_all = "robust_strata_2026-11-15/S0", G_noUD = "robust_strata_2026-11-15/S1", G_noU = "union_strata_2026-11-15/joint_noUnionAny",
            G_noA = "robust_strata_2026-11-15/S2", G_noUDA = "robust_strata_2026-11-15/S3")
# UK Biobank "All CNVs" = the joint model with the single hit (union presence, rotated with the shift) as a REPORTED term (UKBB session 2026-11-15:
#   union_terms_2026-11-15, height, job 55746822; BMI corrected design union_terms_bmi_fix_2026-11-15, job 55747055, used once present, else S0)
# Jon 2026-11-16: BMI genome-wide TWAS measure = all 45 panels (job 55773555, bmi_allpanels_2026-11-16); height stays on gtex_nonbrain
BMI_ALL <- c(G_all = "G_all", G_noUD = "S1", G_noU = "noUnionAny", G_noA = "S2", G_noUDA = "S3")
# BMI S1-S3 with the union single-hit term (Jon 2026-11-16, job bmi_allpanels_union_2026-11-16) once present, else the top-5%-only runs
gpath <- function(cn, t) { if (nzchar(META)) return(file.path(META, "joint", SUBSET[[cn]], t, "grid_summary.tsv"))
  if (t == "bmi") { fu <- file.path(RA, "bmi_allpanels_union_2026-11-16", BMI_ALL[[cn]], "bmi", "grid_summary.tsv")
    if (cn %in% c("G_noUD", "G_noA", "G_noUDA") && file.exists(fu)) return(fu)
    return(file.path(RA, "bmi_allpanels_2026-11-16", BMI_ALL[[cn]], "bmi", "grid_summary.tsv")) }
  if (cn == "G_all") { f <- file.path(RA, if (t == "height") "union_terms_2026-11-15" else "union_terms_bmi_fix_2026-11-15", t, "grid_summary.tsv")
    if (file.exists(f)) return(f) }
  file.path(RA, gfiles[[cn]], t, paste0(t, "_selected"), "grid_summary.tsv") }
for (cn in names(gfiles)) for (t in c("height", "bmi")) {
  g <- fread(gpath(cn, t)); cat(cn, t, gpath(cn, t), "\n")
  ob <- fread(file.path(dirname(gpath(cn, t)), "per_run_z.tsv"))                # analytic t per run (meta: Stouffer t), for the mean-shift tiles
  for (i in seq_len(nrow(g))) { tm <- gmap(g$readout[i]); if (is.na(tm) || !is.finite(g$z_emp[i])) next
    if (grepl("^nall", tm)) { ta <- g$t_analytic[i]; tl <- ob[readout == g$readout[i] & grepl("^loo_", run), t]
      stopifnot(length(tl) == g$loo_n[i])
      add(if (t == "height") "Height" else "BMI", cn, tm, ta, 2 * pnorm(-abs(ta)), all(abs(tl) > 1.96 & sign(tl) == sign(ta))); next }
    add(if (t == "height") "Height" else "BMI", cn, tm, g$z_emp[i], g$P_emp[i], g$loo_held[i] == g$loo_n[i]) } }
D <- rbindlist(out)
GWL <- if (nzchar(META)) "UK Biobank + All of Us\n(genome-wide)" else "UK Biobank\n(genome-wide)"
D[, q := p.adjust(P, "BH")][, lab := fcase(q < 0.05, "**", P < 0.05, "*", default = "")][, box := P < 0.05 & all_held]
full <- CJ(trait = c("Height", "BMI"), col = CL, term = names(TLAB)); D <- D[full, on = c("trait", "col", "term")]
D[, `:=`(trait = factor(trait, c("Height", "BMI")), col = factor(col, CL), term = factor(term, rev(names(TLAB))),
         level = factor(fifelse(grepl("^L", col), "Recurrent CNVs\n(34 intervals)", GWL), c("Recurrent CNVs\n(34 intervals)", GWL)))]
D[is.na(lab), lab := ""]
# cells empty BY DESIGN get a dash: single-hit carriers removed (locus strata), or single hit fitted as a covariate and not estimated (UK Biobank, all CNVs)
D[!is.finite(z) & ((term == "hit_DUP" & col %in% c("L_noUD", "L_noU")) | (term == "hit_DEL" & col == "L_noU") | (grepl("^hit", term) & col == "G_all")), lab := "–"]
D[, axis := factor(AXIS[as.character(term)], AXL)]
fwrite(D[order(trait, col, -as.integer(term))], file.path(FO, "joint_model_heatmap.csv"))
cat(sprintf("%d tiles with a z; %d P < 0.05; %d q < 0.05; %d boxed\n", D[is.finite(z), .N], D[P < 0.05, .N], D[q < 0.05, .N], D[box == TRUE, .N]))
one <- function(tr, top) { d <- D[trait == tr]
  ggplot(d, aes(col, term, fill = pmax(pmin(z, 4), -4))) + geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = lab), size = 2.6, vjust = 0.75) +
    geom_tile(data = d[box == TRUE], fill = NA, colour = "black", linewidth = 0.55) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-4, 4), na.value = "grey93", name = "z", breaks = c(-4, 0, 4)) +
    facet_grid(axis ~ level, scales = "free", space = "free", switch = "y") +
    scale_x_discrete(labels = CLAB) + scale_y_discrete(labels = TLAB) + labs(x = NULL, y = NULL, title = tr) +
    theme_minimal(base_size = 7) + theme(panel.grid = element_blank(), strip.text.x = if (top) element_text(face = "bold", size = 7.5) else element_blank(),
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 7, hjust = 1), strip.placement = "outside",
      axis.text.x = if (!top) element_text(size = 6, angle = 45, hjust = 1, vjust = 1) else element_blank(), axis.text.y = element_text(size = 6.4),
      plot.title = element_text(face = "bold", size = 9), legend.key.width = unit(6, "pt"), legend.key.height = unit(14, "pt"),
      legend.title = element_text(size = 7), legend.text = element_text(size = 6), panel.spacing = unit(3, "pt")) }
library(patchwork)
p <- (one("Height", TRUE) / one("BMI", FALSE)) + plot_layout(guides = "collect", heights = c(1, 1))
for (ext in c("pdf", "png")) ggsave(file.path(FO, paste0("fig_joint_model_heatmap.", ext)), p, width = 125, height = 135, units = "mm", dpi = 300, bg = "white")
cat("ok\n")
