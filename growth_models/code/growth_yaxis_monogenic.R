# growth_yaxis_monogenic.R -- height duplication-axis (y) shift: monogenic (DDG2P) genes vs constraint (top-5% s_het) (Jon 2026-11-14: "make sure you have the ddg2p test.
#   There is a convergence of evidence for a monogenic effect along the dup axis that's only partly related to constraint").
#   Offset of each interval's duplication effect from the mirror line, off = beta_DUP + 0.351 beta_DEL (fixed line through the origin, Figure 3).
#   Gene sets in the duplicated interval (refGene transcript wholly inside; symbol match incl. HGNC previous symbols):
#     top-5% s_het (GeneBayes >= genome 95th percentile); DDG2P monoallelic autosomal definitive/strong: any mechanism; loss of function (absent / decreased
#     product); non-LoF (gain of function, dominant negative, increased / altered product); skeletal panel (any mechanism); DDG2P genes NOT in the top-5% s_het.
#   Tests: REML mean offset in intervals with vs without each set (difference z); joint fit off ~ presence(top-5%) + presence(DDG2P set); every single-interval
#   removal (held = |z| > 1.96, same sign). Height and BMI. Run from the v2 dir. Output yaxis_monogenic*.csv.
suppressMessages(library(data.table)); options(width = 220)
M_ <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly"
ex <- fread("inputs/annot/hg38_refGene_20200708.exon.txt", header = FALSE, col.names = c("chr", "start", "end", "tx", "symbol", "entrez"))[chr %in% paste0("chr", 1:22)][, entrez := as.character(entrez)]
TXS <- ex[, .(ts = min(start), te = max(end)), by = .(tx, chr, entrez, symbol)]
gc_ <- fread("inputs/annot/gene_covariates.tsv", colClasses = list(character = "entrez")); gc_[shet_imputed == 1, shet := NA]; Q95 <- quantile(gc_$shet, .95, na.rm = TRUE)
dd <- fread("inputs/annot/DDG2P.csv"); setnames(dd, gsub(" ", "_", names(dd)))
dd <- dd[allelic_requirement == "monoallelic_autosomal" & confidence %in% c("definitive", "strong")]
dd[, syms := paste(gene_symbol, previous_gene_symbols, sep = ";")]
symset <- function(x) unique(trimws(unlist(strsplit(x$syms, ";"))))
lofm <- function(x) x[grepl("loss of function", molecular_mechanism) | grepl("absent gene product|decreased gene product", variant_consequence)]
SETS <- list(dd_any = symset(dd), dd_lof = symset(lofm(dd)), dd_nonlof = symset(dd[!gene_symbol %in% lofm(dd)$gene_symbol]),
             dd_skeletal = symset(dd[grepl("Skeletal", panel)]))
## dosage-sensitivity scores (downloaded 2026-11-14, Jon OK; "a reviewer is probably going to ask"): Collins et al. 2022 pHaplo / pTriplo (Zenodo 6347673; cut-offs from the paper,
##   pHaplo >= 0.86, pTriplo >= 0.94) and ClinGen dosage curation (score 3 sufficient or 2 emerging evidence)
DS <- "inputs/annot/dosage_sensitivity_dl"
co <- fread(file.path(DS, "Collins_rCNV_2022.dosage_sensitivity_scores.tsv")); setnames(co, 1, "gene")
cg <- fread(file.path(DS, "ClinGen_gene_curation_list_GRCh38.tsv"), skip = 5, sep = "	", quote = ""); setnames(cg, 1, "gene"); setnames(cg, gsub(" ", "_", names(cg)))
SETS$triplo_collins <- co[pTriplo >= .94]$gene; SETS$haplo_collins <- co[pHaplo >= .86]$gene
SETS$ts_clingen <- cg[Triplosensitivity_Score %in% c("2", "3")]$gene; SETS$hi_clingen <- cg[Haploinsufficiency_Score %in% c("2", "3")]$gene
cat(sprintf("Collins pTriplo >= 0.94: %d genes; pHaplo >= 0.86: %d; ClinGen TS 2-3: %d; HI 2-3: %d
", length(SETS$triplo_collins), length(SETS$haplo_collins), length(SETS$ts_clingen), length(SETS$hi_clingen)))
cat(sprintf("DDG2P monoallelic definitive/strong genes: %d (LoF %d, non-LoF only %d, skeletal panel %d); top-5%% s_het cut %.3f\n",
            uniqueN(dd$gene_symbol), uniqueN(lofm(dd)$gene_symbol), uniqueN(dd[!gene_symbol %in% lofm(dd)$gene_symbol]$gene_symbol), uniqueN(dd[grepl("Skeletal", panel)]$gene_symbol), Q95))
CO <- unique(fread(file.path(M_, "pathway_analysis/Final/cnv_genes_annotated.csv"))[, .(NAME = CNV_ID, chr = Chromosome, START = CNV_start_GRCh38, END = CNV_end_GRCh38)])
L0 <- fread("results/eight_traits/growth_conga/growth_conga_loci.csv")
rema <- function(y, s, X) { rl <- function(t2) { w <- 1/(s^2+t2); A <- t(X*w)%*%X; b <- solve(A, t(X*w)%*%y); r <- y - X%*%b; -0.5*(sum(log(s^2+t2)) + determinant(A)$modulus + sum(w*r^2)) }
  t2 <- optimize(function(t) -rl(t), c(0, 4*var(y)))$minimum; w <- 1/(s^2+t2); A <- t(X*w)%*%X; b <- solve(A, t(X*w)%*%y); se <- sqrt(diag(solve(A)))
  data.table(term = colnames(X), b = drop(b), se = se, z = drop(b)/se) }
RES <- list(); GEN <- list(); PERL <- list()
for (TR in c("Height", "BMI")) {
  L <- CO[L0[trait == TR], on = "NAME"][is.finite(b_DEL) & is.finite(b_DUP)]
  G <- lapply(seq_len(nrow(L)), function(i) { x <- L[i]; unique(TXS[chr == x$chr & ts >= x$START & te <= x$END, .(entrez, symbol)]) })
  for (i in seq_len(nrow(L))) { g <- G[[i]]; g[, shet := gc_$shet[match(entrez, gc_$entrez)]]; G[[i]] <- g }
  ind <- function(fun) sapply(G, function(g) as.numeric(any(fun(g))))
  L[, top5 := ind(function(g) !is.na(g$shet) & g$shet >= Q95)]
  for (s in names(SETS)) set(L, j = s, value = ind(function(g) g$symbol %in% SETS[[s]]))
  L[, dd_not_top5 := ind(function(g) g$symbol %in% SETS$dd_any & !(!is.na(g$shet) & g$shet >= Q95))]
  L[, triplo_not_top5 := ind(function(g) g$symbol %in% SETS$triplo_collins & !(!is.na(g$shet) & g$shet >= Q95))]
  L[, top5_or_hi := ind(function(g) (!is.na(g$shet) & g$shet >= Q95) | g$symbol %in% SETS$hi_clingen)]   # union, Jon 2026-11-14
  L[, top5_or_hi_or_ddlof := ind(function(g) (!is.na(g$shet) & g$shet >= Q95) | g$symbol %in% SETS$hi_clingen | g$symbol %in% SETS$dd_lof)]
  L[, off := b_DUP + .351 * b_DEL]; L[, se_off := sqrt(se_DUP^2 + (.351 * se_DEL)^2)]
  VARS <- c("top5", names(SETS), "dd_not_top5", "triplo_not_top5", "top5_or_hi", "top5_or_hi_or_ddlof")
  for (v in VARS) { x <- L[[v]]; if (sum(x) < 2 || sum(1 - x) < 2) next
    f <- rema(L$off, L$se_off, cbind(without = 1 - x, with = x)); d <- rema(L$off, L$se_off, cbind(int = 1, with = x))
    lo <- sapply(seq_len(nrow(L)), function(i) { xi <- x[-i]; if (sum(xi) < 1 || sum(1 - xi) < 1) NA else rema(L$off[-i], L$se_off[-i], cbind(int = 1, with = xi))$z[2] })
    RES[[paste(TR, v)]] <- data.table(trait = TR, set = v, model = "alone", intervals_with = sum(x), intervals = nrow(L), offset_with = f$b[2], offset_without = f$b[1],
      z_without = f$z[1], diff = d$b[2], z_diff = d$z[2], loo_held = sum(abs(lo) > 1.96 & sign(lo) == sign(d$z[2]), na.rm = TRUE)) }
  for (v in c(names(SETS), "dd_not_top5", "triplo_not_top5")) { x <- L[[v]]; if (sum(x) < 2) next
    j <- rema(L$off, L$se_off, cbind(int = 1, top5 = L$top5, dd = x))
    RES[[paste(TR, v, "joint")]] <- data.table(trait = TR, set = v, model = "jointly with top-5% presence", intervals_with = sum(x), intervals = nrow(L),
      z_diff = j$z[3], diff = j$b[3], z_top5_alongside = j$z[2], top5_alongside = j$b[2]) }
  # residual: among intervals without a top-5% gene, does DDG2P presence explain the offset?
  y <- L[top5 == 0]; for (v in c("dd_any", "dd_lof", "dd_nonlof", "dd_skeletal", "triplo_collins", "haplo_collins", "ts_clingen", "hi_clingen")) { x <- y[[v]]; if (sum(x) < 2 || sum(1 - x) < 2) next
    d <- rema(y$off, y$se_off, cbind(int = 1, with = x)); f <- rema(y$off, y$se_off, cbind(without = 1 - x, with = x))
    RES[[paste(TR, v, "nontop5")]] <- data.table(trait = TR, set = v, model = "within intervals without a top-5% gene", intervals_with = sum(x), intervals = nrow(y),
      offset_with = f$b[2], offset_without = f$b[1], z_without = f$z[1], diff = d$b[2], z_diff = d$z[2]) }
  GEN[[TR]] <- rbindlist(lapply(seq_len(nrow(L)), function(i) { g <- G[[i]]
    data.table(trait = TR, NAME = L$NAME[i], off = L$off[i], se_off = L$se_off[i], top5_genes = paste(g[!is.na(shet) & shet >= Q95]$symbol, collapse = " "),
               ddg2p_genes = paste(unique(g[symbol %in% SETS$dd_any]$symbol), collapse = " "), ddg2p_lof = paste(unique(g[symbol %in% SETS$dd_lof]$symbol), collapse = " "),
               ddg2p_nonlof = paste(unique(g[symbol %in% SETS$dd_nonlof]$symbol), collapse = " "), ddg2p_skeletal = paste(unique(g[symbol %in% SETS$dd_skeletal]$symbol), collapse = " "),
               triplo_collins = paste(unique(g[symbol %in% SETS$triplo_collins]$symbol), collapse = " "), ts_clingen = paste(unique(g[symbol %in% SETS$ts_clingen]$symbol), collapse = " "),
               hi_clingen = paste(unique(g[symbol %in% SETS$hi_clingen]$symbol), collapse = " ")) }))
  PERL[[TR]] <- L[, c("trait", "NAME", "off", "se_off", VARS), with = FALSE]
}
R <- rbindlist(RES, fill = TRUE); fwrite(R, "results/eight_traits/growth_conga/yaxis_monogenic.csv"); fwrite(rbindlist(GEN), "results/eight_traits/growth_conga/yaxis_monogenic_genes.csv")
fwrite(rbindlist(PERL), "results/eight_traits/growth_conga/yaxis_monogenic_intervals.csv")
print(R[, .(trait, set, model = substr(model, 1, 24), k = paste0(intervals_with, "/", intervals), with = round(offset_with, 3), without = round(offset_without, 3), z_wo = round(z_without, 2),
            diff = round(diff, 3), z = round(z_diff, 2), loo = loo_held, z_top5 = round(z_top5_alongside, 2))])
print(rbindlist(GEN)[trait == "Height"][order(off), .(NAME, off = round(off, 2), top5_genes, ddg2p_genes, triplo_collins, ts_clingen, hi_clingen)])
