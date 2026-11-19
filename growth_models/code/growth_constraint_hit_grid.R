# growth_constraint_hit_grid.R -- locus level: which gene criterion (broad s_het > 0.029, s_het tertiles, top-5% s_het, union top-5%/ClinGen HI)
#   and which form (count vs presence) carries the constraint and the single-hit terms (Jon 2026-11-16: "Are counts and presence based on the same
#   criteria? Does the grid look at tertiles of constraint or 5% or both?" ... "And single hit").
#   Joint model and ML fit as growth_hit_spec_compare.R (34 intervals, both arms, arm-specific intercept + total genes + TWAS |z|>3 unsigned count
#   + constraint term(s) + single-hit term(s), correlated locus effects). Every specification is a SEPARATE fit; never a term next to its own subset.
#   A: single hit fixed (union presence), constraint varies: none | broad count | tertile counts (mid, high; low = reference via total genes) | top-5% count
#   B: constraint fixed (broad count), single hit varies: none | union presence | union count | top-5% presence | HI presence | broad presence
#   C: 2x2 criterion x form, one dosage term per arm only (no other constraint/hit term): broad count | broad presence | union count | union presence
#   Readouts: AIC (ML) vs the reference of each block, term z (analytic), family and interval LOO held. Usage: Rscript growth_constraint_hit_grid.R <v2 dir>
suppressMessages(library(data.table))
V <- commandArgs(TRUE)[1]; OUT <- file.path(V, "results/eight_traits/growth_conga/constraint_hit_grid_2026-11-16"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
M_ <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly"; C_ <- file.path(M_, "CNV_ASD_growth_analysis/concordance")
ex <- fread(file.path(V, "inputs/annot/hg38_refGene_20200708.exon.txt"), header = FALSE, col.names = c("chr", "start", "end", "tx", "symbol", "entrez"))
ex <- ex[chr %in% paste0("chr", 1:22)][, entrez := as.character(entrez)]
TXS <- ex[, .(ts = min(start), te = max(end)), by = .(tx, chr, entrez)]
U <- ex[, .(pos = min(start)), by = .(entrez, chr)][, .SD[1], by = entrez][, cn := as.integer(sub("chr", "", chr))][order(cn, pos)]; U[, idx := .I]
gc_ <- fread(file.path(V, "inputs/annot/gene_covariates.tsv"), colClasses = list(character = "entrez")); gc_[shet_imputed == 1, shet := NA]
U <- gc_[, .(entrez, shet)][U, on = "entrez"]; setorder(U, idx); Q95 <- quantile(gc_$shet, .95, na.rm = TRUE)
TT <- quantile(gc_$shet, c(1, 2) / 3, na.rm = TRUE)
cat(sprintf("top-5%% cut %.3f; tertile cuts %.4f %.4f; broad cut 0.029; share of genes with s_het > 0.029: %.2f\n", Q95, TT[1], TT[2], mean(gc_$shet > .029, na.rm = TRUE)))
sym2ez <- rbind(gc_[, .(symbol, entrez)], unique(ex[, .(symbol, entrez)]))[!duplicated(symbol)]
cg <- fread(file.path(V, "inputs/annot/dosage_sensitivity_dl/ClinGen_gene_curation_list_GRCh38.tsv"), skip = 5, sep = "\t", quote = ""); setnames(cg, 1, "gene"); setnames(cg, gsub(" ", "_", names(cg)))
hi_sym <- cg[Haploinsufficiency_Score %in% c("2", "3")]$gene
hi_ez <- unique(c(sym2ez$entrez[match(hi_sym, sym2ez$symbol)], unique(ex[symbol %in% hi_sym]$entrez))); U[, hi := entrez %in% hi_ez]
U[, `:=`(t5 = !is.na(shet) & shet >= Q95, br = !is.na(shet) & shet > .029, tmid = !is.na(shet) & shet >= TT[1] & shet < TT[2], thigh = !is.na(shet) & shet >= TT[2])]
U[, un := t5 | hi]
CO <- unique(fread(file.path(M_, "pathway_analysis/Final/cnv_genes_annotated.csv"))[, .(NAME = CNV_ID, chr = Chromosome, START = CNV_start_GRCh38, END = CNV_end_GRCh38)])
S <- fread(file.path(V, "inputs/growth_input_signed.csv"))
TRAITS <- list(Height = list(lab = "Height", tw = "height"), BMI = list(lab = "BMI (TWAS30)", tw = "bmi"))
famof <- function(nm) { b <- sub("_.*", "", nm); b[grepl("^15q13", b)] <- "15q13"; b }
bifit <- function(y, s, arm, loc, X, ml = TRUE) {
  keep <- sapply(colnames(X), function(cn) { rows <- arm == (if (grepl("_DEL$", cn)) 1 else 2); grepl("^int_", cn) || (sum(rows) > 2 && sd(X[rows, cn]) > 0) })
  X <- X[, keep, drop = FALSE]; if (qr(X)$rank < ncol(X)) return(NULL)
  same <- outer(loc, loc, "=="); A1 <- outer(arm, arm, function(i, j) i == 1 & j == 1) & same; A2 <- outer(arm, arm, function(i, j) i == 2 & j == 2) & same
  A12 <- outer(arm, arm, "!=") & same; S2 <- diag(s^2)
  nll <- function(p) { t1 <- exp(p[1]); t2 <- exp(p[2]); r <- tanh(p[3]); Vm <- S2 + t1^2 * A1 + t2^2 * A2 + r * t1 * t2 * A12
    R <- tryCatch(chol(Vm), error = function(e) NULL); if (is.null(R)) return(1e10); Ri <- backsolve(R, diag(nrow(Vm))); Vi <- Ri %*% t(Ri)
    XtVi <- crossprod(X, Vi); B <- XtVi %*% X; RB <- tryCatch(chol(B), error = function(e) NULL); if (is.null(RB)) return(1e10)
    b <- backsolve(RB, forwardsolve(t(RB), XtVi %*% y)); e <- y - X %*% b
    0.5 * (2 * sum(log(diag(R))) + (if (ml) 0 else 2 * sum(log(diag(RB)))) + drop(crossprod(e, Vi %*% e)) + length(y) * log(2 * pi)) }
  o <- optim(c(log(.1), log(.1), -1), nll, method = "Nelder-Mead", control = list(maxit = 1000, reltol = 1e-10))
  o <- optim(o$par, nll, method = "Nelder-Mead", control = list(maxit = 1000, reltol = 1e-10)); p <- o$par
  Vm <- S2 + exp(p[1])^2 * A1 + exp(p[2])^2 * A2 + tanh(p[3]) * exp(p[1]) * exp(p[2]) * A12; Vi <- solve(Vm)
  Bi <- solve(crossprod(X, Vi %*% X)); b <- drop(Bi %*% crossprod(X, Vi %*% y))
  list(b = setNames(b, colnames(X)), z = setNames(b / sqrt(diag(Bi)), colnames(X)), aic = 2 * o$value + 2 * (ncol(X) + 3), rho = tanh(p[3]), k = ncol(X)) }
# specifications: list(block, name, terms); every term is a per-CNV feature column
SPECS <- list(
  list("A_constraint", "none", "p_un"), list("A_constraint", "broad count (current)", c("c_br", "p_un")),
  list("A_constraint", "tertile counts", c("c_tmid", "c_thigh", "p_un")), list("A_constraint", "top-5% count", c("c_t5", "p_un")),
  list("B_single_hit", "none", "c_br"), list("B_single_hit", "union presence (current)", c("c_br", "p_un")), list("B_single_hit", "union count", c("c_br", "c_un")),
  list("B_single_hit", "top-5% presence", c("c_br", "p_t5")), list("B_single_hit", "HI presence", c("c_br", "p_hi")), list("B_single_hit", "broad presence", c("c_br", "p_br")),
  list("C_2x2", "none", character(0)), list("C_2x2", "broad count", "c_br"), list("C_2x2", "broad presence", "p_br"),
  list("C_2x2", "union count", "c_un"), list("C_2x2", "union presence", "p_un"))
RES <- list()
for (TR in names(TRAITS)) {
  d <- S[disorder == TRAITS[[TR]]$lab, .(NAME, b_DEL = beta_DEL, se_DEL, b_DUP = beta_DUP, se_DUP)][CO, on = "NAME", nomatch = 0]
  d[!(is.finite(se_DEL) & se_DEL > 0), `:=`(b_DEL = NA, se_DEL = NA)]; d[!(is.finite(se_DUP) & se_DUP > 0), `:=`(b_DUP = NA, se_DUP = NA)]
  d <- d[!chr %in% c("chrX", "chrY")]; d[, fam := famof(NAME)]
  GDEL <- lapply(seq_len(nrow(d)), function(i) { x <- d[i]; U[entrez %in% unique(ex[chr == x$chr & end >= x$START & start <= x$END]$entrez)]$idx })
  GDUP <- lapply(seq_len(nrow(d)), function(i) { x <- d[i]; U[entrez %in% unique(TXS[chr == x$chr & ts >= x$START & te <= x$END]$entrez)]$idx })
  TW <- fread(file.path(C_, "twas30", paste0(TRAITS[[TR]]$tw, ".txt")))[is.finite(TWAS.Z), .(z = sum(TWAS.Z) / sqrt(.N)), by = GENE]
  TW[, entrez := sym2ez$entrez[match(GENE, sym2ez$symbol)]]; TW <- TW[!is.na(entrez)][!duplicated(entrez)]; zv <- TW$z[match(U$entrez, TW$entrez)]
  feats <- function(G) { f <- data.table(nall = lengths(G), tw = sapply(G, function(g) sum(abs(zv[g]) > 3, na.rm = TRUE)))
    for (a in c("br", "tmid", "thigh", "t5", "hi", "un")) f[[paste0("c_", a)]] <- sapply(G, function(g) sum(U[[a]][g]))
    for (a in c("br", "t5", "hi", "un")) f[[paste0("p_", a)]] <- as.numeric(f[[paste0("c_", a)]] > 0); f }
  fD <- feats(GDEL); fU <- feats(GDUP)
  iD <- which(is.finite(d$b_DEL)); iU <- which(is.finite(d$b_DUP))
  fit <- function(rows, hv) { rD <- intersect(iD, rows); rU <- intersect(iU, rows)
    y <- c(d$b_DEL[rD], d$b_DUP[rU]); s <- c(d$se_DEL[rD], d$se_DUP[rU]); arm <- c(rep(1, length(rD)), rep(2, length(rU))); loc <- c(rD, rU)
    blk <- function(f, r, k) { m <- cbind(int = 1, nall = f$nall[r], tw = f$tw[r]); if (length(hv)) m <- cbind(m, as.matrix(f[r, hv, with = FALSE]))
      colnames(m) <- paste0(colnames(m), "_", k); m }
    XD <- blk(fD, rD, "DEL"); XU <- blk(fU, rU, "DUP")
    X <- rbind(cbind(XD, matrix(0, nrow(XD), ncol(XU), dimnames = list(NULL, colnames(XU)))), cbind(matrix(0, nrow(XU), ncol(XD), dimnames = list(NULL, colnames(XD))), XU))
    bifit(y, s, arm, loc, X) }
  full <- seq_len(nrow(d)); fams <- unique(d$fam); refaic <- list()
  for (sp in SPECS) { blkn <- sp[[1]]; nm <- sp[[2]]; hv <- sp[[3]]; r <- fit(full, hv)
    if (is.null(r)) { RES[[paste(TR, blkn, nm)]] <- data.table(trait = TR, block = blkn, spec = nm, note = "rank-deficient"); next }
    if (nm == "none") refaic[[blkn]] <- r$aic
    row0 <- data.table(trait = TR, block = blkn, spec = nm, AIC = r$aic, dAIC = r$aic - refaic[[blkn]], rho = r$rho, k = r$k)
    if (!length(hv)) { RES[[paste(TR, blkn, nm)]] <- row0; next }
    for (k in c("DEL", "DUP")) for (h in hv) { key <- paste0(h, "_", k); if (!key %in% names(r$z)) next; z0 <- r$z[[key]]
      zf <- sapply(fams, function(f) { q <- fit(which(d$fam != f), hv); if (is.null(q) || !key %in% names(q$z)) NA else q$z[[key]] })
      zi <- sapply(full, function(i) { q <- fit(setdiff(full, i), hv); if (is.null(q) || !key %in% names(q$z)) NA else q$z[[key]] })
      RES[[paste(TR, blkn, nm, key)]] <- cbind(row0, data.table(term = key, n_with = sum((if (k == "DEL") fD else fU)[[sub("^p_", "c_", h)]] > 0, na.rm = TRUE),
        b = r$b[[key]], z = z0, fam_held = sum(abs(zf) > 1.96 & sign(zf) == sign(z0), na.rm = TRUE), fam_n = sum(is.finite(zf)),
        fam_weakest = names(zf)[which.min(abs(zf))], int_held = sum(abs(zi) > 1.96 & sign(zi) == sign(z0), na.rm = TRUE), int_n = sum(is.finite(zi)))) } }
  cat(TR, "done\n"); flush.console()
}
R <- rbindlist(RES, fill = TRUE); fwrite(R, file.path(OUT, "constraint_hit_grid.csv"))
options(width = 250); print(R[, .(trait, block, spec, term, dAIC = round(dAIC, 1), z = round(z, 2), fam = paste0(fam_held, "/", fam_n), int = paste0(int_held, "/", int_n))])
