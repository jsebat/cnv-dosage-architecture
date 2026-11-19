# growth_mirror_decomp_strata.R -- growth_mirror_decomp.R (unsigned TWAS count, union-presence hit) refitted in the three locus-level STRATA of
#   growth_joint_bivariate_strata.R (Jon 2026-11-15: confirm the mirror is robust to the relevant drops): noUnionDUP (DUP arm dropped for intervals whose
#   DUP carries a union gene), noUnionAny, noAsym13; each with every family and single-interval removal inside the stratum. Original header:
# growth_mirror_decomp.R -- the true mirror = TWAS mirror + rho (Jon 2026-11-15: "the Del-Dup correlation explained by TWAS is the convergent evidence for the
#   mirror from common variants. The mirror itself is something bigger and TWAS can't explain all of it. So the true mirror is TWAS mirror + rho").
#   Joint model as growth_joint_bivariate.R (both arms stacked; arm-specific intercept + total genes + constrained-gene count + single hit [+ TWAS count];
#   correlated locus random effects). Two fits per trait, ML:
#     without the TWAS term -> rho0, cov0 = rho0 tau_DEL tau_DUP: the total mirror left after constraint and the single hit;
#     with the TWAS term    -> rho1, cov1: the part TWAS does not explain.  TWAS share = 1 - cov1 / cov0.
#   Tests for rho: LRT against rho = 0 (1 df) in each fit; every single-interval removal and every family removal (rho, LRT z, share).
#   rho is not tested by the genome-shift null: rotating gene annotations does not move the pairing of deletion and duplication effects.
#   Usage: Rscript growth_mirror_decomp.R <v2 dir> [hit spec: union_presence (default) | union_count | top5_presence | top5_count]
suppressMessages(library(data.table))
a <- commandArgs(TRUE); V <- a[1]; HS <- if (length(a) >= 2) a[2] else "union_presence"
OUT <- file.path(V, "results/eight_traits/growth_conga/joint_bivariate_strata_2026-11-15")
M_ <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly"; C_ <- file.path(M_, "CNV_ASD_growth_analysis/concordance")
ex <- fread(file.path(V, "inputs/annot/hg38_refGene_20200708.exon.txt"), header = FALSE, col.names = c("chr", "start", "end", "tx", "symbol", "entrez"))
ex <- ex[chr %in% paste0("chr", 1:22)][, entrez := as.character(entrez)]
TXS <- ex[, .(ts = min(start), te = max(end)), by = .(tx, chr, entrez)]
U <- ex[, .(pos = min(start)), by = .(entrez, chr)][, .SD[1], by = entrez][, cn := as.integer(sub("chr", "", chr))][order(cn, pos)]; U[, idx := .I]
gc_ <- fread(file.path(V, "inputs/annot/gene_covariates.tsv"), colClasses = list(character = "entrez")); gc_[shet_imputed == 1, shet := NA]
U <- gc_[, .(entrez, shet)][U, on = "entrez"]; setorder(U, idx); Q95 <- quantile(gc_$shet, .95, na.rm = TRUE)
sym2ez <- rbind(gc_[, .(symbol, entrez)], unique(ex[, .(symbol, entrez)]))[!duplicated(symbol)]
cg <- fread(file.path(V, "inputs/annot/dosage_sensitivity_dl/ClinGen_gene_curation_list_GRCh38.tsv"), skip = 5, sep = "\t", quote = ""); setnames(cg, 1, "gene"); setnames(cg, gsub(" ", "_", names(cg)))
hi_sym <- cg[Haploinsufficiency_Score %in% c("2", "3")]$gene
hi_ez <- unique(c(sym2ez$entrez[match(hi_sym, sym2ez$symbol)], unique(ex[symbol %in% hi_sym]$entrez))); U[, hi := entrez %in% hi_ez]
U[, t5 := !is.na(shet) & shet >= Q95]; U[, un := t5 | hi]
CO <- unique(fread(file.path(M_, "pathway_analysis/Final/cnv_genes_annotated.csv"))[, .(NAME = CNV_ID, chr = Chromosome, START = CNV_start_GRCh38, END = CNV_end_GRCh38)])
S <- fread(file.path(V, "inputs/growth_input_signed.csv"))
TRAITS <- list(Height = list(lab = "Height", tw = "height"), BMI = list(lab = "BMI (TWAS30)", tw = "bmi"))
famof <- function(nm) { b <- sub("_.*", "", nm); b[grepl("^15q13", b)] <- "15q13"; b }
bifit <- function(y, s, arm, loc, X, rho0 = FALSE) {
  keep <- sapply(colnames(X), function(cn) { rows <- arm == (if (grepl("_DEL$", cn)) 1 else 2); grepl("^int_", cn) || (sum(rows) > 2 && sd(X[rows, cn]) > 0) })
  X <- X[, keep, drop = FALSE]; if (qr(X)$rank < ncol(X)) return(NULL)
  same <- outer(loc, loc, "=="); A1 <- outer(arm, arm, function(i, j) i == 1 & j == 1) & same; A2 <- outer(arm, arm, function(i, j) i == 2 & j == 2) & same
  A12 <- outer(arm, arm, "!=") & same; S2 <- diag(s^2)
  nll <- function(p) { if (rho0) p <- c(p, 0); t1 <- exp(p[1]); t2 <- exp(p[2]); r <- tanh(p[3]); Vm <- S2 + t1^2 * A1 + t2^2 * A2 + r * t1 * t2 * A12
    R <- tryCatch(chol(Vm), error = function(e) NULL); if (is.null(R)) return(1e10); Ri <- backsolve(R, diag(nrow(Vm))); Vi <- Ri %*% t(Ri)
    XtVi <- crossprod(X, Vi); B <- XtVi %*% X; RB <- tryCatch(chol(B), error = function(e) NULL); if (is.null(RB)) return(1e10)
    b <- backsolve(RB, forwardsolve(t(RB), XtVi %*% y)); e <- y - X %*% b
    0.5 * (2 * sum(log(diag(R))) + drop(crossprod(e, Vi %*% e)) + length(y) * log(2 * pi)) }
  st <- if (rho0) c(log(.1), log(.1)) else c(log(.1), log(.1), -1)
  best <- NULL; for (s0 in list(st, st + .7, st - .7)) { o <- optim(s0, nll, method = "Nelder-Mead", control = list(maxit = 2000, reltol = 1e-11))
    o <- optim(o$par, nll, method = "Nelder-Mead", control = list(maxit = 2000, reltol = 1e-11)); if (is.null(best) || o$value < best$value) best <- o }
  p <- if (rho0) c(best$par, 0) else best$par
  list(nll = best$value, tau = exp(p[1:2]), rho = tanh(p[3]), cov = tanh(p[3]) * exp(p[1]) * exp(p[2])) }
RES <- list(); LO <- list()
for (TR in names(TRAITS)) {
  d <- S[disorder == TRAITS[[TR]]$lab, .(NAME, b_DEL = beta_DEL, se_DEL, b_DUP = beta_DUP, se_DUP)][CO, on = "NAME", nomatch = 0]
  d[!(is.finite(se_DEL) & se_DEL > 0), `:=`(b_DEL = NA, se_DEL = NA)]; d[!(is.finite(se_DUP) & se_DUP > 0), `:=`(b_DUP = NA, se_DUP = NA)]
  d <- d[!chr %in% c("chrX", "chrY")]; d[, fam := famof(NAME)]
  GDEL <- lapply(seq_len(nrow(d)), function(i) { x <- d[i]; U[entrez %in% unique(ex[chr == x$chr & end >= x$START & start <= x$END]$entrez)]$idx })
  GDUP <- lapply(seq_len(nrow(d)), function(i) { x <- d[i]; U[entrez %in% unique(TXS[chr == x$chr & ts >= x$START & te <= x$END]$entrez)]$idx })
  TW <- fread(file.path(C_, "twas30", paste0(TRAITS[[TR]]$tw, ".txt")))[is.finite(TWAS.Z), .(z = sum(TWAS.Z) / sqrt(.N)), by = GENE]
  TW[, entrez := sym2ez$entrez[match(GENE, sym2ez$symbol)]]; TW <- TW[!is.na(entrez)][!duplicated(entrez)]; zv <- TW$z[match(U$entrez, TW$entrez)]
  hitf <- function(g) { set <- if (grepl("^union", HS)) U$un[g] else U$t5[g]; if (grepl("presence$", HS)) as.numeric(any(set)) else sum(set) }
  feats <- function(G) data.table(nall = lengths(G), twU = sapply(G, function(g) sum(abs(zv[g]) > 3, na.rm = TRUE)),
    twA = sapply(G, function(g) sum(zv[g] > 3, na.rm = TRUE) - sum(zv[g] < -3, na.rm = TRUE)), con = sapply(G, function(g) sum(U$shet[g] > .029, na.rm = TRUE)), hit = sapply(G, hitf))
  fD <- feats(GDEL); fU <- feats(GDUP); iD <- which(is.finite(d$b_DEL)); iU <- which(is.finite(d$b_DUP))
  fit <- function(rows, tw, rho0 = FALSE) { rD <- intersect(iD, rows$D); rU <- intersect(iU, rows$U)   # rows = list(D, U)
    y <- c(d$b_DEL[rD], d$b_DUP[rU]); s <- c(d$se_DEL[rD], d$se_DUP[rU]); arm <- c(rep(1, length(rD)), rep(2, length(rU))); loc <- c(rD, rU)
    cols <- c("nall", "con", "hit", if (!is.null(tw)) tw)
    blk <- function(f, r, k) { m <- cbind(int = 1, as.matrix(f[r, cols, with = FALSE])); colnames(m) <- paste0(colnames(m), "_", k); m }
    XD <- blk(fD, rD, "DEL"); XU <- blk(fU, rU, "DUP")
    X <- rbind(cbind(XD, matrix(0, nrow(XD), ncol(XU), dimnames = list(NULL, colnames(XU)))), cbind(matrix(0, nrow(XU), ncol(XD), dimnames = list(NULL, colnames(XD))), XU))
    bifit(y, s, arm, loc, X, rho0) }
  dec <- function(rows, tw) { f0 <- fit(rows, NULL); f0z <- fit(rows, NULL, TRUE); f1 <- fit(rows, tw); f1z <- fit(rows, tw, TRUE)
    lrt <- function(a, b) { x <- max(0, 2 * (b$nll - a$nll)); sign(a$rho) * sqrt(x) }
    c(rho_total = f0$rho, cov_total = f0$cov, z_total = lrt(f0, f0z), rho_resid = f1$rho, cov_resid = f1$cov, z_resid = lrt(f1, f1z), twas_share = 1 - f1$cov / f0$cov,
      tau_DEL0 = f0$tau[1], tau_DUP0 = f0$tau[2], tau_DEL1 = f1$tau[1], tau_DUP1 = f1$tau[2]) }
  full <- seq_len(nrow(d)); tw <- "twU"
  ASYM <- fread("C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/UKBB_anthro_CNV/code/resources_local/growth_loci/asymmetric_molly13.tsv")$NAME; stopifnot(all(ASYM %in% d$NAME))
  STR <- list(noUnionDUP = list(D = full, U = which(fU$hit == 0)), noUnionAny = list(D = which(fD$hit == 0), U = which(fU$hit == 0)),
              noAsym13 = list(D = which(!d$NAME %in% ASYM), U = which(!d$NAME %in% ASYM)))
  for (k in names(STR)) { st0 <- STR[[k]]; inn <- union(st0$D, st0$U); fams <- unique(d$fam[inn])
    F <- dec(st0, tw); FA <- t(sapply(fams, function(f) dec(lapply(st0, function(r) r[d$fam[r] != f]), tw)))
    IA <- t(sapply(inn, function(i) dec(lapply(st0, function(r) setdiff(r, i)), tw))); rownames(IA) <- d$NAME[inn]
    held <- function(M, kk) sum(abs(M[, kk]) > 1.96 & sign(M[, kk]) == sign(F[[kk]]), na.rm = TRUE)
    RES[[paste(TR, k)]] <- data.table(trait = TR, stratum = k, twas = tw, hit = HS, t(F), z_total_fam_held = paste0(held(FA, "z_total"), "/", nrow(FA)), z_total_int_held = paste0(held(IA, "z_total"), "/", nrow(IA)),
      z_resid_fam_held = paste0(held(FA, "z_resid"), "/", nrow(FA)), z_resid_int_held = paste0(held(IA, "z_resid"), "/", nrow(IA)),
      z_total_min_int = min(IA[, "z_total"] * sign(F[["z_total"]])), weakest_int_total = rownames(IA)[which.min(IA[, "z_total"] * sign(F[["z_total"]]))],
      z_resid_min_int = min(IA[, "z_resid"] * sign(F[["z_resid"]])), weakest_int_resid = rownames(IA)[which.min(IA[, "z_resid"] * sign(F[["z_resid"]]))],
      share_range_int = sprintf("%.2f to %.2f", min(IA[, "twas_share"]), max(IA[, "twas_share"])))
    LO[[paste(TR, k)]] <- rbind(data.table(trait = TR, stratum = k, removal = "family", dropped = rownames(FA), FA), data.table(trait = TR, stratum = k, removal = "interval", dropped = rownames(IA), IA))
    cat(TR, k, "done\n"); flush.console() }
}
R <- rbindlist(RES); fwrite(R, file.path(OUT, sprintf("mirror_decomposition_strata_%s.csv", HS))); fwrite(rbindlist(LO), file.path(OUT, sprintf("mirror_decomposition_strata_loo_%s.csv", HS)))
options(width = 250); print(R[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)])
