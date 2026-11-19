# growth_conga_axes.R -- height and BMI on the ConGA axis definitions (Jon 2026-11-14: "create a totally new version of the report using the conga defined
#   definitions and axes ... add mean shift as a 3rd axis ... include the presence/absence vs count ad hoc analysis that distinguishes the asymmetric height CNVs").
#   Model = locus_count_mirror.R (the disorder fits) unchanged apart from the inputs:
#     per trait x dosage, REML random-effects meta-regression of the locus effect (SD units) on intercept + total gene count + count_A + none_H (uncentred).
#     count_A = n(z > 3) - n(z < -3) over all genes hit; gene z = GTEx TWAS30 z combined across tissues as sum(z)/sqrt(n tissues) ("scaled", Jon 2026-11-12;
#     the mean-tissue z leaves the |z| > 3 count almost empty). none_H = n genes with GeneBayes s_het > 0.029 (shet_imputed -> NA).
#     CNV-to-gene: DEL = any refGene exon >= 1 bp; DUP = every exon of at least one refGene transcript inside the interval. Intervals = Molly's Sacks intervals (GRCh38).
#     Readouts: mirror (z_DUP - z_DEL)/sqrt2 and rotated mirror (z_DUP + z_DEL)/sqrt2 for count_A; constraint z_DEL, z_DUP for none_H;
#     mean shift: intercept (count-independent) and total-count slope (count-dependent) per dosage; both-dosage shift = smaller |z| if signs agree.
#   Null for count_A / none_H: genomic shift of gene annotations (TWAS z and s_het rotated together), 2,000 replicates, empirical z; every single-locus LOO matched.
#   Mean shift: analytic z (the shift does not move CNV size), every single-locus LOO.
#   Ad hoc presence/absence vs count (height's asymmetric duplication loci): on each dosage, the constraint term none_H replaced by (a) presence of >= 1 top-5% s_het
#     gene, (b) the top-5% gene count, (c) none_H + presence; ML log-likelihoods (fixed effects differ) -> AIC; LOO for the presence term.
#   Public data; runs locally. Usage: Rscript growth_conga_axes.R <v2 dir> [nrep]
suppressMessages(library(data.table))
a <- commandArgs(TRUE); V <- a[1]; NREP <- if (length(a) >= 2) as.integer(a[2]) else 2000L; set.seed(20261114)
OUT <- file.path(V, "results/eight_traits/growth_conga"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
M_ <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly"; C_ <- file.path(M_, "CNV_ASD_growth_analysis/concordance")
## genes (as locus_count_mirror.R)
ex <- fread(file.path(V, "inputs/annot/hg38_refGene_20200708.exon.txt"), header = FALSE, col.names = c("chr", "start", "end", "tx", "symbol", "entrez"))
ex <- ex[chr %in% paste0("chr", 1:22)][, entrez := as.character(entrez)]
TXS <- ex[, .(ts = min(start), te = max(end)), by = .(tx, chr, entrez)]
U <- ex[, .(pos = min(start)), by = .(entrez, chr)][, .SD[1], by = entrez][, cn := as.integer(sub("chr", "", chr))][order(cn, pos)]; U[, idx := .I]
gc_ <- fread(file.path(V, "inputs/annot/gene_covariates.tsv"), colClasses = list(character = "entrez")); gc_[shet_imputed == 1, shet := NA]
U <- gc_[, .(entrez, shet)][U, on = "entrez"]; setorder(U, idx)
Q95 <- quantile(gc_$shet, .95, na.rm = TRUE); cat(sprintf("top-5%% s_het cut %.3f\n", Q95))
sym2ez <- rbind(gc_[, .(symbol, entrez)], unique(ex[, .(symbol, entrez)]))[!duplicated(symbol)]
## loci
CO <- unique(fread(file.path(M_, "pathway_analysis/Final/cnv_genes_annotated.csv"))[, .(NAME = CNV_ID, chr = Chromosome, START = CNV_start_GRCh38, END = CNV_end_GRCh38)])
S <- fread(file.path(V, "inputs/growth_input_signed.csv"))
TRAITS <- list(Height = list(lab = "Height", tw = "height"), BMI = list(lab = "BMI (TWAS30)", tw = "bmi"))
rema <- function(y, s, X, ml = FALSE) { ok <- is.finite(y) & is.finite(s); y <- y[ok]; s <- s[ok]; X <- X[ok, , drop = FALSE]
  keep <- apply(X, 2, function(v) sd(v) > 0); keep[1] <- TRUE; X <- X[, keep, drop = FALSE]
  if (qr(X)$rank < ncol(X) || nrow(X) <= ncol(X) + 1) return(NULL)
  rl <- function(t2) { w <- 1 / (s^2 + t2); XtW <- t(X * w); A <- XtW %*% X; b <- solve(A, XtW %*% y); r <- y - X %*% b
    -0.5 * (sum(log(s^2 + t2)) + (if (ml) 0 else as.numeric(determinant(A)$modulus)) + sum(w * r^2)) }
  t2 <- optimize(function(t) -rl(t), c(0, max(4 * var(y), 1e-4)))$minimum
  w <- 1 / (s^2 + t2); A <- t(X * w) %*% X; b <- solve(A, t(X * w) %*% y); Vb <- solve(A)
  z <- setNames(rep(NA_real_, length(keep)), names(keep)); z[keep] <- b / sqrt(diag(Vb)); bb <- z; bb[keep] <- b
  list(z = z, b = bb, ll = rl(t2) - if (ml) .5 * sum(ok) * log(2 * pi) else 0, k = ncol(X) + 1, tau2 = t2) }
RD <- c("mirror", "rotated", "countA_DEL", "countA_DUP", "con_DEL", "con_DUP", "mirrorU", "countU_DEL", "countU_DUP", "conU_DEL", "conU_DUP")   # *U = sensitivity, unsigned |z| > 3 count in place of count_A
SUMM <- list(); MS <- list(); HIT <- list(); LOCI <- list(); CONT <- list()
for (TR in names(TRAITS)) {
  d <- S[disorder == TRAITS[[TR]]$lab, .(NAME, b_DEL = beta_DEL, se_DEL, b_DUP = beta_DUP, se_DUP)][CO, on = "NAME", nomatch = 0]
  d[!(is.finite(se_DEL) & se_DEL > 0), `:=`(b_DEL = NA, se_DEL = NA)]; d[!(is.finite(se_DUP) & se_DUP > 0), `:=`(b_DUP = NA, se_DUP = NA)]
  d <- d[!chr %in% c("chrX", "chrY")]
  GDEL <- lapply(seq_len(nrow(d)), function(i) { x <- d[i]; U[entrez %in% unique(ex[chr == x$chr & end >= x$START & start <= x$END]$entrez)]$idx })
  GDUP <- lapply(seq_len(nrow(d)), function(i) { x <- d[i]; U[entrez %in% unique(TXS[chr == x$chr & ts >= x$START & te <= x$END]$entrez)]$idx })
  TW <- fread(file.path(C_, "twas30", paste0(TRAITS[[TR]]$tw, ".txt")))[is.finite(TWAS.Z), .(z = sum(TWAS.Z) / sqrt(.N)), by = GENE]
  TW[, entrez := sym2ez$entrez[match(GENE, sym2ez$symbol)]]; TW <- TW[!is.na(entrez)][!duplicated(entrez)]
  zv <- TW$z[match(U$entrez, TW$entrez)]; sv <- U$shet; n <- length(zv)
  i16 <- which(d$NAME %in% c("16p11.2_BP4_BP5", "16p11.2_BP4-BP5"))
  cat(sprintf("%s: %d loci (DEL %d, DUP %d); genes with TWAS z %d\n", TR, nrow(d), sum(is.finite(d$b_DEL)), sum(is.finite(d$b_DUP)), sum(is.finite(zv))))
  feats <- function(G, zz, ss) data.table(nall = lengths(G), countA = sapply(G, function(g) sum(zz[g] > 3, na.rm = TRUE) - sum(zz[g] < -3, na.rm = TRUE)),
                                          noneH = sapply(G, function(g) sum(ss[g] > .029, na.rm = TRUE)), n95 = sapply(G, function(g) sum(ss[g] >= Q95, na.rm = TRUE)),
                                          countU = sapply(G, function(g) sum(abs(zz[g]) > 3, na.rm = TRUE)))
  fitset <- function(zz, ss, rows) { out <- list()
    for (al in c("DEL", "DUP")) { f <- feats((if (al == "DEL") GDEL else GDUP)[rows], zz, ss)
      r <- rema(d[[paste0("b_", al)]][rows], d[[paste0("se_", al)]][rows], cbind(int = 1, nall = f$nall, countA = f$countA, noneH = f$noneH))
      out[[al]] <- if (is.null(r)) c(int = NA, nall = NA, countA = NA, noneH = NA) else r$z
      q <- rema(d[[paste0("b_", al)]][rows], d[[paste0("se_", al)]][rows], cbind(int = 1, nall = f$nall, countU = f$countU, noneH = f$noneH))   # sensitivity: unsigned count
      out[[paste0(al, "U")]] <- if (is.null(q)) c(int = NA, nall = NA, countU = NA, noneH = NA) else q$z }
    c(mirror = (out$DUP[["countA"]] - out$DEL[["countA"]]) / sqrt(2), rotated = (out$DUP[["countA"]] + out$DEL[["countA"]]) / sqrt(2),
      countA_DEL = out$DEL[["countA"]], countA_DUP = out$DUP[["countA"]], con_DEL = out$DEL[["noneH"]], con_DUP = out$DUP[["noneH"]],
      mirrorU = (out$DUPU[["countU"]] - out$DELU[["countU"]]) / sqrt(2), countU_DEL = out$DELU[["countU"]], countU_DUP = out$DUPU[["countU"]],
      conU_DEL = out$DELU[["noneH"]], conU_DUP = out$DUPU[["noneH"]]) }
  full <- seq_len(nrow(d)); sets <- c(list(full = full, no16p = setdiff(full, i16)), setNames(lapply(full, function(i) setdiff(full, i)), paste0("loo_", d$NAME)))
  OBS <- sapply(sets, function(r) fitset(zv, sv, r))
  NUL <- array(NA_real_, c(length(RD), length(sets), NREP), dimnames = list(RD, names(sets), NULL))
  for (k in seq_len(NREP)) { sh <- sample.int(n - 1, 1); ii <- ((seq_len(n) - 1 + sh) %% n) + 1; NUL[, , k] <- sapply(sets, function(r) fitset(zv[ii], sv[ii], r)) }
  EZ <- (OBS - apply(NUL, 1:2, mean, na.rm = TRUE)) / apply(NUL, 1:2, sd, na.rm = TRUE)
  EP <- sapply(seq_along(sets), function(j) sapply(RD, function(rr) { nu <- NUL[rr, j, ]; nu <- nu[is.finite(nu)]; o <- OBS[rr, j]
    if (!is.finite(o) || length(nu) < 10) NA_real_ else (1 + sum(abs(nu - mean(nu)) >= abs(o - mean(nu)))) / (length(nu) + 1) })); dimnames(EP) <- dimnames(OBS)
  SUMM[[TR]] <- rbindlist(lapply(RD, function(rr) { lz <- EZ[rr, grep("^loo_", colnames(EZ))]
    held <- sum(abs(lz) > 1.96 & sign(lz) == sign(EZ[rr, "full"]), na.rm = TRUE); wk <- which.min(ifelse(is.finite(lz), sign(EZ[rr, "full"]) * lz, Inf))
    data.table(disorder = TR, twas = "GTEx TWAS30, scaled tissue z", readout = rr, loci = nrow(d), z_analytic = OBS[rr, "full"], P_analytic = 2 * pnorm(-abs(OBS[rr, "full"])),
               z_emp = EZ[rr, "full"], P_emp = EP[rr, "full"], z_emp_no16p = if (length(i16)) EZ[rr, "no16p"] else NA_real_, P_emp_no16p = if (length(i16)) EP[rr, "no16p"] else NA_real_,
               loo_held = held, loo_n = nrow(d), weakest = sub("^loo_", "", names(lz)[wk]), weakest_z = lz[[wk]], nrep = NREP) }))
  ## mean shift (analytic, LOO) -- same fits
  for (VAR in c("countA", "countU")) {   # countU = sensitivity (unsigned |z| > 3 count in place of count_A)
  ms <- function(rows) { o <- c(); for (al in c("DEL", "DUP")) { f <- feats((if (al == "DEL") GDEL else GDUP)[rows], zv, sv)
    r <- rema(d[[paste0("b_", al)]][rows], d[[paste0("se_", al)]][rows], cbind(int = 1, nall = f$nall, tw = f[[VAR]], noneH = f$noneH))
    o[paste0(c("int_z_", "nall_z_", "int_b_", "nall_b_"), al)] <- if (is.null(r)) NA else c(r$z[["int"]], r$z[["nall"]], r$b[["int"]], r$b[["nall"]] * 10) }; o }
  F <- ms(full); LO <- sapply(full, function(i) ms(setdiff(full, i))); colnames(LO) <- d$NAME
  comb <- function(zD, zU) ifelse(is.finite(zD) & is.finite(zU) & sign(zD) == sign(zU), sign(zD) * pmin(abs(zD), abs(zU)), 0)
  heldf <- function(l, f) sum(abs(l) > 1.96 & sign(l) == sign(f), na.rm = TRUE)
  MS[[paste(TR, VAR)]] <- rbindlist(lapply(c("int", "nall"), function(p) { zD <- F[[paste0(p, "_z_DEL")]]; zU <- F[[paste0(p, "_z_DUP")]]; zm <- comb(zD, zU)
    lD <- LO[paste0(p, "_z_DEL"), ]; lU <- LO[paste0(p, "_z_DUP"), ]
    data.table(disorder = TR, twas_term = VAR, part = if (p == "int") "count-independent (intercept, SD at zero genes)" else "count-dependent (total gene count slope, per 10 genes)",
               b_DEL = F[[paste0(p, "_b_DEL")]], z_DEL = zD, held_DEL = heldf(lD, zD), b_DUP = F[[paste0(p, "_b_DUP")]], z_DUP = zU, held_DUP = heldf(lU, zU),
               mean_shift_z = zm, mean_shift_held = if (zm == 0) 0L else heldf(comb(lD, lU), zm), loci = nrow(d),
               z_no16p_DEL = if (length(i16)) lD[[i16]] else NA_real_, z_no16p_DUP = if (length(i16)) lU[[i16]] else NA_real_) })) }
  ## ad hoc: presence/absence vs count for the constraint term, each dosage (ML fits)
  for (al in c("DEL", "DUP")) { G <- if (al == "DEL") GDEL else GDUP; f <- feats(G, zv, sv); f[, pres95 := as.numeric(n95 > 0)]
    y <- d[[paste0("b_", al)]]; s <- d[[paste0("se_", al)]]; base <- cbind(int = 1, nall = f$nall, countA = f$countA)
    forms <- list(`no constraint term` = base, `none_H count (ConGA)` = cbind(base, noneH = f$noneH), `top-5% s_het count` = cbind(base, n95 = f$n95),
                  `presence of a top-5% s_het gene` = cbind(base, pres95 = f$pres95), `none_H + presence` = cbind(base, noneH = f$noneH, pres95 = f$pres95))
    rs <- lapply(forms, function(X) rema(y, s, X, ml = TRUE)); a0 <- -2 * rs[[1]]$ll + 2 * rs[[1]]$k
    r0 <- rema(y, s, cbind(int = 1, countA = f$countA), ml = TRUE)   # gene content jointly: total count + none_H (2 df), collinear across loci (r ~ 0.9)
    CONT[[paste(TR, al)]] <- data.table(trait = TR, dosage = al, cor_nall_noneH = cor(f$nall[is.finite(y)], f$noneH[is.finite(y)]),
      LRT_2df = 2 * (rs[[2]]$ll - r0$ll), P_2df = pchisq(2 * (rs[[2]]$ll - r0$ll), 2, lower.tail = FALSE),
      z_noneH_alone = rema(y, s, cbind(int = 1, countA = f$countA, noneH = f$noneH))$z[["noneH"]], z_nall_alone = rema(y, s, cbind(int = 1, countA = f$countA, nall = f$nall))$z[["nall"]])
    for (nm in names(forms)) { r <- rs[[nm]]; tm <- setdiff(colnames(forms[[nm]]), colnames(base)); key <- if (length(tm)) tail(tm, 1) else NA
      lo <- if (is.na(key)) NA else sapply(which(is.finite(y)), function(i) { q <- rema(y[-i], s[-i], forms[[nm]][-i, , drop = FALSE], ml = TRUE); if (is.null(q)) NA else q$z[[key]] })
      HIT[[paste(TR, al, nm)]] <- data.table(trait = TR, dosage = al, form = nm, loci = sum(is.finite(y)), loci_with_term = if (is.na(key)) NA else sum(forms[[nm]][is.finite(y), key] > 0),
        term = if (is.na(key)) "" else key, b = if (is.na(key)) NA else r$b[[key]], z = if (is.na(key)) NA else r$z[[key]], P = if (is.na(key)) NA else 2 * pnorm(-abs(r$z[[key]])),
        noneH_z_alongside = if (identical(tm, c("noneH", "pres95"))) r$z[["noneH"]] else NA, logLik_ML = r$ll, dAIC_vs_none = (-2 * r$ll + 2 * r$k) - a0,
        loo_held = if (is.na(key)) NA else sum(abs(lo) > 1.96 & sign(lo) == sign(r$z[[key]]), na.rm = TRUE), tau2 = r$tau2) } }
  fD <- feats(GDEL, zv, sv); fU <- feats(GDUP, zv, sv)
  LOCI[[TR]] <- data.table(trait = TR, NAME = d$NAME, b_DEL = d$b_DEL, se_DEL = d$se_DEL, b_DUP = d$b_DUP, se_DUP = d$se_DUP,
                           nall_DEL = fD$nall, countA_DEL = fD$countA, noneH_DEL = fD$noneH, n95_DEL = fD$n95, nall_DUP = fU$nall, countA_DUP = fU$countA, noneH_DUP = fU$noneH, n95_DUP = fU$n95)
}
fwrite(rbindlist(SUMM), file.path(OUT, "growth_conga_axes_summary.csv")); fwrite(rbindlist(MS), file.path(OUT, "growth_conga_mean_shift.csv"))
fwrite(rbindlist(HIT), file.path(OUT, "growth_conga_presence_vs_count.csv")); fwrite(rbindlist(LOCI), file.path(OUT, "growth_conga_loci.csv")); fwrite(rbindlist(CONT), file.path(OUT, "growth_conga_gene_content.csv"))
options(width = 250); print(rbindlist(SUMM)[, .(disorder, readout, z_emp = round(z_emp, 2), P_emp = signif(P_emp, 2), no16p = round(z_emp_no16p, 2), loo = paste0(loo_held, "/", loo_n), weakest)])
print(rbindlist(MS)[, .(disorder, twas_term, part = substr(part, 1, 16), z_DEL = round(z_DEL, 2), held_DEL, z_DUP = round(z_DUP, 2), held_DUP, ms = round(mean_shift_z, 2), mean_shift_held)])
print(rbindlist(CONT)); print(rbindlist(HIT)[, .(trait, dosage, form, k = loci_with_term, b = round(b, 3), z = round(z, 2), dAIC = round(dAIC_vs_none, 1), loo_held)])
