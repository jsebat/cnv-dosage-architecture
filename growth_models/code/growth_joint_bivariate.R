# growth_joint_bivariate.R -- ONE locus-level joint model per trait (height, BMI) that carries the mirror, constraint and the single hit together
#   (Jon 2026-11-15 via the UKBB session: "a single locus-level joint model per trait that includes the single hit, so the schematic is not stitched
#   from separate fits").
#   Why bivariate: the height single hit is an offset FROM THE MIRROR LINE (conditional on the interval's deletion effect); in the per-arm ConGA fit it is
#   invisible (DUP top-5% presence z -0.67). Using b_DEL as a covariate (the offset model) attenuates the slope (b_DEL is noisy) and is a second model.
#   Model: both arms stacked; arm-specific fixed effects  intercept + total genes + TWAS |z|>3 count + none_H (s_het > 0.029) + union presence
#     (>= 1 gene with top-5% s_het OR ClinGen HI score 2-3), on BOTH arms (joint model quantifies all axes, Jon 2026-11-11);
#     locus random effects (u_DEL, u_DUP) ~ N(0, Sigma) with correlation rho (the mirror line, estimated with sampling error handled); REML.
#   TWAS term: unsigned |z|>3 count (primary for growth, report section G) and signed count_A (ConGA spec, sensitivity).
#   Readouts: mirror (z_DUP - z_DEL)/sqrt2 of the TWAS term; constraint z per arm; hit z per arm; total-gene slope per arm; rho.
#   Check model: DUP ~ intercept + b_DEL + total + TWAS count + none_H + union presence (univariate REML; the UKBB-session offset model with full covariates).
#   Null: genome shift of gene annotations (TWAS z, s_het, ClinGen HI flag rotated together), NREP replicates, empirical z;
#   LOO: removal of each locus FAMILY (cytoband group; UKBB-session caveat: single-interval removal misses overlapping families) and of each interval,
#   both matched in every null replicate. Gene assignment: DEL any refGene exon; DUP a full transcript inside. chrX dropped. Public data, runs locally.
#   Usage (from anywhere): Rscript growth_joint_bivariate.R <v2 dir> [nrep] [cores]
suppressMessages({ library(data.table); library(parallel) })
a <- commandArgs(TRUE); V <- a[1]; NREP <- if (length(a) >= 2) as.integer(a[2]) else 2000L; NC <- if (length(a) >= 3) as.integer(a[3]) else 12L
OUT <- file.path(V, "results/eight_traits/growth_conga/joint_bivariate_2026-11-15"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
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
hi_ez <- unique(c(sym2ez$entrez[match(hi_sym, sym2ez$symbol)], unique(ex[symbol %in% hi_sym]$entrez))); U[, hi := as.numeric(entrez %in% hi_ez)]
cat(sprintf("top-5%% s_het cut %.3f; ClinGen HI genes %d, mapped to refGene %d\n", Q95, length(hi_sym), sum(U$hi)))
CO <- unique(fread(file.path(M_, "pathway_analysis/Final/cnv_genes_annotated.csv"))[, .(NAME = CNV_ID, chr = Chromosome, START = CNV_start_GRCh38, END = CNV_end_GRCh38)])
S <- fread(file.path(V, "inputs/growth_input_signed.csv"))
TRAITS <- list(Height = list(lab = "Height", tw = "height"), BMI = list(lab = "BMI (TWAS30)", tw = "bmi"))
famof <- function(nm) { b <- sub("_.*", "", nm); b[grepl("^15q13", b)] <- "15q13"; b }

## bivariate REML. y, s: stacked effects and SEs; arm: 1 DEL / 2 DUP; loc: locus index; X: block design (columns already arm-specific)
bireml <- function(y, s, arm, loc, X, start = c(log(.1), log(.1), 0)) {
  keep <- sapply(colnames(X), function(cn) { rows <- arm == (if (grepl("_DEL$", cn)) 1 else 2); grepl("^int_", cn) || (sum(rows) > 2 && sd(X[rows, cn]) > 0) })
  X <- X[, keep, drop = FALSE]
  if (qr(X)$rank < ncol(X)) return(NULL)
  same <- outer(loc, loc, "=="); A1 <- outer(arm, arm, function(i, j) i == 1 & j == 1) & same; A2 <- outer(arm, arm, function(i, j) i == 2 & j == 2) & same
  A12 <- outer(arm, arm, "!=") & same; S2 <- diag(s^2)
  nll <- function(p) { t1 <- exp(p[1]); t2 <- exp(p[2]); r <- tanh(p[3])
    Vm <- S2 + t1^2 * A1 + t2^2 * A2 + r * t1 * t2 * A12; R <- tryCatch(chol(Vm), error = function(e) NULL); if (is.null(R)) return(1e10)
    Ri <- backsolve(R, diag(nrow(Vm))); Vi <- Ri %*% t(Ri); XtVi <- crossprod(X, Vi); B <- XtVi %*% X
    RB <- tryCatch(chol(B), error = function(e) NULL); if (is.null(RB)) return(1e10)
    b <- backsolve(RB, forwardsolve(t(RB), XtVi %*% y)); e <- y - X %*% b
    0.5 * (2 * sum(log(diag(R))) + 2 * sum(log(diag(RB))) + drop(crossprod(e, Vi %*% e))) }
  o <- optim(start, nll, method = "Nelder-Mead", control = list(maxit = 600, reltol = 1e-9)); p <- o$par
  t1 <- exp(p[1]); t2 <- exp(p[2]); r <- tanh(p[3]); Vm <- S2 + t1^2 * A1 + t2^2 * A2 + r * t1 * t2 * A12; Vi <- solve(Vm)
  B <- crossprod(X, Vi %*% X); Bi <- solve(B); b <- drop(Bi %*% crossprod(X, Vi %*% y))
  z <- setNames(rep(NA_real_, length(keep)), names(keep)); bb <- z; z[keep] <- b / sqrt(diag(Bi)); bb[keep] <- b
  list(z = z, b = bb, tau = c(t1, t2), rho = r, nll = o$value, par = p) }
unireml <- function(y, s, X) { keep <- apply(X, 2, function(v) sd(v) > 0); keep[1] <- TRUE; X <- X[, keep, drop = FALSE]; if (qr(X)$rank < ncol(X)) return(NULL)
  rl <- function(t2) { w <- 1 / (s^2 + t2); A <- t(X * w) %*% X; b <- solve(A, t(X * w) %*% y); e <- y - X %*% b; -0.5 * (sum(log(s^2 + t2)) + as.numeric(determinant(A)$modulus) + sum(w * e^2)) }
  t2 <- optimize(function(t) -rl(t), c(0, max(4 * var(y), 1e-4)))$minimum; w <- 1 / (s^2 + t2); A <- t(X * w) %*% X; b <- solve(A, t(X * w) %*% y)
  z <- setNames(rep(NA_real_, length(keep)), names(keep)); z[keep] <- b / sqrt(diag(solve(A))); z }

RD <- c("mirror", "tw_DEL", "tw_DUP", "con_DEL", "con_DUP", "hit_DEL", "hit_DUP", "nall_DEL", "nall_DUP", "rho", "chk_hit_DUP", "chk_con_DUP", "chk_tw_DUP", "chk_bDEL")
ALL <- list(); BETA <- list(); FEAT <- list()
for (TR in names(TRAITS)) {
  d <- S[disorder == TRAITS[[TR]]$lab, .(NAME, b_DEL = beta_DEL, se_DEL, b_DUP = beta_DUP, se_DUP)][CO, on = "NAME", nomatch = 0]
  d[!(is.finite(se_DEL) & se_DEL > 0), `:=`(b_DEL = NA, se_DEL = NA)]; d[!(is.finite(se_DUP) & se_DUP > 0), `:=`(b_DUP = NA, se_DUP = NA)]
  d <- d[!chr %in% c("chrX", "chrY")]; d[, fam := famof(NAME)]
  GDEL <- lapply(seq_len(nrow(d)), function(i) { x <- d[i]; U[entrez %in% unique(ex[chr == x$chr & end >= x$START & start <= x$END]$entrez)]$idx })
  GDUP <- lapply(seq_len(nrow(d)), function(i) { x <- d[i]; U[entrez %in% unique(TXS[chr == x$chr & ts >= x$START & te <= x$END]$entrez)]$idx })
  TW <- fread(file.path(C_, "twas30", paste0(TRAITS[[TR]]$tw, ".txt")))[is.finite(TWAS.Z), .(z = sum(TWAS.Z) / sqrt(.N)), by = GENE]
  TW[, entrez := sym2ez$entrez[match(GENE, sym2ez$symbol)]]; TW <- TW[!is.na(entrez)][!duplicated(entrez)]
  zv <- TW$z[match(U$entrez, TW$entrez)]; sv <- U$shet; hv <- U$hi; n <- length(zv)
  feats <- function(G, zz, ss, hh, tw) data.table(nall = lengths(G),
    tw = sapply(G, function(g) if (tw == "U") sum(abs(zz[g]) > 3, na.rm = TRUE) else sum(zz[g] > 3, na.rm = TRUE) - sum(zz[g] < -3, na.rm = TRUE)),
    con = sapply(G, function(g) sum(ss[g] > .029, na.rm = TRUE)), hit = sapply(G, function(g) as.numeric(any(ss[g] >= Q95, na.rm = TRUE) | any(hh[g] > 0))))
  ## stacked rows
  iD <- which(is.finite(d$b_DEL)); iU <- which(is.finite(d$b_DUP))
  fitset <- function(zz, ss, hh, rows, tw, start = NULL) {
    fD <- feats(GDEL, zz, ss, hh, tw); fU <- feats(GDUP, zz, ss, hh, tw); rD <- intersect(iD, rows); rU <- intersect(iU, rows)
    y <- c(d$b_DEL[rD], d$b_DUP[rU]); s <- c(d$se_DEL[rD], d$se_DUP[rU]); arm <- c(rep(1, length(rD)), rep(2, length(rU))); loc <- c(rD, rU)
    blk <- function(f, r, k) { m <- cbind(int = 1, nall = f$nall[r], tw = f$tw[r], con = f$con[r], hit = f$hit[r]); colnames(m) <- paste0(colnames(m), "_", k); m }
    X <- rbind(cbind(blk(fD, rD, "DEL"), blk(fU, rD, "DUP") * 0), cbind(blk(fD, rU, "DEL") * 0, blk(fU, rU, "DUP")))
    r <- bireml(y, s, arm, loc, X, start = if (is.null(start)) c(log(.1), log(.1), 0) else start)
    both <- intersect(rD, rU)   # check model: DUP conditional on the observed deletion effect
    ck <- unireml(d$b_DUP[both], d$se_DUP[both], cbind(int = 1, bDEL = d$b_DEL[both], nall = fU$nall[both], tw = fU$tw[both], con = fU$con[both], hit = fU$hit[both]))
    g <- function(v, k) if (is.null(v) || !k %in% names(v)) NA_real_ else v[[k]]
    zz_ <- if (is.null(r)) NULL else r$z
    list(v = c(mirror = (g(zz_, "tw_DUP") - g(zz_, "tw_DEL")) / sqrt(2), tw_DEL = g(zz_, "tw_DEL"), tw_DUP = g(zz_, "tw_DUP"), con_DEL = g(zz_, "con_DEL"), con_DUP = g(zz_, "con_DUP"),
              hit_DEL = g(zz_, "hit_DEL"), hit_DUP = g(zz_, "hit_DUP"), nall_DEL = g(zz_, "nall_DEL"), nall_DUP = g(zz_, "nall_DUP"), rho = if (is.null(r)) NA else r$rho,
              chk_hit_DUP = g(ck, "hit"), chk_con_DUP = g(ck, "con"), chk_tw_DUP = g(ck, "tw"), chk_bDEL = g(ck, "bDEL")), fit = r) }
  full <- seq_len(nrow(d)); fams <- unique(d$fam)
  sets <- c(list(full = full), setNames(lapply(fams, function(f) which(d$fam != f)), paste0("fam_", fams)), setNames(lapply(full, function(i) setdiff(full, i)), paste0("int_", d$NAME)))
  for (TWT in c("U", "A")) {
    t0 <- Sys.time(); F0 <- fitset(zv, sv, hv, full, TWT); st <- F0$fit$par
    OBS <- sapply(sets, function(r) fitset(zv, sv, hv, r, TWT, st)$v)
    BETA[[paste(TR, TWT)]] <- data.table(trait = TR, twas = TWT, term = names(F0$fit$b), b = F0$fit$b, z_analytic = F0$fit$z, tau_DEL = F0$fit$tau[1], tau_DUP = F0$fit$tau[2], rho = F0$fit$rho)
    cat(sprintf("%s %s: observed fits done (%.0f s); rho %.2f; running %d null replicates on %d cores\n", TR, TWT, as.numeric(Sys.time() - t0, units = "secs"), F0$fit$rho, NREP, NC)); flush.console()
    set.seed(20261115); SH <- sample.int(n - 1, NREP, replace = TRUE)
    cl <- makeCluster(NC); invisible(clusterEvalQ(cl, suppressMessages(library(data.table)))); clusterExport(cl, c("fitset", "feats", "bireml", "unireml", "zv", "sv", "hv", "n", "sets", "TWT", "st", "GDEL", "GDUP", "d", "iD", "iU", "Q95"), envir = environment())
    NL <- parLapply(cl, SH, function(sh) { ii <- ((seq_len(n) - 1 + sh) %% n) + 1; sapply(sets, function(r) fitset(zv[ii], sv[ii], hv[ii], r, TWT, st)$v) }); stopCluster(cl)
    NUL <- simplify2array(NL)   # RD x sets x NREP
    mu <- apply(NUL, 1:2, mean, na.rm = TRUE); sdv <- apply(NUL, 1:2, sd, na.rm = TRUE); EZ <- (OBS - mu) / sdv
    EP <- OBS * NA; for (i in seq_len(nrow(OBS))) for (j in seq_len(ncol(OBS))) { nu <- NUL[i, j, ]; nu <- nu[is.finite(nu)]; o <- OBS[i, j]
      if (is.finite(o) && length(nu) >= 10) EP[i, j] <- (1 + sum(abs(nu - mean(nu)) >= abs(o - mean(nu)))) / (length(nu) + 1) }
    held <- function(rr, pre) { lz <- EZ[rr, grep(pre, colnames(EZ))]; f <- EZ[rr, "full"]; c(n = sum(is.finite(lz)), held = sum(abs(lz) > 1.96 & sign(lz) == sign(f), na.rm = TRUE),
      minz = if (all(!is.finite(lz))) NA else unname(lz)[which.min(ifelse(is.finite(lz), sign(f) * lz, Inf))], weakest = sub(pre, "", names(lz)[which.min(ifelse(is.finite(lz), sign(f) * lz, Inf))])) }
    ALL[[paste(TR, TWT)]] <- rbindlist(lapply(setdiff(RD, "rho"), function(rr) { hf <- held(rr, "^fam_"); hi_ <- held(rr, "^int_")
      data.table(trait = TR, twas = if (TWT == "U") "unsigned |z|>3 count" else "signed count_A", readout = rr, z_analytic = OBS[rr, "full"], z_emp = EZ[rr, "full"], P_emp = EP[rr, "full"],
        fam_held = paste0(hf["held"], "/", hf["n"]), fam_min_z = as.numeric(hf["minz"]), fam_weakest = hf["weakest"],
        int_held = paste0(hi_["held"], "/", hi_["n"]), int_min_z = as.numeric(hi_["minz"]), int_weakest = hi_["weakest"], nrep = NREP) }))
    fwrite(data.table(set = colnames(EZ), t(EZ)), file.path(OUT, sprintf("ez_by_set_%s_%s.csv", TR, TWT)))
    cat(sprintf("%s %s done in %.1f min\n", TR, TWT, as.numeric(Sys.time() - t0, units = "mins"))); flush.console()
  }
  fD <- feats(GDEL, zv, sv, hv, "U"); fU <- feats(GDUP, zv, sv, hv, "U")
  FEAT[[TR]] <- data.table(trait = TR, NAME = d$NAME, family = d$fam, b_DEL = d$b_DEL, se_DEL = d$se_DEL, b_DUP = d$b_DUP, se_DUP = d$se_DUP,
    nall_DEL = fD$nall, twU_DEL = fD$tw, con_DEL = fD$con, hit_DEL = fD$hit, nall_DUP = fU$nall, twU_DUP = fU$tw, con_DUP = fU$con, hit_DUP = fU$hit)
}
R <- rbindlist(ALL); fwrite(R, file.path(OUT, "joint_bivariate_summary.csv")); fwrite(rbindlist(BETA), file.path(OUT, "joint_bivariate_coefficients.csv"))
fwrite(rbindlist(FEAT), file.path(OUT, "joint_bivariate_loci.csv"))
options(width = 250); print(R[, .(trait, twas = substr(twas, 1, 8), readout, za = round(z_analytic, 2), z = round(z_emp, 2), P = signif(P_emp, 2), fam_held, fam_min = round(fam_min_z, 2), fam_weakest, int_held, int_weakest)])
print(rbindlist(BETA)[, .(trait, twas, term, b = signif(b, 3), z = round(z_analytic, 2), rho = round(rho, 2))])
