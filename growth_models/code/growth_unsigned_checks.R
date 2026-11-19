# growth_unsigned_checks.R -- why does the height / BMI mirror track the UNSIGNED |z| > 3 count? (Jon 2026-11-14: "Yes run checks. Growth ... will be published
#   separately, so we can do a deeper dive with BMI and height"). Exploratory; the mirror definition stays the signed count (Jon 2026-11-14).
#   Three readings: (1) TWAS sign unreliable across GTEx tissues (magnitude survives the tissue sum, sign does not); (2) growth genes mostly push one way
#   (dosage prior independent of TWAS sign); (3) artefact: the unsigned count tracks density of well-imputed / broadly expressed genes.
#   Same model as growth_conga_axes.R: per dosage REML meta-regression of locus effect on intercept + total gene count + X (+ none_H), X per check; genomic-shift
#   null (all gene-level vectors rotated together), NREP replicates; every single-locus LOO (analytic z); mirror-style readout (z_DUP - z_DEL)/sqrt2 per term.
#   Checks:
#     A  reference: unsigned own count (countU), signed own count (countA)
#     B  negative controls (reading 3): unsigned count from the OTHER growth trait (same GTEx panel and pipeline), from SCZ and MD PECfull (brain), and density measures:
#        genes with any own-trait TWAS model, genes with >= 10 tissue models; then own countU jointly with each control
#     C  sign split (readings 1 vs 2): n(z > 3) and n(z < -3) as separate terms, jointly
#     D  sign reliability (reading 1): per gene, fraction of tissue z with the sign of the combined z; reliable = >= 0.8; signed count over reliable genes and
#        unsigned count over unreliable genes, jointly
#     E  descriptives: sign agreement and fraction positive among |z| > 3 genes, genome-wide vs in the loci
#   Usage: Rscript growth_unsigned_checks.R <v2 dir> [nrep]
suppressMessages(library(data.table))
a <- commandArgs(TRUE); V <- a[1]; NREP <- if (length(a) >= 2) as.integer(a[2]) else 1000L; set.seed(20261115)
OUT <- file.path(V, "results/eight_traits/growth_conga/unsigned_checks"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
M_ <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly"; C_ <- file.path(M_, "CNV_ASD_growth_analysis/concordance")
ex <- fread(file.path(V, "inputs/annot/hg38_refGene_20200708.exon.txt"), header = FALSE, col.names = c("chr", "start", "end", "tx", "symbol", "entrez"))
ex <- ex[chr %in% paste0("chr", 1:22)][, entrez := as.character(entrez)]
TXS <- ex[, .(ts = min(start), te = max(end)), by = .(tx, chr, entrez)]
U <- ex[, .(pos = min(start)), by = .(entrez, chr)][, .SD[1], by = entrez][, cn := as.integer(sub("chr", "", chr))][order(cn, pos)]; U[, idx := .I]
gc_ <- fread(file.path(V, "inputs/annot/gene_covariates.tsv"), colClasses = list(character = "entrez")); gc_[shet_imputed == 1, shet := NA]
U <- gc_[, .(entrez, shet)][U, on = "entrez"]; setorder(U, idx)
sym2ez <- rbind(gc_[, .(symbol, entrez)], unique(ex[, .(symbol, entrez)]))[!duplicated(symbol)]
gtex <- function(tw) { T <- fread(file.path(C_, "twas30", paste0(tw, ".txt")))[is.finite(TWAS.Z)]
  G <- T[, .(z = sum(TWAS.Z) / sqrt(.N), ntis = .N, agree = mean(sign(TWAS.Z) == sign(sum(TWAS.Z)))), by = GENE]
  G[, entrez := sym2ez$entrez[match(GENE, sym2ez$symbol)]]; G <- G[!is.na(entrez)][!duplicated(entrez)]; G[match(U$entrez, entrez)] }
pec <- function(tw) { T <- fread(file.path(V, "inputs/annot/twas", paste0(tw, ".PECfull.tsv")))[is.finite(z), .(entrez = as.character(entrez), z)][!duplicated(entrez)]; T$z[match(U$entrez, T$entrez)] }
GT <- list(Height = gtex("height"), BMI = gtex("bmi")); ZS <- pec("SCZ"); ZM <- pec("MD")
CO <- unique(fread(file.path(M_, "pathway_analysis/Final/cnv_genes_annotated.csv"))[, .(NAME = CNV_ID, chr = Chromosome, START = CNV_start_GRCh38, END = CNV_end_GRCh38)])
S <- fread(file.path(V, "inputs/growth_input_signed.csv"))
rema <- function(y, s, X) { ok <- is.finite(y) & is.finite(s); y <- y[ok]; s <- s[ok]; X <- X[ok, , drop = FALSE]
  keep <- apply(X, 2, function(v) sd(v) > 0); keep[1] <- TRUE; X <- X[, keep, drop = FALSE]
  if (qr(X)$rank < ncol(X) || nrow(X) <= ncol(X) + 1) return(NULL)
  rl <- function(t2) { w <- 1 / (s^2 + t2); XtW <- t(X * w); A <- XtW %*% X; b <- solve(A, XtW %*% y); r <- y - X %*% b
    -0.5 * (sum(log(s^2 + t2)) + as.numeric(determinant(A)$modulus) + sum(w * r^2)) }
  t2 <- optimize(function(t) -rl(t), c(0, max(4 * var(y), 1e-4)))$minimum
  w <- 1 / (s^2 + t2); A <- t(X * w) %*% X; b <- solve(A, t(X * w) %*% y)
  z <- setNames(rep(NA_real_, length(keep)), names(keep)); z[keep] <- b / sqrt(diag(solve(A))); z }
RES <- list(); DESC <- list()
for (TR in c("Height", "BMI")) {
  OT <- setdiff(c("Height", "BMI"), TR)
  d <- S[disorder == c(Height = "Height", BMI = "BMI (TWAS30)")[[TR]], .(NAME, b_DEL = beta_DEL, se_DEL, b_DUP = beta_DUP, se_DUP)][CO, on = "NAME", nomatch = 0][!chr %in% c("chrX", "chrY")]
  d[!(is.finite(se_DEL) & se_DEL > 0), `:=`(b_DEL = NA, se_DEL = NA)]
  GDEL <- lapply(seq_len(nrow(d)), function(i) { x <- d[i]; U[entrez %in% unique(ex[chr == x$chr & end >= x$START & start <= x$END]$entrez)]$idx })
  GDUP <- lapply(seq_len(nrow(d)), function(i) { x <- d[i]; U[entrez %in% unique(TXS[chr == x$chr & ts >= x$START & te <= x$END]$entrez)]$idx })
  g <- GT[[TR]]; o <- GT[[OT]]
  # gene-level matrix (universe order); every column rotated together in the null
  GM <- cbind(z = g$z, ntis = g$ntis, agree = g$agree, zo = o$z, zs = ZS, zm = ZM, shet = U$shet)
  ## locus features from a gene matrix
  feats <- function(M, G) { f <- function(fun) sapply(G, function(i) fun(M[i, , drop = FALSE]))
    data.table(nall = lengths(G), noneH = f(function(m) sum(m[, "shet"] > .029, na.rm = TRUE)),
      countU = f(function(m) sum(abs(m[, "z"]) > 3, na.rm = TRUE)), countA = f(function(m) sum(m[, "z"] > 3, na.rm = TRUE) - sum(m[, "z"] < -3, na.rm = TRUE)),
      npos = f(function(m) sum(m[, "z"] > 3, na.rm = TRUE)), nneg = f(function(m) sum(m[, "z"] < -3, na.rm = TRUE)),
      ctrl_other = f(function(m) sum(abs(m[, "zo"]) > 3, na.rm = TRUE)), ctrl_scz = f(function(m) sum(abs(m[, "zs"]) > 3, na.rm = TRUE)), ctrl_md = f(function(m) sum(abs(m[, "zm"]) > 3, na.rm = TRUE)),
      dens_any = f(function(m) sum(is.finite(m[, "z"]))), dens_10 = f(function(m) sum(m[, "ntis"] >= 10, na.rm = TRUE)),
      rel_signed = f(function(m) { k <- abs(m[, "z"]) > 3 & m[, "agree"] >= .8; sum(k & m[, "z"] > 0, na.rm = TRUE) - sum(k & m[, "z"] < 0, na.rm = TRUE) }),
      unrel_unsigned = f(function(m) sum(abs(m[, "z"]) > 3 & m[, "agree"] < .8, na.rm = TRUE)), rel_unsigned = f(function(m) sum(abs(m[, "z"]) > 3 & m[, "agree"] >= .8, na.rm = TRUE))) }
  MODELS <- list(A_unsigned = "countU", A_signed = "countA", B_other_trait = "ctrl_other", B_scz = "ctrl_scz", B_md = "ctrl_md", B_dens_any = "dens_any", B_dens_10 = "dens_10",
                 B_own_plus_other = c("countU", "ctrl_other"), B_own_plus_scz = c("countU", "ctrl_scz"), B_own_plus_md = c("countU", "ctrl_md"), B_own_plus_dens10 = c("countU", "dens_10"),
                 C_pos_neg = c("npos", "nneg"), D_reliable_signed_unreliable_unsigned = c("rel_signed", "unrel_unsigned"), D_reliable_unsigned = "rel_unsigned")
  fitall <- function(M, rows) { FD <- feats(M, GDEL[rows]); FU <- feats(M, GDUP[rows]); out <- c()
    for (mn in names(MODELS)) { tm <- MODELS[[mn]]
      zD <- rema(d$b_DEL[rows], d$se_DEL[rows], as.matrix(cbind(int = 1, FD[, c("nall", tm, "noneH"), with = FALSE])))
      zU <- rema(d$b_DUP[rows], d$se_DUP[rows], as.matrix(cbind(int = 1, FU[, c("nall", tm, "noneH"), with = FALSE])))
      for (t in tm) { a1 <- if (is.null(zD)) NA else zD[[t]]; a2 <- if (is.null(zU)) NA else zU[[t]]
        out[paste(mn, t, c("DEL", "DUP", "mirror"), sep = "|")] <- c(a1, a2, (a2 - a1) / sqrt(2)) } }
    out }
  full <- seq_len(nrow(d)); n <- nrow(GM)
  OBS <- fitall(GM, full); LOO <- sapply(full, function(i) fitall(GM, setdiff(full, i)))
  NUL <- sapply(seq_len(NREP), function(k) { sh <- sample.int(n - 1, 1); ii <- ((seq_len(n) - 1 + sh) %% n) + 1; fitall(GM[ii, , drop = FALSE], full) })
  mu <- rowMeans(NUL, na.rm = TRUE); sdv <- apply(NUL, 1, sd, na.rm = TRUE)
  ep <- sapply(seq_along(OBS), function(j) { nu <- NUL[j, ]; nu <- nu[is.finite(nu)]; if (!is.finite(OBS[j]) || length(nu) < 10) NA else (1 + sum(abs(nu - mean(nu)) >= abs(OBS[j] - mean(nu)))) / (length(nu) + 1) })
  kk <- do.call(rbind, strsplit(names(OBS), "|", fixed = TRUE))
  RES[[TR]] <- data.table(trait = TR, model = kk[, 1], term = kk[, 2], readout = kk[, 3], z_analytic = OBS, z_emp = (OBS - mu) / sdv, P_emp = ep,
    loo_held = rowSums(abs(LOO) > 1.96 & sign(LOO) == sign(OBS), na.rm = TRUE), loci = nrow(d), nrep = NREP,
    z_analytic_no16p = LOO[, which(d$NAME == "16p11.2_BP4_BP5")])
  ## descriptives
  inl <- unique(unlist(c(GDEL, GDUP))); big <- which(abs(g$z) > 3)
  DESC[[TR]] <- rbind(data.table(trait = TR, set = "genome-wide", genes_absz_gt3 = length(big), frac_positive = mean(g$z[big] > 0), median_sign_agreement = median(g$agree[big]),
                                 frac_agreement_ge_0.8 = mean(g$agree[big] >= .8), median_tissue_models = median(g$ntis[big])),
                      data.table(trait = TR, set = "in the CNV loci", genes_absz_gt3 = length(intersect(big, inl)), frac_positive = mean(g$z[intersect(big, inl)] > 0),
                                 median_sign_agreement = median(g$agree[intersect(big, inl)]), frac_agreement_ge_0.8 = mean(g$agree[intersect(big, inl)] >= .8),
                                 median_tissue_models = median(g$ntis[intersect(big, inl)])),
                      data.table(trait = TR, set = "all genes with a TWAS model", genes_absz_gt3 = sum(is.finite(g$z)), frac_positive = mean(g$z > 0, na.rm = TRUE),
                                 median_sign_agreement = median(g$agree, na.rm = TRUE), frac_agreement_ge_0.8 = mean(g$agree >= .8, na.rm = TRUE), median_tissue_models = median(g$ntis, na.rm = TRUE)))
  # per-locus features for the report (public data)
  fwrite(cbind(d[, .(NAME, b_DEL, se_DEL, b_DUP, se_DUP)], setnames(feats(GM, GDUP), function(x) paste0(x, "_DUP"))), file.path(OUT, sprintf("loci_features_%s.csv", TR)))
}
R <- rbindlist(RES); fwrite(R, file.path(OUT, "unsigned_checks.csv")); fwrite(rbindlist(DESC), file.path(OUT, "unsigned_checks_descriptives.csv"))
options(width = 250); print(R[, .(trait, model, term, readout, z_an = round(z_analytic, 2), z_emp = round(z_emp, 2), P_emp = signif(P_emp, 2), loo = paste0(loo_held, "/", loci))], nrows = 300)
print(rbindlist(DESC))
