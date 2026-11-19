# growth_yaxis_presence_only.R -- is height's duplication-axis (y) shift restricted to intervals carrying a top-5% s_het gene? (Jon 2026-11-14)
#   Duplication offset from the mirror line = beta_DUP - (-rho * beta_DEL), rho = 0.351 (fixed line through the origin, as in Figure 3); REML random-effects mean of the
#   offset in carriers and non-carriers; then within non-carriers, does the offset track other constraint measures (none_H, top-10% / top-20% s_het counts)?
#   Also the free-slope version (DUP ~ DEL with group intercepts). Run from the v2 dir.
suppressMessages(library(data.table)); options(width = 200)
L <- fread("results/eight_traits/growth_conga/growth_conga_loci.csv")[trait == "Height"]
ex <- fread("inputs/annot/hg38_refGene_20200708.exon.txt", header = FALSE, col.names = c("chr", "start", "end", "tx", "symbol", "entrez"))[chr %in% paste0("chr", 1:22)][, entrez := as.character(entrez)]
gc_ <- fread("inputs/annot/gene_covariates.tsv", colClasses = list(character = "entrez")); gc_[shet_imputed == 1, shet := NA]
q <- quantile(gc_$shet, c(.8, .9, .95), na.rm = TRUE)
CO <- unique(fread("C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/pathway_analysis/Final/cnv_genes_annotated.csv")[, .(NAME = CNV_ID, chr = Chromosome, START = CNV_start_GRCh38, END = CNV_end_GRCh38)])
TXS <- ex[, .(ts = min(start), te = max(end)), by = .(tx, chr, entrez)]
L <- CO[L, on = "NAME"]
sh <- function(i) { x <- L[i]; gc_[entrez %in% unique(TXS[chr == x$chr & ts >= x$START & te <= x$END]$entrez)]$shet }
L[, `:=`(n80 = sapply(seq_len(.N), function(i) sum(sh(i) >= q[1], na.rm = TRUE)), n90 = sapply(seq_len(.N), function(i) sum(sh(i) >= q[2], na.rm = TRUE)),
         maxshet = sapply(seq_len(.N), function(i) max(c(0, sh(i)), na.rm = TRUE)))]
L[, carrier := n95_DUP > 0]; L[, off := b_DUP + .351 * b_DEL]; L[, se_off := sqrt(se_DUP^2 + (.351 * se_DEL)^2)]
rema <- function(y, s, X) { rl <- function(t2) { w <- 1/(s^2+t2); A <- t(X*w)%*%X; b <- solve(A, t(X*w)%*%y); r <- y - X%*%b; -0.5*(sum(log(s^2+t2)) + determinant(A)$modulus + sum(w*r^2)) }
  t2 <- optimize(function(t) -rl(t), c(0, 4*var(y)))$minimum; w <- 1/(s^2+t2); A <- t(X*w)%*%X; b <- solve(A, t(X*w)%*%y); se <- sqrt(diag(solve(A)))
  data.table(term = colnames(X), b = drop(b), se = se, z = drop(b)/se, P = 2*pnorm(-abs(drop(b)/se)), tau = sqrt(t2), n = length(y)) }
R <- list()
for (g in c(TRUE, FALSE)) { x <- L[carrier == g]; R[[paste("mean offset", g)]] <- rema(x$off, x$se_off, cbind(mean = rep(1, nrow(x))))[, group := ifelse(g, "carriers", "non-carriers")] }
x <- L[carrier == FALSE]
for (v in c("noneH_DUP", "n90", "n80", "nall_DUP", "maxshet")) R[[v]] <- rema(x$off, x$se_off, cbind(int = 1, X = x[[v]]))[, `:=`(group = "non-carriers", measure = v)]
R[["free"]] <- rema(L$b_DUP, L$se_DUP, cbind(noncarrier = as.numeric(!L$carrier), carrier = as.numeric(L$carrier), del = L$b_DEL))[, group := "all, free slope"]
RR <- rbindlist(R, fill = TRUE); fwrite(RR, "results/eight_traits/growth_conga/yaxis_presence_only.csv")
print(RR[, .(group, measure, term, b = round(b, 3), se = round(se, 3), z = round(z, 2), P = signif(P, 2), n)])
print(L[order(carrier, off), .(NAME, carrier, n95 = n95_DUP, n90, off = round(off, 2), se = round(se_off, 2), maxshet = round(maxshet, 2))])
