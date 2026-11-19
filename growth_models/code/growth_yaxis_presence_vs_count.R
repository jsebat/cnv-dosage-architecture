# growth_yaxis_presence_vs_count.R -- does the height duplication offset fit presence/absence of a dosage-sensitive gene as well as or
#   better than a count? (Jon 2026-11-15). Same model as growth_yaxis_hit.R: duplication-arm REML meta-regression on the mirror line
#   (intercept + the deletion effect) plus one gene term; ML AIC. Gene sets: top-5% s_het, ClinGen HI (score 2-3), their union, and the
#   constrained-gene count (s_het > 0.029). Each set as presence (>= 1 gene) and as count; head-to-head fits with both; leave-one-out by
#   single interval and by locus family (all 16p11.2 intervals, all 22q11.2 intervals, both 7q11.23 intervals together).
#   Run from the v2 dir: Rscript code/growth_yaxis_presence_vs_count.R
suppressMessages(library(data.table)); options(width = 200)
L <- fread("results/eight_traits/growth_conga/growth_conga_loci.csv")
G <- fread("results/eight_traits/growth_conga/yaxis_monogenic_genes.csv")
rema0 <- function(y, s, X) { ok <- is.finite(y) & is.finite(s) & apply(is.finite(X), 1, all); y <- y[ok]; s <- s[ok]; X <- X[ok, , drop = FALSE]
  rl <- function(t2, ml) { w <- 1/(s^2+t2); A <- t(X*w)%*%X; b <- solve(A, t(X*w)%*%y); r <- y - X%*%b; -0.5*(sum(log(s^2+t2)) + (if (ml) 0 else determinant(A)$modulus) + sum(w*r^2)) }
  t2 <- optimize(function(t) -rl(t, FALSE), c(0, 4*var(y)))$minimum; w <- 1/(s^2+t2); A <- t(X*w)%*%X; b <- solve(A, t(X*w)%*%y)
  t2m <- optimize(function(t) -rl(t, TRUE), c(0, 4*var(y)))$minimum
  list(z = drop(b/sqrt(diag(solve(A)))), b = drop(b), aic = -2*(rl(t2m, TRUE) - .5*length(y)*log(2*pi)) + 2*(ncol(X)+1), n = length(y)) }
rema <- function(y, s, X) tryCatch(rema0(y, s, X), error = function(e) list(z = c(int = NA, del = NA, term = NA, pres = NA, count = NA), b = c(term = NA), aic = NA, n = length(y)))
ng <- function(v) sapply(strsplit(ifelse(is.na(v), "", v), " "), function(z) length(z[nzchar(z)]))
un <- function(a, b) mapply(function(p, q) length(unique(c(strsplit(ifelse(is.na(p), "", p), " ")[[1]], strsplit(ifelse(is.na(q), "", q), " ")[[1]])) |> setdiff("")), a, b)
OUT <- list(); LOO <- list()
for (tr in c("Height", "BMI")) {
  x <- L[trait == tr]; g <- G[trait == tr]; x <- g[, .(NAME, top5_genes, hi_clingen)][x, on = "NAME"]
  x[, `:=`(c_top5 = ng(top5_genes), c_hi = ng(hi_clingen), c_union = un(top5_genes, hi_clingen), c_con = noneH_DUP)]
  x[, `:=`(p_top5 = as.numeric(c_top5 > 0), p_hi = as.numeric(c_hi > 0), p_union = as.numeric(c_union > 0), p_con = as.numeric(c_con > 0))]
  x <- x[is.finite(b_DUP) & is.finite(b_DEL)]
  y <- x$b_DUP; s <- x$se_DUP; base <- cbind(int = 1, del = x$b_DEL); r0 <- rema(y, s, base)
  fam <- fifelse(grepl("^16p11", x$NAME), "16p11.2", fifelse(grepl("^22q11", x$NAME), "22q11.2", fifelse(grepl("^7q11.23", x$NAME), "7q11.23", x$NAME)))
  for (set in c("top5", "hi", "union", "con")) for (form in c("p", "c")) {
    v <- x[[paste0(form, "_", set)]]; if (sd(v) == 0) next
    r <- rema(y, s, cbind(base, term = v)); zf <- sapply(unique(fam), function(f) { k <- fam != f; if (sd(v[k]) == 0) NA else rema(y[k], s[k], cbind(base[k, ], term = v[k]))$z[["term"]] })
    zi <- sapply(seq_along(y), function(i) rema(y[-i], s[-i], cbind(base[-i, ], term = v[-i]))$z[["term"]])
    OUT[[paste(tr, set, form)]] <- data.table(trait = tr, set = set, form = ifelse(form == "p", "presence", "count"), intervals_with = sum(v > 0), loci = length(y),
      b = r$b[["term"]], z = r$z[["term"]], AIC = r$aic, dAIC_vs_mirror_only = r$aic - r0$aic,
      loo_interval_held = sum(abs(zi) > 1.96 & sign(zi) == sign(r$z[["term"]])), loo_family_min_abs_z = min(abs(zf), na.rm = TRUE),
      loo_family_weakest = names(zf)[which.min(abs(zf))])
  }
  for (set in c("top5", "hi", "union")) {   # head to head: presence and count of the same set in one model
    r <- rema(y, s, cbind(base, pres = x[[paste0("p_", set)]], count = x[[paste0("c_", set)]]))
    OUT[[paste(tr, set, "both")]] <- data.table(trait = tr, set = set, form = "presence + count", intervals_with = sum(x[[paste0("p_", set)]]), loci = length(y),
      b = NA, z = NA, AIC = r$aic, dAIC_vs_mirror_only = r$aic - r0$aic, z_presence_joint = r$z[["pres"]], z_count_joint = r$z[["count"]])
  }
  cat("\n==", tr, ": mirror line only AIC", round(r0$aic, 1), "\n")
}
R <- rbindlist(OUT, fill = TRUE); print(R[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 2) else v)])
fwrite(R, "results/eight_traits/growth_conga/yaxis_presence_vs_count.csv")
