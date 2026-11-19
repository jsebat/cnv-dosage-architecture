# growth_yaxis_hit.R -- height asymmetric CNVs shift along the duplication (y) axis (Jon 2026-11-14). Duplication-arm REML fits with a mirror term (the deletion effect, i.e.
#   the mirror line; or the ConGA count) plus presence / count of top-5% s_het genes in the duplicated interval. ML AIC; LOO for presence. Run from the v2 dir.
suppressMessages(library(data.table)); options(width=200)
L <- fread("results/eight_traits/growth_conga/growth_conga_loci.csv")
F <- rbind(fread("results/eight_traits/growth_conga/unsigned_checks/loci_features_Height.csv")[, trait := "Height"], fread("results/eight_traits/growth_conga/unsigned_checks/loci_features_BMI.csv")[, trait := "BMI"], fill = TRUE)
L <- F[, .(trait, NAME, countU_DUP)][L, on = .(trait, NAME)]
rema <- function(y, s, X) { ok <- is.finite(y) & is.finite(s) & apply(is.finite(X),1,all); y <- y[ok]; s <- s[ok]; X <- X[ok, , drop = FALSE]
  rl <- function(t2, ml) { w <- 1/(s^2+t2); A <- t(X*w)%*%X; b <- solve(A, t(X*w)%*%y); r <- y - X%*%b; -0.5*(sum(log(s^2+t2)) + (if (ml) 0 else determinant(A)$modulus) + sum(w*r^2)) }
  t2 <- optimize(function(t) -rl(t, FALSE), c(0, 4*var(y)))$minimum; w <- 1/(s^2+t2); A <- t(X*w)%*%X; b <- solve(A, t(X*w)%*%y)
  t2m <- optimize(function(t) -rl(t, TRUE), c(0, 4*var(y)))$minimum
  list(z = round(drop(b/sqrt(diag(solve(A)))), 2), b = round(drop(b), 3), aic = -2*(rl(t2m, TRUE) - .5*length(y)*log(2*pi)) + 2*(ncol(X)+1), n = length(y)) }
OUT <- list(); LOOR <- list()
for (tr in c("Height", "BMI")) { x <- L[trait == tr]; x[, pres := as.numeric(n95_DUP > 0)]; y <- x$b_DUP; s <- x$se_DUP
  M <- list(
    `ConGA DUP fit + presence` = cbind(int=1, nall=x$nall_DUP, countA=x$countA_DUP, noneH=x$noneH_DUP, pres=x$pres),
    `ConGA DUP fit with unsigned count + presence` = cbind(int=1, nall=x$nall_DUP, countU=x$countU_DUP, noneH=x$noneH_DUP, pres=x$pres),
    `mirror line from the deletion effect + presence` = cbind(int=1, del=x$b_DEL, pres=x$pres),
    `mirror line + presence + top-5% count` = cbind(int=1, del=x$b_DEL, pres=x$pres, n95=x$n95_DUP),
    `mirror line + top-5% count` = cbind(int=1, del=x$b_DEL, n95=x$n95_DUP),
    `mirror line + none_H count` = cbind(int=1, del=x$b_DEL, noneH=x$noneH_DUP),
    `mirror line only` = cbind(int=1, del=x$b_DEL))
  cat("\n==", tr, "\n"); for (m in names(M)) { r <- rema(y, s, M[[m]]); { cat(sprintf("%-62s AIC %6.1f  z: %s
", m, r$aic, paste(names(r$z), r$z, collapse = " "))); OUT[[paste(tr, m)]] <- data.table(trait = tr, model = m, AIC = r$aic, term = names(r$z), b = r$b, z = r$z, loci = r$n) } }
  # LOO for presence in mirror-line model
  lo <- sapply(seq_len(nrow(x)), function(i) rema(y[-i], s[-i], cbind(int=1, del=x$b_DEL[-i], pres=x$pres[-i]))$z[["pres"]]); LOOR[[tr]] <- data.table(trait = tr, held = sum(abs(lo) > 1.96 & sign(lo) == sign(rema(y, s, cbind(int=1, del=x$b_DEL, pres=x$pres))$z[["pres"]])), loci = length(lo), weakest_abs_z = min(abs(lo))); cat("LOO presence (mirror line): held", sum(abs(lo) > 1.96 & sign(lo) == sign(rema(y, s, cbind(int=1, del=x$b_DEL, pres=x$pres))$z[["pres"]])), "/", length(lo), " weakest z", round(min(abs(lo)),2), "\n") }
fwrite(rbindlist(OUT), "results/eight_traits/growth_conga/yaxis_hit.csv"); fwrite(rbindlist(LOOR), "results/eight_traits/growth_conga/yaxis_hit_loo.csv")
