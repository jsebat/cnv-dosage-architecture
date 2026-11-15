## Independent re-implementation of Omar Shanta's probe-level CNV GWAS (per-platform
## logistic / clogit LRT, then N-weighted Stouffer meta) at 4 SCZ lead probes.
## Compares to Omar's META_N_noSex (current, Nov 2025) and META (= Table S3, Nov 2024).
suppressMessages(library(survival))
D <- "/expanse/projects/sebat1/oshanta/home/IndividualDisorder/Freeze_2024/inputs/data_10perc/ALL/SCZ/platform"
probes <- data.frame(name=c("22q11.21","16p11.2prox","3q29","15q11.2","22q11.21_S3"),
  type=c(1,3,1,1,1), chr=c(22,16,3,15,22),
  bp=c(20069935,29897071,197136643,22990339,20600473))
plat_noFam <- c("Affy5","Affy500K","Affy6","Axiom","Illum300K","Illum550","Illum610","IllumOmniExpress","Psych")
plat_Fam <- c("Affy6","Illum610","Psych")
pcs <- paste0("C",1:10)
fit_one <- function(phe, cnv, p, fam, sex) {
  cc <- cnv[cnv$TYPE==p$type & cnv$CHR==p$chr & cnv$BP1<=p$bp & cnv$BP2>=p$bp, ]
  phe$CNV_at_SNP <- as.integer(phe$IID %in% cc$IID)
  nC <- sum(phe$CNV_at_SNP); cas <- sum(phe$CNV_at_SNP*phe$aff)
  if (nC==0) return(c(n=0,cases=0,ctrl=0,beta=0,se=0,p=1))
  cov <- c(if (sex) "SEX", pcs)
  if (!fam) {
    f0 <- as.formula(paste("aff ~", paste(cov, collapse="+")))
    f1 <- as.formula(paste("aff ~ CNV_at_SNP +", paste(cov, collapse="+")))
    m0 <- glm(f0, data=phe, family=binomial); m1 <- glm(f1, data=phe, family=binomial)
    pv <- anova(m0, m1, test="Chisq")[2, ncol(anova(m0,m1,test="Chisq"))]
    b <- coef(m1)["CNV_at_SNP"]; se <- summary(m1)$coefficients[2,2]
  } else {
    f0 <- as.formula(paste("aff ~", paste(c(cov,"strata(FID)"), collapse="+")))
    f1 <- as.formula(paste("aff ~ CNV_at_SNP +", paste(c(cov,"strata(FID)"), collapse="+")))
    m0 <- clogit(f0, data=phe); m1 <- clogit(f1, data=phe)
    a <- anova(m0, m1); pv <- a[2, ncol(a)]
    b <- coef(m1)["CNV_at_SNP"]; se <- sqrt(vcov(m1)[1,1])
  }
  if (is.na(b)) b <- 0
  c(n=nC, cases=cas, ctrl=nC-cas, beta=unname(b), se=se, p=pv)
}
run_version <- function(sex) {
  res <- list()
  for (fam in c(FALSE, TRUE)) for (pl in if (fam) plat_Fam else plat_noFam) {
    pre <- sprintf("%s/all_platforms_10kb_10probes_%s_SCZ_ALL_%s", D, if (fam) "Fams" else "noFams", pl)
    phe <- read.table(paste0(pre,".phe"), h=T, stringsAsFactors=F); phe$aff <- phe$AFF - 1
    cnv <- read.table(paste0(pre,".cnv"), h=T, stringsAsFactors=F)
    Nca <- sum(phe$AFF==2); Nco <- sum(phe$AFF==1)
    for (k in seq_len(nrow(probes))) {
      r <- fit_one(phe, cnv, probes[k,], fam, sex)
      res[[length(res)+1]] <- data.frame(probe=probes$name[k], platform=pl, fam=fam, Ncas=Nca, Ncon=Nco, t(r))
    }
  }
  res <- do.call(rbind, res)
  res$z <- sign(res$beta) * -qnorm(res$p/2)
  res$w <- sqrt(4/(1/res$Ncas + 1/res$Ncon))
  meta <- do.call(rbind, lapply(split(res, res$probe), function(d)
    data.frame(probe=d$probe[1], cases=sum(d$cases), ctrl=sum(d$ctrl),
               Z=sum(d$z*d$w)/sqrt(sum(d$w^2)))))
  meta$P <- 2*pnorm(-abs(meta$Z)); meta$version <- if (sex) "SEX+PCs (=META, Table S3)" else "PCs only (=META_N_noSex)"
  list(per=res, meta=meta)
}
out <- "/home/jsebat/ukbb_anthro/replication/out"; dir.create(out, showWarnings=F, recursive=T)
for (sex in c(FALSE, TRUE)) {
  r <- run_version(sex); tag <- if (sex) "sex" else "noSex"
  write.table(r$per, sprintf("%s/per_platform_%s.tsv", out, tag), sep="\t", quote=F, row.names=F)
  write.table(r$meta, sprintf("%s/meta_%s.tsv", out, tag), sep="\t", quote=F, row.names=F)
  print(r$meta, digits=10)
}
