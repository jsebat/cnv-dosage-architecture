## Independent re-implementation of Omar's allele-level association (Table S9) for SCZ at 4 alleles.
## Allele calls rebuilt from raw PLINK CNVs with his overlap rules:
##   cytoband loci (e.g. 3q29): bedtools -r -f 0.5 (50% reciprocal);  BP loci: -F 0.5 (CNV covers >=50% of locus);
##   compound alleles per CNV_genotype_combine.py (22q11.21 A-B & B-C & C-D -> A-D; 15q11.2 BP1-BP2 dropped if BP1-BP3).
## Per platform: glm(aff ~ allele + C1..C10) LRT (noFam) / clogit(+strata(FID)) (Fam); N-weighted Stouffer meta.
suppressMessages(library(survival))
D <- "/expanse/projects/sebat1/oshanta/home/IndividualDisorder/Freeze_2024/inputs/data_10perc/ALL/SCZ/platform"
W <- "/expanse/projects/sebat1/oshanta/home/IndividualDisorder/Freeze_2024/inputs/results/r2_burden_noSex/wdir/platform/SCZ"
plat_noFam <- c("Affy5","Affy500K","Affy6","Axiom","Illum300K","Illum550","Illum610","IllumOmniExpress","Psych")
plat_Fam <- c("Affy6","Illum610","Psych")
pcs <- paste0("C",1:10)
ov <- function(s1,e1,s2,e2) pmax(0, pmin(e1,e2) - pmax(s1,s2))
recip <- function(c, s, e) { o <- ov(c$BP1,c$BP2,s,e); o >= 0.5*(c$BP2-c$BP1) & o >= 0.5*(e-s) }
coverF <- function(c, s, e) ov(c$BP1,c$BP2,s,e) >= 0.5*(e-s)
call_alleles <- function(cnv, ids) {
  del <- cnv[cnv$TYPE==1,]; dup <- cnv[cnv$TYPE==3,]
  has <- function(x, sel) ids %in% x$matchID[sel]
  g <- data.frame(matchID=ids)
  d3 <- del[del$CHR==3,]; g$X_3q29_DEL <- as.integer(has(d3, recip(d3,196011534,197617305)))
  d15 <- del[del$CHR==15,]
  b12 <- has(d15, coverF(d15,22778537,23067754)); b13 <- has(d15, coverF(d15,22778537,28145193))
  g$X_15q11.2_BP1.BP2_DEL <- as.integer(b12 & !b13)
  u16 <- dup[dup$CHR==16,]; g$X_16p11.2_BP4.BP5_DUP <- as.integer(has(u16, coverF(u16,29639519,30189452)))
  d22 <- del[del$CHR==22,]
  AB <- has(d22, coverF(d22,19049819,20250000)); BC <- has(d22, coverF(d22,20400000,20665000))
  CD <- has(d22, coverF(d22,20750000,21112437))
  BD <- BC & CD; AD <- AB & BD
  g$X_22q11.21_A.D_DEL <- as.integer(AD)
  g
}
alleles <- c("X_3q29_DEL","X_15q11.2_BP1.BP2_DEL","X_16p11.2_BP4.BP5_DUP","X_22q11.21_A.D_DEL")
res <- list()
for (fam in c(FALSE,TRUE)) for (pl in if (fam) plat_Fam else plat_noFam) {
  pre <- sprintf("%s/all_platforms_10kb_10probes_%s_SCZ_ALL_%s", D, if (fam) "Fams" else "noFams", pl)
  phe <- read.table(paste0(pre,".phe"), h=T, stringsAsFactors=F); phe$aff <- phe$AFF - 1
  phe$matchID <- paste0(phe$FID, ":", phe$IID)
  cnv <- read.table(paste0(pre,".cnv"), h=T, stringsAsFactors=F); cnv$matchID <- paste0(cnv$FID, ":", cnv$IID)
  g <- call_alleles(cnv, phe$matchID)
  om <- tryCatch(read.table(sprintf("%s/%s_%s/genotype_allele.indiv", W, if (fam) "Fam" else "noFam", pl), h=T, check.names=F, stringsAsFactors=F), error=function(e) NULL)
  dat <- merge(phe, g, by="matchID")
  Nca <- sum(phe$AFF==2); Nco <- sum(phe$AFF==1)
  for (a in alleles) {
    n_omar <- if (!is.null(om) && a %in% names(om)) { m <- merge(phe[,c("matchID","aff")], om[,c("matchID",a)], by="matchID"); sum(m[[a]]>0) } else NA
    agree <- if (!is.null(om) && a %in% names(om)) { m <- merge(g[,c("matchID",a)], om[,c("matchID",a)], by="matchID"); mean(m[[2]]==(m[[3]]>0)) } else NA
    n <- sum(dat[[a]]); cas <- sum(dat[[a]]*dat$aff)
    b <- 0; se <- NA; p <- 1
    if (n > 0) {
      if (!fam) {
        m0 <- glm(as.formula(paste("aff ~", paste(pcs,collapse="+"))), data=dat, family=binomial)
        m1 <- glm(as.formula(paste("aff ~", a, "+", paste(pcs,collapse="+"))), data=dat, family=binomial)
        an <- anova(m0, m1, test="Chisq"); p <- an[2,5]; b <- coef(m1)[2]; se <- summary(m1)$coefficients[2,2]
      } else {
        m0 <- clogit(as.formula(paste("aff ~", paste(pcs,collapse="+"), "+ strata(FID)")), data=dat)
        m1 <- clogit(as.formula(paste("aff ~", a, "+", paste(pcs,collapse="+"), "+ strata(FID)")), data=dat)
        an <- anova(m0, m1); p <- an[2, ncol(an)]; b <- coef(m1)[1]; se <- sqrt(vcov(m1)[1,1])
      }
      if (is.na(b)) b <- 0
      if (is.na(p)) p <- 1
    }
    res[[length(res)+1]] <- data.frame(allele=a, platform=pl, fam=fam, Ncas=Nca, Ncon=Nco, n=n, cases=cas, ctrl=n-cas,
                                       n_omar_matrix=n_omar, agree_with_omar=agree, beta=unname(b), se=se, p=p)
  }
}
res <- do.call(rbind, res)
res$z <- sign(res$beta) * -qnorm(res$p/2); res$w <- sqrt(4/(1/res$Ncas + 1/res$Ncon))
meta <- do.call(rbind, lapply(split(res, res$allele), function(d) data.frame(allele=d$allele[1],
  cases=sum(d$cases), ctrl=sum(d$ctrl), Z=sum(d$z*d$w)/sqrt(sum(d$w^2)))))
meta$P <- 2*pnorm(-abs(meta$Z))
out <- "/home/jsebat/ukbb_anthro/replication/out"
write.table(res, file.path(out,"s9_per_platform.tsv"), sep="\t", quote=F, row.names=F)
write.table(meta, file.path(out,"s9_meta.tsv"), sep="\t", quote=F, row.names=F)
print(res[,c("allele","platform","fam","n","n_omar_matrix","agree_with_omar","cases","ctrl","p")])
print(meta, digits=12)
