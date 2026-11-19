# Rebuild the height/BMI per-interval predictors with the corrected gene -> s_het mapping (2026-11-09) and add the dosage-sensitive gene counts.
# hgnc_map's newer GRCh38 Ensembl ids for some genes (all of 17q12, HNF1B, LHX1, RBM8A, ACACA, TXNIP, ...) are not used by GeneBayes, so Ensembl-
# only matching left them without s_het and out of M, K and the tiers. Fallback: HGNC id (current, previous and alias symbols).
# Recomputed (TWAS30, mean z across tissues per gene; s_het cut 0.029 as in growth_two_axis_input.py):
#   M = sum z over less-constrained TWAS genes (s_het < 0.029); M_abs = sum |z| over them; K_twas = sum |z| over constrained TWAS genes
#   (s_het >= 0.029); K_n = their number; S_con = sum z over them; S_shet = sum z over all TWAS genes with s_het.
# Added (all genes in the interval, no TWAS; GeneBayes s_het percentile genome-wide):
#   D67, D80, D90, D95 = number of genes above the 67th/80th/90th/95th percentile; DDG2P_dosage = monoallelic DDG2P genes (definitive/strong) with
#   a loss-of-function mechanism.
# The old-mapping values are recomputed first and checked against the current inputs; old files are kept as *_v1_ensembl_only.csv.
# Usage: Rscript fix_growth_inputs.R <v2 dir> <dir with DDG2P.csv.gz>
suppressMessages(library(data.table))
a <- commandArgs(TRUE); V <- a[1]; RF <- a[2]; I <- file.path(V, "inputs")
M_ <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly"; C_ <- file.path(M_, "CNV_ASD_growth_analysis/concordance")
G_ <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByMolly/cnv-growth-bmi/concordance"
genes <- unique(fread(file.path(M_, "pathway_analysis/Final/cnv_genes_annotated.csv"))[, .(CNV_ID, Gene)])
hg <- fread(file.path(G_, "hgnc_map.tsv"), colClasses = "character"); sh <- fread(file.path(G_, "GeneBayes_shet_Zeng2024.tsv"))
sh[, gpct := frank(post_mean, na.last = "keep") / sum(!is.na(post_mean))]
al <- rbind(hg[, .(sym = symbol, hgnc_id, ensg = ensembl_gene_id)],
            hg[!is.na(prev_symbol) & prev_symbol != "", .(sym = unlist(strsplit(prev_symbol, "[|]"))), by = .(hgnc_id, ensg = ensembl_gene_id)],
            hg[!is.na(alias_symbol) & alias_symbol != "", .(sym = unlist(strsplit(alias_symbol, "[|]"))), by = .(hgnc_id, ensg = ensembl_gene_id)])[!duplicated(sym)]
look <- function(sym, fix) { e <- al$ensg[match(sym, al$sym)]; i <- match(e, sh$ensg)
  if (fix) { j <- is.na(i); i[j] <- match(al$hgnc_id[match(sym[j], al$sym)], sh$hgnc) }; i }
tw <- function(t) fread(file.path(C_, "twas30", paste0(t, ".txt")))[, .(z = mean(TWAS.Z, na.rm = TRUE)), by = GENE]
sums <- function(fix) rbindlist(lapply(c(Height = "height", "BMI (TWAS30)" = "bmi"), function(t) {
  z <- tw(t); g <- copy(genes); g[, z := z$z[match(Gene, z$GENE)]]; g[, shet := sh$post_mean[look(Gene, fix)]]
  g[!is.na(z), .(M = sum(z[!is.na(shet) & shet < .029]), M_abs = sum(abs(z[!is.na(shet) & shet < .029])), K_twas = sum(abs(z[!is.na(shet) & shet >= .029])),
                 K_n = sum(!is.na(shet) & shet >= .029), S_con = sum(z[!is.na(shet) & shet >= .029]), S_shet = sum(z[!is.na(shet)])), by = .(NAME = CNV_ID)] }), idcol = "disorder")
old <- sums(FALSE); new <- sums(TRUE)
S <- fread(file.path(I, "growth_input_signed.csv")); T <- fread(file.path(I, "growth_trait_input.csv"))
chk <- old[S, on = .(NAME, disorder), nomatch = 0]
cat("old-mapping recompute vs current inputs: max |dM| =", signif(max(abs(chk$M - chk$i.M), na.rm = TRUE), 3), " max |dK| =", signif(max(abs(chk$K_twas - chk$i.K_twas), na.rm = TRUE), 3), "\n")
# dosage-sensitive gene counts (all genes, fixed mapping) and DDG2P dosage genes
dd <- fread(text = readLines(gzfile(file.path(RF, "DDG2P.csv.gz"))), sep = ","); setnames(dd, tolower(gsub("[^A-Za-z0-9]+", "_", names(dd))))
ddm <- dd[grepl("definitive|strong", tolower(confidence)) & grepl("monoallelic", tolower(allelic_requirement))]
dd_lof <- unique(ddm$gene_symbol[grepl("loss of function", ddm$molecular_mechanism, ignore.case = TRUE) | grepl("absent gene product", ddm$variant_consequence, ignore.case = TRUE)])
g <- copy(genes); g[, gp := sh$gpct[look(Gene, TRUE)]]
D <- g[, .(D67 = sum(gp > .667, na.rm = TRUE), D80 = sum(gp > .8, na.rm = TRUE), D90 = sum(gp > .9, na.rm = TRUE), D95 = sum(gp > .95, na.rm = TRUE),
           DDG2P_dosage = sum(Gene %in% dd_lof)), by = .(NAME = CNV_ID)]
fwrite(D, file.path(I, "growth_D_counts.csv"))
# write corrected inputs (back up the old ones once)
for (f in c("growth_input_signed.csv", "growth_trait_input.csv")) { b <- file.path(I, sub(".csv", "_v1_ensembl_only.csv", f)); if (!file.exists(b)) file.copy(file.path(I, f), b) }
cols <- c("M", "M_abs", "K_twas", "K_n", "S_con", "S_shet")
for (cc in cols) { S[new, on = .(NAME, disorder), (cc) := get(paste0("i.", cc))] }
S[disorder %in% c("Height", "BMI (TWAS30)"), `:=`(M_check = M, K_check = K_twas)]
S <- D[S, on = "NAME"]; fwrite(S, file.path(I, "growth_input_signed.csv"))
T[, disorder := fifelse(trait == "Height", "Height", "BMI (TWAS30)")]
for (cc in intersect(cols, names(T))) T[new, on = .(NAME, disorder), (cc) := get(paste0("i.", cc))]
T[, disorder := NULL]; T <- D[T, on = "NAME"]; fwrite(T, file.path(I, "growth_trait_input.csv"))
st <- fread(file.path(I, "growth_strata_input.csv")); st <- D[st[, setdiff(names(st), names(D)[-1]), with = FALSE], on = "NAME"]; fwrite(st, file.path(I, "growth_strata_input.csv"))
d <- merge(old, new, by = c("disorder", "NAME"), suffixes = c(".old", ".new"))
print(d[abs(M.old - M.new) + abs(K_twas.old - K_twas.new) > 1e-9, .(disorder, NAME, M.old = round(M.old, 2), M.new = round(M.new, 2), K.old = round(K_twas.old, 2), K.new = round(K_twas.new, 2))])
cat("D counts:\n"); print(D[order(-D95)])
