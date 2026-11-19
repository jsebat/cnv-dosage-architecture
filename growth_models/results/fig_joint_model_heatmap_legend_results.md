# Joint-model summary figure: legend and results text (draft, 2026-11-15)

Figure: `fig_joint_model_heatmap.pdf` (code `v2_2026-11-08/code/fig_joint_model_heatmap.R`; cell values `joint_model_heatmap.csv`).

## Figure legend

**The gene dosage architecture of height and BMI, and its robustness.** Each column is one joint model, in which every term below is fitted together. Left, the 34 recurrent CNV intervals with published height and BMI effects (Sacks et al.); right, all CNVs called genome-wide in the UK Biobank. Within each block, the first column uses all CNVs. The other columns refit the same model after removing (i) duplications that contain a single dosage-sensitive gene, (ii) all CNVs that contain one, (iii) the 13 asymmetric loci that lie off the mirror line, or (iv) both (i) and (iii). Colour is the z score of each term. Grey means the term is not part of that model. A dash means the term cannot be estimated by design: the CNVs that define it were removed, or (UK Biobank, all CNVs) it was fitted as an adjustment and not estimated.

Rows are grouped by the axis of the dose response that each term captures (labels on the left): **Mirror**, **Constraint**, **Single hit** and **Mean shift**. Within each axis, "Deletion" and "Duplication" are the terms for each type of CNV.

Terms (rows), by axis:
- *Mirror:* deletion and duplication of the same genes move the trait in opposite directions.
- *Mirror, total:* the negative correlation between the deletion and the duplication effects of the same loci.
- *Mirror explained by TWAS:* the part predicted by genes whose expression is associated with the trait in common-variant studies (transcriptome-wide association, TWAS). It is the difference between the duplication and the deletion terms of the TWAS gene count, (z_DUP − z_DEL)/√2. Its sign shows the direction: positive when duplication raises and deletion lowers the trait.
- *Mirror not explained by TWAS:* the negative correlation that remains once the TWAS term is in the model.
- *Explained by TWAS genes, deletions / duplications:* the number of TWAS genes in the CNV (|TWAS z| > 3, regardless of sign; UK Biobank: constrained genes only. Height uses the summed |z| over GTEx non-brain panels. BMI uses the count of genes with |mean z| > 4 over all 45 panels, per Jon 2026-11-16, after the stratified cross-validation found no tissue set better than all panels.)
- *Single hit:* the CNV contains at least one gene in the top 5% of selective constraint (s_het) or haploinsufficient in ClinGen. The same union criterion is used in every column (UK Biobank BMI removal columns rerun 2026-11-16).
- *Constraint:* the number of constrained genes in the CNV (s_het > 0.029).
- *Mean shift:* the effect of the total number of genes in the CNV, a shift in the mean effect that grows with CNV size. The size-independent shift (the intercepts) was not significant for either trait at locus level (|z| ≤ 1.04).

Statistics: z against a null in which all gene annotations are shifted together along the genome, so that every CNV keeps its size (1,000–2,000 replicates; mirror total and not explained by TWAS: likelihood-ratio test). * P < 0.05; ** false discovery rate < 0.05 across all tiles. Black box: P < 0.05 and the term stays significant, with the same sign, when each recurrent locus is removed in turn (34 intervals; 24 recurrent units in the UK Biobank), each removal tested against a null with the same locus removed.

## Results text

**The mirror and the single dosage-sensitive gene describe the recurrent CNVs.** In the 34 recurrent intervals, the deletion and duplication effects of the same loci were strongly anticorrelated for both height (z = 6.1) and BMI (z = 5.9), and the anticorrelation held when each locus was removed in turn. For height, genes associated with the trait in common-variant studies explained about half of this mirror (58%; z = 3.7 for the TWAS term). Most of it came through deletions (z = −3.6), and it remained significant after removing duplications of a single dosage-sensitive gene (z = 2.4), all CNVs carrying one (z = 2.7), or the 13 asymmetric loci (z = 3.2). The remaining height effect was a lower duplication effect at loci that carry a single dosage-sensitive gene (z = −3.4), which places those loci on a second line below the mirror line. For BMI, TWAS genes explained little of the mirror (17%), and most of it (z = 5.5) was not explained by TWAS. Neither the number of constrained genes nor the total number of genes had an effect at locus level for either trait.

**Constraint and CNV size act genome-wide.** Across all CNVs in the UK Biobank, the number of constrained genes lowered height for both deletions (z = −2.6) and duplications (z = −3.6) and raised BMI for duplications (z = 3.1). These effects held when every recurrent locus was removed in turn and after every removal of CNV classes, so they do not depend on the recurrent CNVs or on single dosage-sensitive genes. The genome-wide BMI mirror explained by TWAS was strong (z = −8.1; TWAS genes from all 45 expression panels). It survived removal of duplications that carry a single dosage-sensitive gene (z = −6.4), through the deletion arm (z = 7.8). It did not survive removal of all CNVs that carry one (z = +0.7). Genome-wide, the height mirror was weak (z = 2.2) and depended on one locus, 16p11.2 BP4–BP5.

## Open before finalising
- Done 2026-11-16: BMI genome-wide TWAS measure = all 45 panels (stratified CV: no tissue set beat all panels out of sample, BMI gain z +0.34, P .30); BMI removal columns rerun with the union single-hit term (bmi_allpanels_union_2026-11-16). All of Us replication: code in GitHub AllOfUs/, handoff to Milad.
- The tissue grid is exploratory and goes to the supplement, with the caveat that no tissue was validated out of sample.
- Significance of the total and residual mirror within the robustness strata (locus level) is not yet computed (grey cells).
- All of Us replication: AllOfUs/ code + handoff (2026-11-16); meta-analysis panels pending Milad's export.
