#!/bin/bash
# Run on an Expanse compute node (small): UKBB reference probe positions + shared genome shifts per trait (2,000) for the AoU meta-analysis.
set -e
P=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv; R=$P/code/resources; PY=/home/jsebat/micromamba/envs/dnmt3a-py/bin/python
mkdir -p $R/shared_shifts
$PY - <<'PY'
import pandas as pd
P='/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/results/pipeline_ukbb/taskA_gwas'
L=[]
for t in ['height','bmi']:
    for x in ['del','dup']:
        d=pd.read_csv(f'{P}/{t}/sumstats_{x}.tsv.gz',sep='\t',usecols=['CHR','BP']); d=d[d.CHR<=22]; d['TYPE']=x.upper(); L.append(d)
R=pd.concat(L).drop_duplicates().sort_values(['TYPE','CHR','BP'])
R.to_csv('/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/code/resources/ukbb_ref_probes.tsv.gz',sep='\t',index=False)
print('ref probes',R.groupby('TYPE').size().to_dict())
PY
cd $P/code/cnv_anthro
for t in height bmi; do
  $PY conga_grid_qt.py --pheno $P/data/ukbb_pheno_std.tsv --trait ${t}_rint --cnv $P/data/ukbb_pgc_xd.cnv --twas $R/twas/$t.txt --panels gtex_nonbrain \
    --exons $R/conga/hg38_refGene_20200708.exon.txt --gene_cov $R/conga/gene_covariates.tsv --clingen $R/conga/ClinGen_gene_curation_list_GRCh38.tsv \
    --loci $R/growth_loci/loci_sacks34.tsv --loo_runs $R/growth_loci/loo_runs_sacks34.tsv --design designs/grid_univariate.json --only "none_H|adj" \
    --nrep 2000 --seed 2026 --perms_file $R/shared_shifts/${t}_shifts.npy --perms_only --out $P/results/shared_shifts_info/$t
done
chmod -R g+rX $R
