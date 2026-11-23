#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH -t 00:30:00
#SBATCH -J cnv_filter_ukbb
#SBATCH -o /expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23/logs/%x.%j.log
# Rare-CNV filter for Fig. 3 (Jon 2026-11-23), UK Biobank: per-position carrier frequency of the PGC-freeze call set and the calls / people the
# filter removes at 0.1%, 0.5%, 1%, 2% and 5% (cnv_freq_filter.py; a call is dropped when >= half of its length has same-dosage carrier
# frequency above the threshold). Tabulation only: no model is refitted. MAXF = threshold of the filtered table written (default 0.02).
# null model: none (deterministic filter). Unique output dir: $B/ukbb_cnv_filter_$SLURM_JOB_ID.
set -e
B=/expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23; P=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv
PY=${PY:-/home/jsebat/micromamba/envs/dnmt3a-py/bin/python}; REPO=${REPO:-$B/code/AllOfUs}; W=$B/ukbb_cnv_filter_$SLURM_JOB_ID; mkdir -p $W
cp $REPO/cnv_anthro/cnv_freq_filter.py $REPO/cnv_anthro/common.py $W/; cd $W; export OPENBLAS_NUM_THREADS=2
$PY cnv_freq_filter.py --cnv $P/data/ukbb_pgc_xd.cnv --pheno $P/data/ukbb_pheno_std.tsv --max_freq ${MAXF:-0.02} --report 0.001,0.005,0.01,0.02,0.05 \
   --out_cnv $W/ukbb_pgc_xd_rare.cnv --out_dir $W
chmod -R g+rX $W; chmod o-rwx $W/ukbb_pgc_xd_rare.cnv
# which of the 34 recurrent intervals lose calls (aggregate counts per locus; UK Biobank, stays at UCSD)
cp $REPO/cnv_anthro/conga_axes_qt.py $W/; $PY - <<PYEOF
import pandas as pd, numpy as np
from common import read_plink_cnv
from conga_axes_qt import assign_loci
R = "$P/code/resources/growth_loci/loci_sacks34.tsv"; ids = pd.read_csv("$P/data/ukbb_pheno_std.tsv", sep="\t", usecols=["IID"], dtype=str).IID
a = read_plink_cnv("$P/data/ukbb_pgc_xd.cnv", keep_iids=ids); a = a[a.CHR <= 22].reset_index(drop=True); b = pd.read_csv("$W/ukbb_pgc_xd_rare.cnv", sep="\t", dtype={"IID": str})
loci = pd.read_csv(R, sep="\t"); a["locus"] = assign_loci(a, loci); key = ["IID", "CHR", "BP1", "BP2", "TYPE"]
a = a.merge(b[key].drop_duplicates().assign(kept=True), on=key, how="left"); a["kept"] = a.kept.fillna(False).astype(bool)
t = a[a.locus != ""].groupby("locus").agg(calls=("kept", "size"), kept=("kept", "sum")); t["dropped"] = t.calls - t.kept; t["freq_pct"] = (100 * t.calls / len(ids)).round(3)
print("calls assigned to the 34 recurrent intervals:", int(t.calls.sum()), "| dropped by the filter:", int(t.dropped.sum()))
print(t[t.dropped > 0].sort_values("dropped", ascending=False).to_string()); t.to_csv("$W/recurrent_loci_calls_dropped.tsv", sep="\t")
PYEOF
