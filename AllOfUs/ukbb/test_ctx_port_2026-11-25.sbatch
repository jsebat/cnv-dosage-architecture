#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH -t 02:00:00
#SBATCH -J ctx_port_test
#SBATCH -o /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/ctx_port_test_2026-11-25.%j.log
# Tests of the All of Us port of the person-level ConGA context test:
#  (1) UK Biobank reproduction: the ported conga_ctx_qt.py on the real UKBB data for a handful of sets (smallest P signal groups and
#      contexts) must give the same z, z_emp and P_emp as the UKBB runs behind Supplementary Fig. 1 (full run: height DEL, BMI DUP;
#      locus-removal runs: height DEL);
#  (2) end to end on synthetic cohorts: run_pipeline.py --steps ctx,export on two synthetic cohorts, export checks, meta_ctx.py (A+A and
#      A+B) checks, and conga_contrasts.py stage B drawing the figure from the meta output.
# Args: $1 = repo snapshot (holding AllOfUs/, recurrentCNVgenotyping/, cnv_gwas/meta_aou_2026-11-24/conga_contrasts.py). New output dir.
# null model: genome shift of the gene annotation with the shared UKBB shifts (resources/shared_shifts); matched in removal runs
set -euo pipefail
REPO=$1
P=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv; RS=$P/code/resources; R=$P/results/cnv_gwas_meta_2026-11-24
O=$P/results/ctx_port_test_${SLURM_JOB_ID}; mkdir -p $O
CA=$REPO/AllOfUs/cnv_anthro; PY=/home/jsebat/micromamba/envs/dnmt3a-py/bin/python
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 MPLBACKEND=Agg
SETS=$R/conga_supp_56063298/sets/all_sets_membership.tsv.gz; FAM=$R/conga_supp_56063298/sets/all_sets_families.tsv
LOOU=$P/code/meta_aou_2026-11-24/loo_units_meta12.tsv
COMMON="--exons $RS/conga/hg38_refGene_20200708.exon.txt --loci34 $RS/growth_loci/loci_sacks34.tsv --shifts $RS/shared_shifts/{trait}_shifts.npy"
UK="--pheno $P/data/ukbb_pheno_std.tsv --cnv $P/data/ukbb_pgc_xd.cnv --families $FAM"
cat $R/conga_supp_56063298/contrasts/loo_set_list.txt $R/conga_fig16_56072429/loo_set_list.txt | sort -u > $O/ukbb_loo_sets.txt

echo "=== (1) UK Biobank reproduction"
for TD in height_DEL bmi_DUP; do
  t=${TD%_*}; d=${TD#*_}
  $PY $CA/tests/test_ctx_port.py pick --ukbb_full $R/conga_supp_56063298/full_$TD,$R/conga_ctx_56044342/full_$TD --sets $SETS --n 3 --out $O/pick_$TD
  $PY $CA/conga_ctx_qt.py $UK $COMMON --trait $t --dosage $d --sets $O/pick_$TD/membership.tsv.gz --out $O/port_full_$TD
  $PY $CA/tests/test_ctx_port.py repro --mode full --port $O/port_full_$TD --ukbb $R/conga_supp_56063298/full_$TD,$R/conga_ctx_56044342/full_$TD
done
grep -Fxf $O/ukbb_loo_sets.txt $O/pick_height_DEL/sets.txt > $O/loo_sets_height_DEL.txt
$PY $CA/conga_ctx_qt.py $UK $COMMON --trait height --dosage DEL --sets $O/pick_height_DEL/membership.tsv.gz --out $O/port_loo_height_DEL \
   --runs loo --set_list $O/loo_sets_height_DEL.txt --loo_units $LOOU
$PY $CA/tests/test_ctx_port.py repro --mode loo --port $O/port_loo_height_DEL \
   --ukbb $R/conga_supp_56063298/loo_height_DEL,$R/conga_fig16_56072429/loo_height_DEL

echo "=== (2) synthetic cohorts end to end"
S=$O/synthetic; mkdir -p $S
$PY - "$SETS" "$S" <<'PYEOF'
import sys, pandas as pd
sm = pd.read_csv(sys.argv[1], sep="\t", dtype=str); out = sys.argv[2]
ids = pd.Series(sm.set_id.unique())
keep = list(ids[ids.str.startswith("SG:") & ids.str.endswith(("|HC", "|LC"))])
for pre in ["gtex:", "hpa_sc:", "limb:", "brain:BSpan", "brain:Wang"]:
    u = sorted({s.rsplit("|", 1)[0] for s in ids if s.startswith(pre)})[:4]
    keep += [f"{x}|{c}" for x in u for c in ("HC", "LC")]
sm[sm.set_id.isin(keep)].to_csv(f"{out}/membership.tsv.gz", sep="\t", index=False)
pd.Series(["SG:height:2|HC", "SG:height:2|LC", "SG:bmi:4|HC", "SG:bmi:4|LC"]).to_csv(f"{out}/loo_sets.txt", index=False, header=False)
print("synthetic membership sets:", len(set(keep) & set(ids)))
PYEOF
for C in A:1 B:2; do
  c=${C%:*}; seed=${C#*:}
  $PY $CA/tests/make_synthetic_ctx.py --exons $RS/conga/hg38_refGene_20200708.exon.txt --loci34 $RS/growth_loci/loci_sacks34.tsv \
     --sets $S/membership.tsv.gz --signal_set "SG:height:2|HC" --n 4000 --seed $seed --out $S/data_$c
  $PY - "$S" "$c" "$RS" "$FAM" "$LOOU" <<'PYEOF'
import json, sys
S, c, RS, FAM, LOOU = sys.argv[1:]
pcs = [f"PC{i}" for i in range(1, 17)]
C = {"cohort": f"SYNTHETIC_{c}", "results_dir": f"{S}/results_{c}", "export_dir": f"{S}/export_{c}", "min_cell": 20,
     "pheno_table": f"{S}/data_{c}/pheno.tsv", "cnv_table": f"{S}/data_{c}/calls.cnv", "prep_commands": [],
     "traits": {"height": {"column": "height_rint", "numeric": ["age", "age2"] + pcs, "factors": ["sex"]},
                "bmi": {"column": "bmi_rint", "numeric": ["age_bmi", "age_bmi2"] + pcs, "factors": ["sex"]}},
     "ctx": {"strata": ["all"], "traits": ["height", "bmi"], "dosages": ["DEL", "DUP"], "sets": f"{S}/membership.tsv.gz",
             "families": FAM, "exons": f"{RS}/conga/hg38_refGene_20200708.exon.txt", "loci34": f"{RS}/growth_loci/loci_sacks34.tsv",
             "shifts": f"{RS}/shared_shifts/{{trait}}_shifts.npy", "loo_units": LOOU, "loo_set_list": f"{S}/loo_sets.txt",
             "nrep": 200, "jobs": 4}}
json.dump(C, open(f"{S}/config_{c}.json", "w"), indent=1)
PYEOF
  $PY $CA/run_pipeline.py --config $S/config_$c.json --steps ctx,export --threads 8
  $PY $CA/tests/test_ctx_port.py export --results $S/results_$c/ctx --export $S/export_$c/ctx
done
$PY $CA/meta_ctx.py --cohort A=$S/export_A/ctx --cohort A2=$S/export_A/ctx --out $S/meta_AA
$PY $CA/meta_ctx.py --cohort A=$S/export_A/ctx --cohort B=$S/export_B/ctx --out $S/meta_AB
$PY $CA/tests/test_ctx_port.py meta --a $S/export_A/ctx --aa $S/meta_AA --ab $S/meta_AB
$PY $REPO/cnv_gwas/meta_aou_2026-11-24/conga_contrasts.py --stage B --ctx_run $S/meta_AB --loo_run $S/meta_AB \
   --groups $R/conga_supp_56063298/groups --out $S/contrasts_AB --cohort_label "synthetic A + B"
ls $S/contrasts_AB
echo "signal set SG:height:2|HC in the synthetic meta (planted dose-aligned effect on height):"
zcat $S/meta_AB/full_height_DEL/conga_ctx_full.tsv.gz | grep -F "SG:height:2|" || true
zcat $S/meta_AB/full_height_DUP/conga_ctx_full.tsv.gz | grep -F "SG:height:2|" || true
chmod -R g+rX $O; echo ALL_DONE
