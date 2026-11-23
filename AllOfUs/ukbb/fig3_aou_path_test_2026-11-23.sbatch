#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH -t 05:00:00
#SBATCH --nice=100
#SBATCH -J fig3_path_test
#SBATCH -o /expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23/logs/%x.%j.log
# Sacks 2026 dosage paper, Fig. 3 joint model (Jon 2026-11-23): end-to-end test of the code Milad runs in All of Us, on UK Biobank data.
# run_pipeline.py --steps fig3 (= run_fig3.py: 8 joint-model fits + 2 gene-mirror runs, from the committed spec) -> export -> UKBB-only meta ->
# comparison with the meta of the finished UK Biobank runs (REF, from fig3_ukbb_export_meta_2026-11-23.sbatch). Also writes a second export
# with the All of Us small-cell rule (min_cell 20) to exercise the suppression and its audit.
# null model: genome shift of the gene annotation, 2,000 shared UK Biobank shifts (sha1 checked by run_fig3.py); gene mirror: gene-circle shift.
# Unique output dir: $B/aou_path_test_$SLURM_JOB_ID (code snapshot inside).
#   sbatch --export=ALL,REF=<meta_ukbb_only_<jobid>> fig3_aou_path_test_2026-11-23.sbatch
set -e
B=/expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23; REF=${REF:?set REF=<UKBB-only meta dir of the finished runs>}
PY=${PY:-/home/jsebat/micromamba/envs/dnmt3a-py/bin/python}; REPO=${REPO:-$B/code/AllOfUs}; export PATH=$(dirname $PY):$PATH
W=$B/aou_path_test_$SLURM_JOB_ID; mkdir -p $W; S=$W/code_snapshot; mkdir -p $S; cp -r $REPO/cnv_anthro $REPO/meta $REPO/ukbb $S/; cd $S/cnv_anthro
$PY - <<EOF
import json
C = json.load(open("configs/ukbb_expanse.json"))
C["results_dir"] = "$W/results"; C["export_dir"] = "$W/export"; C["min_cell"] = 0; C["strata"] = {"all": {}}
C["cohort"] = "UKBB_via_fig3_code_path_TEST"
json.dump(C, open("$W/test_config.json", "w"), indent=1)
C["min_cell"] = 20; C["fig3"]["export_dir"] = "$W/export_mincell20"; json.dump(C, open("$W/test_config_mincell20.json", "w"), indent=1)
EOF
$PY run_pipeline.py --config $W/test_config.json --steps fig3,fig3_export --strata all --jobs 5
$PY run_fig3.py --config $W/test_config_mincell20.json --steps export
$PY ../meta/meta_fig3.py --cohort UKBB=$W/export_fig3/all/fig3 --out $W/meta_ukbb_only --code $S/cnv_anthro
set +e
$PY ../meta/compare_fig3_tables.py --a $REF --b $W/meta_ukbb_only > $W/compare_with_finished_runs.txt 2>&1
echo "compare exit $?"; cat $W/compare_with_finished_runs.txt
chmod -R g+rX $W
