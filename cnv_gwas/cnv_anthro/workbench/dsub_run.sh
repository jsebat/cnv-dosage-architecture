#!/bin/bash
# dsub launcher for the AoU Workbench (run from a Workbench terminal). One task per step / trait x panel set.
# Code + resources are staged in ${WORKSPACE_BUCKET}/cnv_anthro_code; run outputs go to ${WORKSPACE_BUCKET}/cnv_anthro_run.
# Uses the Workbench's dsub helper conventions (provider google-cls-v2, user project = GOOGLE_PROJECT).
set -euo pipefail
STEP=${1:?usage: dsub_run.sh <steps> [trait] [panels]}; TRAIT=${2:-}; PANELS=${3:-}
IMAGE=python:3.11-slim
dsub --provider google-cls-v2 --project "${GOOGLE_PROJECT}" --user-project "${GOOGLE_PROJECT}" \
  --regions us-central1 --machine-type n2-standard-8 --disk-size 200 --image "$IMAGE" \
  --network network --subnetwork subnetwork --service-account "$(gcloud config get-value account)" \
  --logging "${WORKSPACE_BUCKET}/dsub/logs/" \
  --env STEP="$STEP" --env TRAIT="$TRAIT" --env PANELS="$PANELS" --env WORKSPACE_CDR="${WORKSPACE_CDR}" --env GOOGLE_PROJECT="${GOOGLE_PROJECT}" \
  --input-recursive CODE="${WORKSPACE_BUCKET}/cnv_anthro_code" \
  --input-recursive RUNIN="${WORKSPACE_BUCKET}/cnv_anthro_run" \
  --output-recursive RUNOUT="${WORKSPACE_BUCKET}/cnv_anthro_run" \
  --command 'pip install -q numpy scipy pandas matplotlib google-cloud-bigquery db-dtypes && \
             mkdir -p ${HOME}/cnv_anthro_run && cp -r ${RUNIN}/. ${HOME}/cnv_anthro_run/ 2>/dev/null || true; \
             cd ${CODE}/cnv_anthro && python run_pipeline.py --config configs/aou_workbench.json --steps ${STEP} \
               ${TRAIT:+--only_trait ${TRAIT}} ${PANELS:+--only_panels ${PANELS}} && cp -r ${HOME}/cnv_anthro_run/. ${RUNOUT}/'
