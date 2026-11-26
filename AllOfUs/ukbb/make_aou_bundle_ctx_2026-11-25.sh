#!/bin/bash
# Add-on resources bundle for the All of Us run of the person-level ConGA context test (run_pipeline.py step 'ctx'): the context and
# signal-group gene sets, their BH families, the locus-removal units and set list, the loci file and the two UK Biobank shift files the
# UKBB ctx runs used (code/resources/shared_shifts/{height,bmi}_shifts.npy; NOT in the 2026-11-16 bundle). Gene lists and shift offsets
# only; no individual data. Unpacks into AllOfUs/resources/ next to the 2026-11-16 bundle (which supplies conga/hg38_refGene exons).
# null model: n/a (packaging; the shift files define the genome-shift null shared by UKBB and All of Us)
#   bash make_aou_bundle_ctx_2026-11-25.sh        (Expanse login node; ~10 MB)
set -euo pipefail
P=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv; R=$P/results/cnv_gwas_meta_2026-11-24
B=/expanse/projects/sebat1/jsebat/sacks2026_aou/aou_bundle; N=cnv_resources_ctx_2026-11-25
if [ -e $B/$N.tar.gz ]; then echo "$B/$N.tar.gz exists; remove it first if you mean to rebuild"; exit 1; fi
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1      # login-node memory cap
T=$B/$N.tmp; rm -rf $T; mkdir -p $T/resources/conga_ctx $T/resources/shared_shifts
cp $R/conga_supp_56063298/sets/all_sets_membership.tsv.gz $R/conga_supp_56063298/sets/all_sets_families.tsv $T/resources/conga_ctx/
cp $P/code/meta_aou_2026-11-24/loo_units_meta12.tsv $P/code/resources/growth_loci/loci_sacks34.tsv $T/resources/conga_ctx/
cat $R/conga_supp_56063298/contrasts/loo_set_list.txt $R/conga_fig16_56072429/loo_set_list.txt | sort -u > $T/resources/conga_ctx/loo_set_list.txt
cp $P/code/resources/shared_shifts/height_shifts.npy $P/code/resources/shared_shifts/bmi_shifts.npy $T/resources/shared_shifts/
# the shift files must be the ones the UKBB runs used (info.json shifts_sha1 = first 16 hex of sha1 of the int64 array)
PY=/home/jsebat/micromamba/envs/dnmt3a-py/bin/python
$PY - $T/resources/shared_shifts $R <<'PYEOF'
import hashlib, json, sys, numpy as np
sh, R = sys.argv[1:]
for t in ["height", "bmi"]:
    h = hashlib.sha1(np.load(f"{sh}/{t}_shifts.npy").astype(np.int64).tobytes()).hexdigest()[:16]
    ref = {json.load(open(f"{R}/{run}/full_{t}_{d}/info.json"))["shifts_sha1"] for run in ("conga_ctx_56044342", "conga_supp_56063298") for d in ("DEL", "DUP")}
    print(t, h, "UKBB runs:", ref)
    if ref != {h}:
        sys.exit(f"{t}: shift file differs from the UKBB ctx runs")
open(f"{sh}/SHIFTS_MANIFEST_ctx.tsv", "w").write("file\tshifts_sha1_16\tmatch_ukbb_ctx_runs\n" + "".join(
    f"{t}_shifts.npy\t{hashlib.sha1(np.load(f'{sh}/{t}_shifts.npy').astype(np.int64).tobytes()).hexdigest()[:16]}\tyes\n" for t in ["height", "bmi"]))
PYEOF
cd $T && tar czf $B/$N.tar.gz resources && cd $B && rm -rf $T
sha256sum $N.tar.gz > $N.sha256 && ls -la $N.* && tar tzf $N.tar.gz
chmod -R g+rX $B
