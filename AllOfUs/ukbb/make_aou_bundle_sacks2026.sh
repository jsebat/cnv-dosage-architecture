#!/bin/bash
# Resources bundle for the All of Us run (public annotations, locus files, UKBB shift lists, UKBB probe positions; no individual data).
# Run on Expanse after make_ukbb_shifts.sbatch finished with every SHIFTS_MANIFEST row 'yes' (light work: tar of ~15 MB, fine on a login node).
#   OUT_ROOT=<same as for make_ukbb_shifts.sbatch> bash make_aou_bundle_sacks2026.sh
set -e
P=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv; OUT_ROOT=${OUT_ROOT:-/expanse/projects/sebat1/$USER/sacks2026_aou}
B=$OUT_ROOT/aou_bundle; N=cnv_resources_sacks2026_2026-11-16
SH=$OUT_ROOT/shared_shifts_sacks2026_2026-11-16/shared_shifts
if grep -q -P '\tNO$' $SH/SHIFTS_MANIFEST.tsv; then echo "SHIFTS_MANIFEST has mismatches; stop"; exit 1; fi
if [ -e $B/$N.tar.gz ]; then echo "$B/$N.tar.gz exists; remove it first if you mean to rebuild"; exit 1; fi
mkdir -p $B/$N.tmp/resources/shared_shifts
cd $P/code/resources
cp -r conga omar_alleles twas growth_loci segDups.bed ukbb_ref_probes.tsv.gz Collins_rCNV_2022.dosage_sensitivity_scores.tsv $B/$N.tmp/resources/
cp $SH/*.npy $SH/SHIFTS_MANIFEST.tsv $B/$N.tmp/resources/shared_shifts/
cd $B/$N.tmp && tar czf $B/$N.tar.gz resources && cd $B && rm -rf $B/$N.tmp
sha256sum $N.tar.gz > $N.sha256 && ls -la $N.* && tar tzf $N.tar.gz | sed 's#/[^/]*$##' | sort | uniq -c
chmod -R g+rX $B
