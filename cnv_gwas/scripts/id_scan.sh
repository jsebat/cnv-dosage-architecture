#!/bin/bash
# Pre-push scan for the ukbb_aou subtree: fails on individual-level files, ID-like columns, study IDs or credentials.
# Usage: bash cnv_gwas/scripts/id_scan.sh   (from the repository root)
set -u
D=${1:-cnv_gwas}
PAT='(MONT|UCLA|PENN|TORO)[0-9]{3,4}|NDAR[A-Z0-9]{8}|SP0[0-9]{6}|SF0[0-9]{5}|REACH[0-9]{4,}|SDSM-[A-Z0-9]{4}|AKIA[0-9A-Z]{16}|aws_secret_access_key|BEGIN (RSA|OPENSSH) PRIVATE KEY'
fail=0
hits=$(grep -rnE "$PAT" "$D" --exclude=id_scan.sh 2>/dev/null); [ -n "$hits" ] && { echo "ID/credential patterns:"; echo "$hits"; fail=1; }
indiv=$(find "$D" -type f \( -name '*.cnv' -o -name '*.phe' -o -name '*.indiv' -o -name '*.vcf*' -o -name '*.fam' -o -name '*.bed' -o -name '*.bim' -o -name '*pheno_std*' \)); [ -n "$indiv" ] && { echo "Individual-level file types:"; echo "$indiv"; fail=1; }
for f in $(find "$D/aggregates" -type f \( -name '*.tsv' -o -name '*.csv' -o -name '*.tsv.gz' \) 2>/dev/null); do
  h=$( (zcat -f "$f" | head -1) 2>/dev/null)
  echo "$h" | tr '\t,' '\n\n' | grep -qixE 'IID|FID|eid|person_id|research_id|sample|sample_id' && { echo "ID column in $f"; fail=1; }
done
# AoU aggregates: no count from 1 to 20 (export_aggregate.py writes '<=20')
for f in $(find "$D/aggregates/aou" -type f 2>/dev/null); do echo "check AoU small cells by hand: $f"; done
[ $fail = 0 ] && echo "id_scan: clean" || exit 1
