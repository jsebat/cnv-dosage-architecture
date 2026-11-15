#!/bin/bash
# Resources bundle for the AoU run (aggregate / public files only): annotations, locus files, shared UKBB shifts, UKBB probe positions.
set -e
P=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv; B=$P/aou_bundle; rm -rf $B; mkdir -p $B
cd $P/code
tar czf $B/cnv_gwas_resources_2026-11-15.tar.gz resources/conga resources/omar_alleles resources/twas resources/growth_loci resources/shared_shifts resources/segDups.bed resources/ukbb_ref_probes.tsv.gz
cd $B && sha256sum cnv_gwas_resources_2026-11-15.tar.gz > cnv_gwas_resources_2026-11-15.sha256 && ls -la && tar tzf cnv_gwas_resources_2026-11-15.tar.gz | sed 's#/[^/]*$##' | sort | uniq -c
chmod -R g+rX $B
