"""Regenerate the UK Biobank genome-shift lists used by every UKBB run that enters the UKBB + All of Us meta-analysis, and verify each one.

The shifts of conga_grid_qt.py are drawn per replicate from the seed (seed * 100003 + rep) and kept when the burden of the shifted
annotation is within 10% of the observed burden. The burden comes from the CNVs of the people analysed, so the shift list depends on
the trait (who has a phenotype), the TWAS panel set (which genes count as TWAS genes) and the carrier subset (who was removed).
For each combination used in a UKBB result, this script reruns conga_grid_qt.py --perms_only with the original inputs, saves
<trait>_<panels>_<subset>.npy (2,000 integer offsets; public, no individual data) and checks its sha1 against the perms_sha1 that the
original UKBB run recorded in run_info.json. All of Us loads these files, so replicate i is the same annotation shift in both cohorts.

Run on Expanse inside a job (make_ukbb_shifts.sbatch). Writes <out>/<trait>_<panels>_<subset>.npy and <out>/SHIFTS_MANIFEST.tsv.
Person-level strata files are written to <work>/strata and deleted at the end.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import numpy as np
import pandas as pd

P = "/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv"
R = f"{P}/code/resources"
RES = f"{P}/results"
# (trait, panels, subset) -> UKBB run_info.json files whose perms_sha1 the regenerated list must reproduce (any nrep; checked on offs[:nrep])
SPECS = [
    ("height", "gtex_nonbrain", "all", ["union_terms_2026-11-15/height", "growthloci_2026-11-15/grid/height/none_H__adj",
                                         "growthloci_2026-11-15/joint_sacks34/height/height_selected"]),
    ("bmi", "gtex_nonbrain", "all", ["union_terms_2026-11-15/bmi", "union_terms_bmi_fix_2026-11-15/bmi", "growthloci_2026-11-15/grid/bmi/none_H__adj",
                                      "growthloci_2026-11-15/joint_sacks34/bmi/bmi_selected"]),
    ("bmi", "all", "all", ["bmi_allpanels_2026-11-16/G_all/bmi"]),
    ("height", "gtex_nonbrain", "noUnionDUP", ["robust_strata_2026-11-15/S1/height/height_selected"]),
    ("height", "gtex_nonbrain", "noAsym13", ["robust_strata_2026-11-15/S2/height/height_selected"]),
    ("height", "gtex_nonbrain", "noUnionDUP_noAsym13", ["robust_strata_2026-11-15/S3/height/height_selected"]),
    ("height", "gtex_nonbrain", "noUnionAny", ["union_strata_2026-11-15/joint_noUnionAny/height/height_selected"]),
    ("bmi", "all", "noUnionDUP", ["bmi_allpanels_2026-11-16/S1/bmi", "bmi_allpanels_union_2026-11-16/S1/bmi"]),
    ("bmi", "all", "noAsym13", ["bmi_allpanels_2026-11-16/S2/bmi", "bmi_allpanels_union_2026-11-16/S2/bmi"]),
    ("bmi", "all", "noUnionDUP_noAsym13", ["bmi_allpanels_2026-11-16/S3/bmi", "bmi_allpanels_union_2026-11-16/S3/bmi"]),
    ("bmi", "all", "noUnionAny", ["bmi_allpanels_2026-11-16/noUnionAny/bmi"]),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--code", required=True, help="frozen code snapshot (conga_grid_qt.py, make_union_strata.py, designs/)")
    ap.add_argument("--out", required=True); ap.add_argument("--work", required=True)
    ap.add_argument("--python", default=sys.executable); ap.add_argument("--nrep", type=int, default=2000)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); os.makedirs(a.work, exist_ok=True)
    G = ["--exons", f"{R}/conga/hg38_refGene_20200708.exon.txt", "--gene_cov", f"{R}/conga/gene_covariates.tsv",
         "--clingen", f"{R}/conga/ClinGen_gene_curation_list_GRCh38.tsv"]
    sd = f"{a.work}/strata"
    subprocess.run([a.python, "make_union_strata.py", "--pheno", f"{P}/data/ukbb_pheno_std.tsv", "--cnv", f"{P}/data/ukbb_pgc_xd.cnv", *G,
                    "--asym_loci", f"{R}/growth_loci/asymmetric_molly13.tsv", "--out", sd], check=True, cwd=a.code)
    rows = []
    for trait, panels, subset, refs in SPECS:
        pheno = f"{P}/data/ukbb_pheno_std.tsv" if subset == "all" else f"{sd}/pheno_{subset}.tsv"
        npy = f"{a.out}/{trait}_{panels}_{subset}.npy"
        if os.path.exists(npy):
            os.remove(npy)                                   # always recompute; never verify a stale file
        subprocess.run([a.python, "conga_grid_qt.py", "--pheno", pheno, "--trait", f"{trait}_rint", "--cnv", f"{P}/data/ukbb_pgc_xd.cnv",
                        "--twas", f"{R}/twas/{trait}.txt", "--panels", panels, *G, "--loci", f"{R}/growth_loci/loci_sacks34.tsv",
                        "--loo_runs", f"{R}/growth_loci/loo_runs_sacks34.tsv", "--design", "designs/grid_univariate.json", "--only", "none_H|adj",
                        "--nrep", str(a.nrep), "--seed", "2026", "--perms_file", npy, "--perms_only", "--out", f"{a.work}/info/{trait}_{panels}_{subset}"],
                       check=True, cwd=a.code)
        offs = np.load(npy).astype(np.int64)
        for ref in refs:
            f = f"{RES}/{ref}/run_info.json"
            if not os.path.exists(f):
                rows.append(dict(file=os.path.basename(npy), trait=trait, panels=panels, subset=subset, ref=ref, ref_N="", ref_nrep="",
                                 ref_sha1="", regenerated_sha1="", match="ref run not found")); continue
            info = json.load(open(f)); n = int(info.get("nrep", a.nrep)); ref_sha = info.get("perms_sha1", "")
            sha = hashlib.sha1(offs[:n].tobytes()).hexdigest()[:16]
            rows.append(dict(file=os.path.basename(npy), trait=trait, panels=panels, subset=subset, ref=ref, ref_N=info.get("N"), ref_nrep=n,
                             ref_sha1=ref_sha, regenerated_sha1=sha, match=("yes" if sha == ref_sha else ("no sha1 recorded" if not ref_sha else "NO"))))
        print(pd.DataFrame(rows).tail(len(refs)).to_string(index=False), flush=True)
    M = pd.DataFrame(rows); M.to_csv(f"{a.out}/SHIFTS_MANIFEST.tsv", sep="\t", index=False)
    for f in os.listdir(sd):
        if f.startswith("pheno_"):
            os.remove(os.path.join(sd, f))                  # person-level subsets; rebuilt by make_union_strata.py when needed
    bad = M[M.match == "NO"]
    print(M.to_string(index=False))
    if len(bad):
        raise SystemExit(f"{len(bad)} reference runs do NOT match their regenerated shifts; do not ship these files")


if __name__ == "__main__":
    main()
