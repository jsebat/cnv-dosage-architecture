"""Lay out the existing UK Biobank results behind the Sacks 2026 figures in the same folder structure that run_pipeline.py writes for
All of Us, then export them with the same export_aggregate.py (min_cell 0: the AoU small-cell rule does not apply to UKBB counts). The UKBB and AoU exports then have identical form and the meta-analysis
(meta/meta_sacks2026.py) reads both the same way.

Nothing is recomputed: the results tree is built from per-file symlinks to the original UKBB result folders (no extra disk space; the
project file system is shared and nearly full). Only the export folder holds copies, and it holds aggregates only.

  python assemble_ukbb_reference.py --code <cnv_anthro dir> --out /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/results/ukbb_reference_sacks2026
Writes <out>/results/all/... (symlinks), <out>/export/all/... (aggregate copies, EXPORT_MANIFEST.txt) and <out>/ASSEMBLY_MANIFEST.tsv
(which UKBB folder each piece came from, and which pieces are missing).
"""
import argparse
import json
import os
import subprocess
import sys
import pandas as pd

RES = "/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/results"
SKIP_DIRS = {"code_snapshot", "smoke", "strata", "work", "logs"}


def safe(n):
    return n.replace(":", "_").replace("|", "__").replace("~", "-")


def plan():
    """(destination under results/all, source folder under RES) for every UKBB piece the figures use."""
    P = [("locus_effects", "export_ukbb_strata/all/locus_effects")]           # CNV GWAS, S9 alleles, growth loci: the UKBB stratified export
    for t in ["height", "bmi"]:
        P.append((f"taskA_gwas/{t}", f"export_ukbb_strata/all/taskA_gwas/{t}"))
        P.append((f"taskA_alleles/{t}", f"export_ukbb_strata/all/taskA_alleles/{t}"))
    for t in ["height", "bmi"]:
        P.append((f"burden/{t}", f"dup_burden_2026-11-15/{t}"))                                         # Fig. 2A-B, 3O-P
        P.append((f"additivity/{t}/all", f"additivity_2026-11-15/{t}"))                                  # Fig. 2I-J
        for g in ["onlyasym", "onlymirror"]:
            for d in ["both", "DUP"]:
                P.append((f"additivity/{t}/{g}_{d}", f"additivity_{g}_2026-11-15/{t}/only_{d}"))         # Fig. 3H-I
        for h in ["top5_hi", "phaplo", "union"]:
            for ex in ["all", "noSacks34"]:
                if (h, ex) != ("top5_hi", "all"):
                    P.append((f"additivity/{t}/hit_{h}_{ex}", f"additivity_hit3_2026-11-15/{t}/{h}_{ex}"))  # model figure
        v = f"{t}_selected"
        P.append((f"scores/{t}/{v}", f"score_axes_2026-11-15/{t}/{v}"))                                  # Fig. 2L
        P.append((f"grid/joint_sacks34/{t}/{v}", f"growthloci_2026-11-15/joint_sacks34/{t}/{v}"))        # model figure F-G
        P.append((f"grid/joint_drivers/{t}/{v}", f"joint_drivers_2026-11-15/{t}/{v}"))                   # model figure F-G (minor_all)
        uni = json.load(open(os.path.join(CODE, "designs/grid_univariate_adj.json")))
        for var in uni:
            P.append((f"grid/grid_univariate/{t}/{safe(var['name'])}", f"growthloci_2026-11-15/grid/{t}/{safe(var['name'])}"))   # Fig. 2G, 2K
        for var in json.load(open(os.path.join(CODE, "designs/grid_callsize.json"))):
            P.append((f"grid/grid_callsize/{t}/{safe(var['name'])}", f"growthloci_2026-11-15/callsize/{t}/{safe(var['name'])}"))  # Fig. 2H
    # joint-model heatmap, UK Biobank columns (fig_joint_model_heatmap.R)
    J = [("height", "all", "height_selected_union", "union_terms_2026-11-15/height"),
         ("height", "noUnionDUP", "height_selected", "robust_strata_2026-11-15/S1/height/height_selected"),
         ("height", "noUnionAny", "height_selected", "union_strata_2026-11-15/joint_noUnionAny/height/height_selected"),
         ("height", "noAsym13", "height_selected", "robust_strata_2026-11-15/S2/height/height_selected"),
         ("height", "noUnionDUP_noAsym13", "height_selected", "robust_strata_2026-11-15/S3/height/height_selected"),
         ("bmi", "all", "bmi_selected_union", "bmi_allpanels_2026-11-16/G_all/bmi"),
         ("bmi", "noUnionDUP", "bmi_selected_union_noDUPunion", "bmi_allpanels_union_2026-11-16/S1/bmi"),
         ("bmi", "noUnionAny", "bmi_selected", "bmi_allpanels_2026-11-16/noUnionAny/bmi"),
         ("bmi", "noAsym13", "bmi_selected_union", "bmi_allpanels_union_2026-11-16/S2/bmi"),
         ("bmi", "noUnionDUP_noAsym13", "bmi_selected_union_noDUPunion", "bmi_allpanels_union_2026-11-16/S3/bmi")]
    for t, sub, var, src in J:
        P.append((f"joint/{sub}/{t}/{safe(var)}", src))
    # BMI univariate grid on all 45 panels: only once the UKBB run exists (folder name set by its job)
    for var in json.load(open(os.path.join(CODE, "designs/grid_univariate_adj.json"))):
        P.append((f"grid/grid_univariate_bmi_allpanels/bmi/{safe(var['name'])}", f"grid_univariate_bmi_allpanels_2026-11-16/bmi/{safe(var['name'])}"))
    return P


def main():
    global CODE
    ap = argparse.ArgumentParser()
    ap.add_argument("--code", required=True); ap.add_argument("--out", required=True); ap.add_argument("--python", default=sys.executable)
    a = ap.parse_args(); CODE = a.code
    root = os.path.join(a.out, "results", "all"); rows = []
    for dst, src in plan():
        s = os.path.join(RES, src); d = os.path.join(root, dst)
        files = [f for f in os.listdir(s) if os.path.isfile(os.path.join(s, f))] if os.path.isdir(s) else []
        files = [f for f in files if f != "null_full.npz" and not f.endswith(".log")]
        rows.append(dict(destination=dst, source=src, files=len(files), status="ok" if files else "MISSING"))
        if not files:
            continue
        os.makedirs(d, exist_ok=True)
        for f in files:
            lk = os.path.join(d, f)
            if os.path.lexists(lk):
                os.remove(lk)
            os.symlink(os.path.join(s, f), lk)
    M = pd.DataFrame(rows); M.to_csv(os.path.join(a.out, "ASSEMBLY_MANIFEST.tsv"), sep="\t", index=False)
    print(M.status.value_counts().to_string()); print(M[M.status != "ok"].to_string(index=False))
    subprocess.run([a.python, "export_aggregate.py", "--results", os.path.join(a.out, "results"), "--export", os.path.join(a.out, "export"),
                    "--min_cell", "0", "--strata", "all", "--partitions", "{}"], check=True, cwd=a.code)


if __name__ == "__main__":
    main()
