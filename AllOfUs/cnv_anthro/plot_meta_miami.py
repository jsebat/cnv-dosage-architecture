"""Miami plot (DUP up, DEL down) of a probe-level CNV meta-analysis (meta_sumstats.py output), same layout as plots.py.

Probes: the reference cohort's grid (meta_sumstats.py --ref_label). Points: fixed-effect (default) or random-effects P; probes
tested in only one cohort grey-blue; probes with heterogeneity P_Q < --het_p ringed in orange. Genome-wide line: 0.05 / the
reference cohort's independent-test count (its run_info.json), since the meta has no carrier sets of its own to merge.

  python plot_meta_miami.py --meta_dir meta/height_all --ref_run_info <UKBB export>/all/taskA_gwas/height/run_info.json \\
      --title "Height, UKBB + AoU" --out meta/height_all/miami_meta_height.png
"""
import argparse
import json
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator

GRCH38 = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
          114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468]
OFFS = pd.Series(np.concatenate([[0], np.cumsum(GRCH38[:-1])]), index=range(1, 23))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta_dir", required=True); ap.add_argument("--ref_run_info", required=True)
    ap.add_argument("--effects", default="FE", choices=["FE", "RE"]); ap.add_argument("--het_p", type=float, default=1e-3)
    ap.add_argument("--title", default=""); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    info = json.load(open(a.ref_run_info))
    fig, ax = plt.subplots(figsize=(11, 5.2)); offs = OFFS; ntot = {}
    for t, sgn in [("dup", 1), ("del", -1)]:
        S = pd.read_csv(os.path.join(a.meta_dir, f"meta_probes_{t}.tsv.gz"), sep="\t")
        S["CHR"] = pd.to_numeric(S.CHR, errors="coerce"); S["BP"] = pd.to_numeric(S.BP, errors="coerce"); S = S.dropna(subset=["CHR", "BP"])
        S = S[S.CHR <= 22]; ntot[t] = len(S)
        x = S.BP + S.CHR.map(offs); p = S[f"P_{a.effects}"].clip(lower=1e-300); lp = -np.log10(p) * sgn
        one = S.k < 2
        col = np.where(S.CHR % 2 == 0, "#7a7a7a", "#1f1f1f")
        ax.scatter(x[~one], lp[~one], s=4, c=col[(~one).to_numpy()], lw=0, rasterized=True)
        ax.scatter(x[one], lp[one], s=3, c="#9db4cc", lw=0, rasterized=True)
        het = (S.P_Q < a.het_p).fillna(False)
        ax.scatter(x[het], lp[het], s=14, facecolors="none", edgecolors="#e67e22", lw=0.6)
        thr = info.get(f"{t.upper()}_gws_threshold")
        if thr and not isinstance(thr, str):
            ax.axhline(-np.log10(thr) * sgn, color="#c0392b", lw=0.8, ls="--")
    mids = OFFS + np.array(GRCH38) / 2
    ax.set_xticks(mids.values); ax.set_xticklabels([str(c) for c in OFFS.index], fontsize=7); ax.set_xlim(0, OFFS.iloc[-1] + GRCH38[-1])
    yt = ax.get_yticks(); ax.yaxis.set_major_locator(FixedLocator(yt)); ax.set_yticklabels([f"{abs(v):g}" for v in yt]); ax.axhline(0, color="black", lw=0.6)
    ax.set_ylabel(f"-log10 P ({a.effects} meta)   (DUP above, DEL below)")
    ax.set_title(a.title + f"   probes DUP {ntot['dup']:,}, DEL {ntot['del']:,}; threshold from the reference cohort's independent tests", fontsize=9.5)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(a.out, dpi=200); print("wrote", a.out)


if __name__ == "__main__":
    main()
