"""Miami plot (DUP up, DEL down) of a probe-level CNV GWAS, in Omar's style: kept probes black/grey by chromosome,
filtered probes light red, genome-wide line at 0.05 / n_independent_tests per CNV type."""
import argparse
import json
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gwas_dir", required=True); ap.add_argument("--title", default=""); ap.add_argument("--out", required=True)
    ap.add_argument("--min_ncnv_plot", type=int, default=1, help="hide probes with fewer carriers (AoU small-cell rule)")
    a = ap.parse_args()
    info = json.load(open(os.path.join(a.gwas_dir, "run_info.json")))
    fig, ax = plt.subplots(figsize=(11, 5.2))
    offs = None
    for t, sgn in [("dup", 1), ("del", -1)]:
        S = pd.read_csv(os.path.join(a.gwas_dir, f"sumstats_{t}.tsv.gz"), sep="\t")
        S = S[S.NCNVs >= a.min_ncnv_plot]
        if offs is None:
            ln = S.groupby("CHR").BP.max()
            offs = pd.Series(np.concatenate([[0], np.cumsum(ln.values[:-1])]), index=ln.index)
        x = S.BP + S.CHR.map(offs)
        lp = -np.log10(S.P.clip(lower=1e-300)) * sgn
        f = ~S.KEEP
        ax.scatter(x[f], lp[f], s=3, c="#f4a6a6", lw=0, rasterized=True)
        col = np.where(S.CHR % 2 == 0, "#7a7a7a", "#1f1f1f")
        ax.scatter(x[~f], lp[~f], s=4, c=col[~f.to_numpy()], lw=0, rasterized=True)
        thr = info.get(f"{t.upper()}_gws_threshold")
        if thr:
            ax.axhline(-np.log10(thr) * sgn, color="#c0392b", lw=0.8, ls="--")
    mids = offs + pd.Series(offs).diff().shift(-1).fillna(offs.iloc[-1] * 0.02) / 2
    ax.set_xticks(mids.values); ax.set_xticklabels([("X" if c == 23 else str(c)) for c in offs.index], fontsize=7)
    yt = ax.get_yticks(); ax.set_yticklabels([f"{abs(v):g}" for v in yt])
    ax.axhline(0, color="black", lw=0.6)
    ax.set_ylabel("-log10 P   (DUP above, DEL below)")
    ax.set_title(a.title + f"   N = {info['N']:,}; tests DUP {info.get('DUP_independent_tests')}, DEL {info.get('DEL_independent_tests')}", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(a.out, dpi=200)


if __name__ == "__main__":
    main()
