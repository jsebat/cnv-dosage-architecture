# ConGA supplementary figure: per context (signal group, cell type, tissue, brain region) the DEL and DUP tiles and the two dosage
# readouts of the lab's dosage analyses - mirror, (z_DUP - z_DEL)/sqrt2, and dose-aligned, (z_DEL + z_DUP)/sqrt2 - each in LC and
# in HC genes, against the genome-shift null, with BH FDR and locus-robust marks (12 meta loci + minor_all + no_34). Also writes the
# signal-group table and the mirror-direction table (DUP up / DEL down vs DUP down / DEL up).
# null model: genome shift of the gene annotation (2,000 shared shifts, the same replicate for DEL and DUP, so each contrast has its
#   own null distribution); matched nulls in every locus-removal run
"""Readouts as in AllOfUs/cnv_anthro (conga_grid_qt.py readouts: mirror = (t_DUP - t_DEL)/sqrt2; the dosage paper's dose-aligned
contrast) and the ConGA paper's fig4_dosage_contrast.R. Here the two arms come from separate score tests (conga_ctx_qt.py; each model
keeps the person's total number of genes hit of that CNV type), combined per null replicate.
Stage A: statistics for every set, figure rows, and the list of sets that need locus-removal runs (shown units with any tile FDR < 0.10).
Stage B: adds the removal runs, writes the table and the figure.

  python conga_contrasts.py --stage A --ctx_run <conga_ctx run> --grp_run <group full run> --groups <signal groups dir> --out <dir>
  python conga_contrasts.py --stage B ... --loo_run <loo run>
"""
import argparse
import glob
import json
import os
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control as fdc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

S2 = np.sqrt(2)
TILES = [(c, s) for c in ("LC", "HC") for s in ("DEL", "DUP", "mirror", "aligned")]
TN = {"height": "Height", "bmi": "BMI"}
LIMB = {"B": "limb B cells", "PrehyperChon": "prehypertrophic chondrocytes", "DistalMes": "distal limb mesenchyme",
        "ArtiChon": "articular chondrocytes", "RestingChon": "resting chondrocytes", "ProlifChon": "proliferating chondrocytes",
        "OsteoB": "osteoblasts", "ChondroProg": "chondrocyte progenitors", "PeriChon": "perichondrium", "OCP": "osteochondral progenitors"}


def category(sid):
    if sid.startswith("SG:"):
        return "signal group"
    if sid.startswith("gtex:"):
        return "tissue"
    if sid.startswith("brain:BSpan"):
        return "brain region"
    return "cell type"


def unit_of(sid):
    return sid.rsplit("|", 1)[0]


def load1(run, t, d, mode):
    p = f"{run}/{mode}_{t}_{d}"
    D = pd.read_csv(f"{p}/conga_ctx_{mode}.tsv.gz", sep="\t")
    nul = {r: np.load(f) for r in D.run.unique() for f in [f"{p}/null_z_{r}.npy"] if os.path.exists(f)}
    return D, nul


def load(runs, t, d, mode):
    """one run dir, or several comma-separated (removal runs for extra sets): per run, rows and null columns are concatenated
    in the same order"""
    parts = [load1(r, t, d, mode) for r in runs.split(",") if r]
    if len(parts) == 1:
        return parts[0]
    rn = list(parts[0][1])
    D = pd.concat([x[x.run == r] for r in rn for x, _ in parts], ignore_index=True)
    nul = {r: np.hstack([n[r] for _, n in parts]) for r in rn}
    return D, nul


def stats_for(zdel, zdup, ndel, ndup):
    """obs arrays (K) and null arrays (R x K) -> dict stat -> (obs, null)"""
    return {"DEL": (zdel, ndel), "DUP": (zdup, ndup), "mirror": ((zdup - zdel) / S2, (ndup - ndel) / S2),
            "aligned": ((zdel + zdup) / S2, (ndel + ndup) / S2)}


def emp(obs, nul):
    ok = np.isfinite(nul)
    ge = ((np.abs(nul) >= np.abs(obs)[None, :] - 1e-12) & ok).sum(0); R = ok.sum(0)
    mu = np.nanmean(nul, 0); sd = np.nanstd(nul, 0)
    return (obs - mu) / np.where(sd > 0, sd, np.nan), (1 + ge) / (R + 1)


def full_table(runs):
    rows = []
    for t in ["height", "bmi"]:
        for run in runs:
            Dd, Nd = load(run, t, "DEL", "full"); Du, Nu = load(run, t, "DUP", "full")
            assert list(Dd.set_id) == list(Du.set_id)
            st = stats_for(Dd.z.to_numpy(), Du.z.to_numpy(), Nd["full"].astype(float), Nu["full"].astype(float))
            for s, (o, n) in st.items():
                ze, pe = emp(o, n)
                rows.append(pd.DataFrame(dict(trait=t, set_id=Dd.set_id, stat=s, obs=o, z_emp=ze, P_emp=pe, n_genes=Dd.n_genes)))
    T = pd.concat(rows, ignore_index=True)
    T["unit"] = T.set_id.map(unit_of); T["constraint"] = T.set_id.str.rsplit("|", n=1).str[1]; T["category"] = T.set_id.map(category)
    T["FDR_BH"] = np.nan
    for _, ix in T.groupby(["trait", "stat", "constraint", "category"]).groups.items():
        T.loc[ix, "FDR_BH"] = fdc(T.loc[ix, "P_emp"].to_numpy(), method="bh")
    return T


def tile_table(T):
    """keep the 16 figure tiles: DEL, DUP, mirror, dose-aligned in LC and in HC, per trait"""
    keep = np.zeros(len(T), bool)
    for con, s in TILES:
        keep |= ((T.constraint == con) & (T.stat == s)).to_numpy()
    return T[keep].copy()


def choose_rows(Tt, n_per=10):
    best = Tt.groupby(["category", "unit"]).P_emp.min().reset_index()
    rows = [best[best.category == "signal group"]]
    for c in ["cell type", "tissue", "brain region"]:
        rows.append(best[best.category == c].sort_values("P_emp").head(n_per))
    return pd.concat(rows)


def label(u, names):
    if u.startswith("SG:"):
        _, t, k = u.split(":"); nm = names.get(u, "")
        nm = nm.replace("GO:", "").replace("R:", "").replace("KEGG:", "").title()
        return f"{TN[t]} SG{k}: {nm[:52]}"
    fam, rest = u.split(":", 1)
    rest = rest.replace("WangEtal_", "").replace("BSpan_", "").replace("_Z1", "").replace("_", " ")
    rest = LIMB.get(rest, rest)
    src = {"gtex": "GTEx", "hpa_sc": "HPA", "limb": "limb", "brain": "brain"}[fam]
    return f"{rest} ({src})"


def figure(Tt, rows, names, out, cohort="UK Biobank"):
    cols = [(t, con, s) for t in ["height", "bmi"] for con, s in TILES]
    R = rows.sort_values(["category", "P_emp"]).reset_index(drop=True)
    order = ["signal group", "cell type", "tissue", "brain region"]
    R["c"] = R.category.map({c: i for i, c in enumerate(order)}); R = R.sort_values(["c", "P_emp"]).reset_index(drop=True)
    Z = np.full((len(R), len(cols)), np.nan); star = np.full(Z.shape, "", object)
    idx = Tt.set_index(["trait", "unit", "constraint", "stat"])
    for i, r in R.iterrows():
        for j, (t, con, s) in enumerate(cols):
            k = (t, r.unit, con, s)
            if k in idx.index:
                x = idx.loc[k]; x = x.iloc[0] if isinstance(x, pd.DataFrame) else x
                Z[i, j] = x.z_emp
                star[i, j] = ("■" if bool(x.get("locus_robust", False)) else "*") if x.FDR_BH < 0.05 else ""
    h = max(6, 0.17 * len(R) + 2.6)
    fig, ax = plt.subplots(figsize=(12.5, min(h, 14.5)))
    lim = np.nanmax(np.abs(Z)) if np.isfinite(Z).any() else 1
    im = ax.imshow(Z, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim), interpolation="nearest")
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            if star[i, j]:
                ax.text(j, i, star[i, j], ha="center", va="center", fontsize=6 if star[i, j] == "*" else 5)
    ax.set_yticks(range(len(R))); ax.set_yticklabels([label(u, names) for u in R.unit], fontsize=5.3)
    for tl, c in zip(ax.get_yticklabels(), R.category):   # signal-group rows prominent
        if c == "signal group":
            tl.set_fontweight("bold"); tl.set_fontsize(5.8)
    ax.set_xticks(range(len(cols)))
    SL = {"mirror": "mirror (DUP−DEL)/√2", "aligned": "dose-aligned (DEL+DUP)/√2", "DEL": "DEL", "DUP": "DUP"}
    ax.set_xticklabels([con + " " + SL[s] for t, con, s in cols], fontsize=6, rotation=60, ha="right")
    nt = len(TILES)
    for x in [nt / 2 - 0.5, nt + nt / 2 - 0.5]:   # LC | HC within each trait
        ax.axvline(x, color="white", lw=1.5)
    ax.axvline(nt - 0.5, color="black", lw=1.2)
    for b in np.flatnonzero(np.diff(R.c.to_numpy())) + 0.5:
        ax.axhline(b, color="black", lw=0.8)
    ax.text(nt / 2 - 0.5, -1.2, "Height", ha="center", fontsize=9); ax.text(nt + nt / 2 - 0.5, -1.2, "BMI", ha="center", fontsize=9)
    ax.set_title(f"Person-level ConGA, {cohort}: mirror and dose-aligned readouts of CNV dosage in LC and HC genes", fontsize=9, pad=18)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01); cb.set_label("z vs genome-shift null (red: trait higher)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    fig.text(0.01, 0.005, "* BH FDR < 5% (within trait, readout, constraint stratum and category);  ■ also holds in all 14 locus-removal runs. "
             "Rows: every signal group (bold), then the 10 most significant cell types, tissues and brain regions.", fontsize=6)
    fig.subplots_adjust(left=0.31, right=0.94, top=0.94, bottom=0.15); fig.savefig(out, dpi=300); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=["A", "B"]); ap.add_argument("--ctx_run", required=True)
    ap.add_argument("--grp_run", default="", help="separate signal-group run (UKBB); empty when --ctx_run holds every set (meta_ctx.py output)"); ap.add_argument("--groups", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--loo_run", default=""); ap.add_argument("--cohort_label", default="UK Biobank")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    names = {}
    for t in ["height", "bmi"]:
        g = pd.read_csv(f"{a.groups}/{t}.groups.tsv", sep="\t"); names.update(dict(zip(g.set_group, g.name)))
    T = full_table([r for r in (a.ctx_run, a.grp_run) if r]); Tt = tile_table(T); rows = choose_rows(Tt)
    if a.stage == "A":
        u = set(rows.unit) & set(Tt[Tt.FDR_BH < 0.10].unit)
        sel = sorted({f"{x}|{c}" for x in u for c in ("LC", "HC")})
        pd.Series(sel).to_csv(f"{a.out}/loo_set_list.txt", index=False, header=False)
        T.to_csv(f"{a.out}/conga_contrasts_full.tsv.gz", sep="\t", index=False, float_format="%.5g")
        rows.to_csv(f"{a.out}/figure_rows.tsv", sep="\t", index=False)
        print("rows", len(rows), "sets for removal runs", len(sel)); return
    # stage B: locus-removal runs (both dosages, same selected sets)
    hold = []
    for t in ["height", "bmi"]:
        Dd, Nd = load(a.loo_run, t, "DEL", "loo"); Du, Nu = load(a.loo_run, t, "DUP", "loo")
        for rn in Dd.run.unique():
            xd = Dd[Dd.run == rn].reset_index(drop=True); xu = Du[Du.run == rn].reset_index(drop=True)
            assert list(xd.set_id) == list(xu.set_id)
            st = stats_for(xd.z.to_numpy(), xu.z.to_numpy(), Nd[rn].astype(float), Nu[rn].astype(float))
            for s, (o, n) in st.items():
                ze, pe = emp(o, n)
                hold.append(pd.DataFrame(dict(trait=t, set_id=xd.set_id, stat=s, run=rn, obs_run=o, z_emp_run=ze, P_emp_run=pe)))
    H = pd.concat(hold, ignore_index=True).merge(T[["trait", "set_id", "stat", "obs"]], on=["trait", "set_id", "stat"])
    H["holds"] = (np.sign(H.obs_run) == np.sign(H.obs)) & (H.P_emp_run < 0.05)
    nrun = H.run.nunique()
    S = H.groupby(["trait", "set_id", "stat"]).agg(n_runs=("run", "nunique"), n_hold=("holds", "sum"),
                                                   worst_run=("P_emp_run", lambda s: H.loc[s.idxmax(), "run"])).reset_index()
    T = T.merge(S, on=["trait", "set_id", "stat"], how="left")
    T["locus_robust"] = (T.FDR_BH < 0.05) & (T.n_runs == nrun) & (T.n_hold == nrun)
    T.to_csv(f"{a.out}/conga_contrasts_full.tsv.gz", sep="\t", index=False, float_format="%.5g")
    H.to_csv(f"{a.out}/conga_contrasts_removal_runs.tsv.gz", sep="\t", index=False, float_format="%.5g")
    Tt = tile_table(T); Tt.to_csv(f"{a.out}/conga_contrasts_tiles.tsv", sep="\t", index=False, float_format="%.4g")
    figure(Tt, rows, names, f"{a.out}/SuppFig_ConGA_contexts.png", a.cohort_label)
    # summary for the report
    summ = Tt[Tt.FDR_BH < 0.05].groupby(["trait", "constraint", "stat", "category"]).agg(n_fdr05=("set_id", "size"),
                                                                                         n_robust=("locus_robust", "sum")).reset_index()
    summ.to_csv(f"{a.out}/conga_contrasts_summary.tsv", sep="\t", index=False); print(summ.to_string())
    top = Tt[Tt.locus_robust].sort_values("P_emp").groupby(["trait", "constraint", "stat"]).head(5)
    top.to_csv(f"{a.out}/conga_contrasts_top_robust.tsv", sep="\t", index=False, float_format="%.4g")
    json.dump(dict(removal_runs=int(nrun)), open(f"{a.out}/info.json", "w"))
    tables(T, names, a.out)


def tables(T, names, out):
    """signal-group table (every group, both traits, all 16 tiles) and mirror-direction table (every context, LC and HC)"""
    T = T[T.constraint.isin(["LC", "HC"])].copy()
    T["name"] = T.unit.map(lambda u: label(u, names))
    W = T.pivot_table(index=["trait", "category", "unit", "name", "constraint"], columns="stat",
                      values=["z_emp", "P_emp", "FDR_BH", "locus_robust", "n_hold", "n_runs"], aggfunc="first")
    W.columns = [f"{s}_{v}" for v, s in W.columns]; W = W.reset_index()
    sg = W[W.category == "signal group"].copy()
    sg["group_trait"] = sg.unit.str.split(":").str[1]
    sg.sort_values(["trait", "group_trait", "unit", "constraint"]).to_csv(f"{out}/conga_signal_groups_all_tiles.tsv", sep="\t",
                                                                         index=False, float_format="%.4g")
    M = W[["trait", "category", "unit", "name", "constraint", "DEL_z_emp", "DUP_z_emp", "mirror_z_emp", "mirror_P_emp",
           "mirror_FDR_BH", "mirror_locus_robust", "mirror_n_hold", "mirror_n_runs"]].copy()
    M["direction"] = np.where(M.mirror_z_emp > 0, "DUP up / DEL down", "DUP down / DEL up")
    M["locus_robust"] = M.mirror_locus_robust.fillna(False).astype(bool).map({True: "yes", False: "no"})
    M["rank_in_direction"] = M.groupby(["trait", "direction"]).mirror_P_emp.rank(method="first")
    M.sort_values(["trait", "direction", "rank_in_direction"]).to_csv(f"{out}/conga_mirror_direction.tsv", sep="\t", index=False,
                                                                     float_format="%.4g")


if __name__ == "__main__":
    main()
