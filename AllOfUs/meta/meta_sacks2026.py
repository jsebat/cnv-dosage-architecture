"""UK Biobank + All of Us meta-analysis of every genome-wide result behind the Sacks 2026 figures, for height and BMI.

Inputs: the two aggregate exports, which have identical form (AoU: run_pipeline.py --steps export; UKBB: ukbb/assemble_ukbb_reference.py).
Output: figure-ready tables with the SAME folder names, file names and columns as the UKBB tables the figure scripts read
(UKBB_anthro_CNV/results_aggregate/<analysis>_<date>/...), so figures/ redraws every UK Biobank panel from the meta-analysis by changing
one root folder (GENOMEWIDE_ROOT). The joint-model heatmap reads <out>/joint/<carrier subset>/<trait>/grid_summary.tsv.
Each meta table also keeps the per-cohort values and the heterogeneity test.

  python meta_sacks2026.py --ukbb <UKBB export>/all --aou <AoU export>/all --out <meta dir> \
      [--ukbb_gwas <UKBB CNV-GWAS export>/all] [--code ../cnv_anthro]

What is combined, and how
  CNV GWAS, S9 alleles, growth loci   meta_sumstats.py (fixed and random effects, Cochran Q, I2) + plot_meta_miami.py
  burden bins, dose-response bins     inverse-variance (IVW) mean per bin; carriers summed (Fig. 2A-B, 2I-J, 3H, 3O)
  burden / additivity models          IVW per (group or model, term) (U-shape P, Fig. 3P slopes, model figure single-hit tests)
  sub-additivity, asym vs mirror      IVW of the 2+ CNVs x constrained-gene interaction in each class, then the class difference (Fig. 3I)
  axis-score additivity               IVW of each cohort's ratio (slope in 2+ CNV carriers / slope in single-CNV carriers), SE = SE_diff/slope
                                      (Fig. 2L); written as slope_single = 1, slope_multi_minus_single = ratio - 1, SE_diff = SE(ratio),
                                      so the figure code reads it unchanged
  univariate grid, call size, joint   meta_grid.py: sqrt(N)-weighted Stouffer of t, replicate by replicate under the shared UKBB shifts,
  models                              empirical z and P, held across the 24 single-unit removals, no_loci, minor_all (Fig. 2G, 2H, 2K,
                                      model figure E-G, joint-model heatmap). Rows whose cohorts did not share shifts are flagged.
"""
import argparse
import glob
import json
import os
import subprocess
import sys
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.normpath(os.path.join(HERE, "..", "cnv_anthro"))


def ivw(b, s):
    b = np.asarray(b, float); s = np.asarray(s, float); ok = np.isfinite(b) & np.isfinite(s) & (s > 0)
    b, s = b[ok], s[ok]; k = len(b)
    if k == 0:
        return dict(k=0, est=np.nan, se=np.nan, P=np.nan, Q=np.nan, P_Q=np.nan, I2=np.nan)
    w = 1 / s ** 2; bf = np.sum(w * b) / w.sum(); sf = np.sqrt(1 / w.sum())
    Q = float(np.sum(w * (b - bf) ** 2)) if k > 1 else 0.0
    return dict(k=k, est=bf, se=sf, P=2 * stats.norm.sf(abs(bf / sf)), Q=Q, P_Q=stats.chi2.sf(Q, k - 1) if k > 1 else np.nan,
                I2=max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0)


def rd(p):
    return pd.read_csv(p, sep="\t") if os.path.exists(p) else None


def num(x):
    return pd.to_numeric(x, errors="coerce")


def meta_rows(tabs, key, bcol, scol, ncols=(), keep=()):
    """tabs: {cohort: DataFrame}. IVW of bcol/scol per key; carriers in ncols summed (suppressed AoU counts '<=20' -> not added);
    per-cohort estimates kept as <bcol>_<cohort>; extra columns in keep are taken from the first cohort that has the row."""
    labs = [c for c in tabs if tabs[c] is not None]
    allk = pd.concat([tabs[c][key] for c in labs]).drop_duplicates()
    idx = {c: tabs[c].drop_duplicates(key).set_index(key) for c in labs}
    out = []
    for kv in allk.itertuples(index=False):
        kv = tuple(kv); row = dict(zip(key, kv)); b, s = [], []
        for c in labs:
            t = idx[c]; r = t.loc[kv if len(key) > 1 else kv[0]] if (kv if len(key) > 1 else kv[0]) in t.index else None
            bb = float(num(pd.Series([r[bcol]])).iloc[0]) if r is not None else np.nan
            ss = float(num(pd.Series([r[scol]])).iloc[0]) if r is not None else np.nan
            b.append(bb); s.append(ss); row[f"{bcol}_{c}"] = bb; row[f"{scol}_{c}"] = ss
            for nc in ncols:
                if r is not None and nc in r.index:
                    row[f"{nc}_{c}"] = r[nc]
            for kc in keep:
                if r is not None and kc in r.index and kc not in row:
                    row[kc] = r[kc]
        m = ivw(b, s)
        row.update({bcol: m["est"], scol: m["se"], "P": m["P"], "k": m["k"], "Q": m["Q"], "P_Q": m["P_Q"], "I2": m["I2"],
                    "cohorts": ";".join(c for c, bb, ss in zip(labs, b, s) if np.isfinite(bb) and np.isfinite(ss))})
        for nc in ncols:
            v = [num(pd.Series([row.get(f"{nc}_{c}")])).iloc[0] for c in labs]
            row[nc] = np.nansum(v) if np.isfinite(v).any() else np.nan
        out.append(row)
    return pd.DataFrame(out)


def wr(d, p):
    os.makedirs(os.path.dirname(p), exist_ok=True); d.to_csv(p, sep="\t", index=False, float_format="%.6g"); print("wrote", p, len(d))


def burden(C, out):
    for t in ["height", "bmi"]:
        B = {c: rd(f"{r}/burden/{t}/signed_burden_bins.tsv") for c, r in C.items()}
        if all(v is not None for v in B.values()):
            d = meta_rows(B, ["group", "bin"], "mean_SD", "SE", ncols=["n"], keep=["side", "order"])
            wr(d.sort_values(["group", "order"]), f"{out}/dup_burden_2026-11-15/{t}__signed_burden_bins.tsv")
        for f in ["signed_burden_models", "dup_burden_models"]:
            M = {c: rd(f"{r}/burden/{t}/{f}.tsv") for c, r in C.items()}
            key = ["group", "term"] if f == "signed_burden_models" else ["model", "term"]
            if all(v is not None for v in M.values()):
                wr(meta_rows(M, key, "beta_SD", "SE", ncols=["N"]), f"{out}/dup_burden_2026-11-15/{t}__{f}.tsv")


def additivity(C, out):
    for t in ["height", "bmi"]:
        for run, name in [("all", f"{t}__dose_response_bins.tsv")] + [(f"{g}_{d}", f"{t}__{g}_only_{d}__dose_response_bins.tsv")
                                                                       for g in ["onlyasym", "onlymirror"] for d in ["both", "DUP"]]:
            B = {c: rd(f"{r}/additivity/{t}/{run}/dose_response_bins.tsv") for c, r in C.items()}
            if all(v is not None for v in B.values()):
                wr(meta_rows(B, ["dosage", "group", "HC_bin"], "mean_SD", "SE", ncols=["n"]), f"{out}/additivity_2026-11-15/{name}")
        for run in sorted({os.path.basename(p) for r in C.values() for p in glob.glob(f"{r}/additivity/{t}/*")}):
            M = {c: rd(f"{r}/additivity/{t}/{run}/additivity_models.tsv") for c, r in C.items()}
            if all(v is not None for v in M.values()):
                d = meta_rows(M, ["model", "term"], "beta", "SE"); d["T"] = d.beta / d.SE
                tgt = f"{out}/additivity_2026-11-15/{t}__additivity_models.tsv" if run == "all" else \
                      f"{out}/additivity_hit3_2026-11-15/{t}/{run.replace('hit_', '')}/additivity_models.tsv" if run.startswith("hit_") else \
                      f"{out}/additivity_2026-11-15/{t}__{run}__additivity_models.tsv"
                d = d[["term", "beta", "SE", "T", "P", "model"] + [c for c in d.columns if c not in ("term", "beta", "SE", "T", "P", "model")]]
                wr(d, tgt)
                if run == "all":                                  # the 'all' run is also the top5_hi / all-CNVs run of the single-hit set
                    wr(d, f"{out}/additivity_hit3_2026-11-15/{t}/top5_hi_all/additivity_models.tsv")
    # Fig. 3I: sub-additivity in asymmetric vs mirror CNV carriers (same layout as jobs/compare_subadditivity.py)
    rows = []
    for t in ["height", "bmi"]:
        for dos in ["both", "DUP"]:
            est = {}
            for grp in ["onlyasym", "onlymirror"]:
                M = {c: rd(f"{r}/additivity/{t}/{grp}_{dos}/additivity_models.tsv") for c, r in C.items()}
                if not all(v is not None for v in M.values()):
                    continue
                d = meta_rows(M, ["model", "term"], "beta", "SE").set_index(["model", "term"])
                for mod, mn in [("linear HC + multi interactions", "linear"), ("concave HC + multi interactions", "concave")]:
                    for term in ["multi_x_HC_DEL", "multi_x_HC_DUP"]:
                        if (mod, term) in d.index:
                            x = d.loc[(mod, term)]; est[(grp, mn, term)] = (x.beta, x.SE, x.P, x.cohorts)
            for mn in ["linear", "concave"]:
                for term in ["multi_x_HC_DEL", "multi_x_HC_DUP"]:
                    a_, b_ = est.get(("onlyasym", mn, term)), est.get(("onlymirror", mn, term))
                    if a_ and b_:
                        z = (a_[0] - b_[0]) / np.sqrt(a_[1] ** 2 + b_[1] ** 2)
                        rows.append(dict(trait=t, carriers="both dosages" if dos == "both" else "DUP carriers", model=mn, term=term,
                                         asym=f"{a_[0]:+.3f} (SE {a_[1]:.3f}, P {a_[2]:.2g}; meta {a_[3]})",
                                         mirror=f"{b_[0]:+.3f} (SE {b_[1]:.3f}, P {b_[2]:.2g}; meta {b_[3]})",
                                         z_diff=round(float(z), 2), P_diff=float(2 * stats.norm.sf(abs(z)))))
    if rows:
        wr(pd.DataFrame(rows), f"{out}/additivity_2026-11-15/subadditivity_asym_vs_mirror_2026-11-15.tsv")


def scores(C, out):
    for t in ["height", "bmi"]:
        for v in sorted({os.path.basename(p) for r in C.values() for p in glob.glob(f"{r}/scores/{t}/*")}):
            A = {c: rd(f"{r}/scores/{t}/{v}/additivity.tsv") for c, r in C.items()}
            if not all(x is not None for x in A.values()):
                continue
            rows = []
            for ax in A[list(A)[0]].axis:
                row = dict(axis=ax); r_, s_ = [], []
                for c, d in A.items():
                    x = d[d.axis == ax]
                    if not len(x):
                        r_.append(np.nan); s_.append(np.nan); continue
                    x = x.iloc[0]; sl = float(num(pd.Series([x.slope_single])).iloc[0])
                    dd = float(num(pd.Series([x.slope_multi_minus_single])).iloc[0]); se = float(num(pd.Series([x.SE_diff])).iloc[0])
                    r_.append((sl + dd) / sl); s_.append(se / abs(sl))
                    row[f"ratio_{c}"] = r_[-1]; row[f"SE_ratio_{c}"] = s_[-1]
                    for nc in ["n_single", "n_multi"]:
                        row[f"{nc}_{c}"] = x.get(nc)
                m = ivw(r_, s_)
                row.update(slope_single=1.0, slope_multi_minus_single=m["est"] - 1, SE_diff=m["se"],
                           P_interaction=2 * stats.norm.sf(abs((m["est"] - 1) / m["se"])) if m["k"] else np.nan, ratio=m["est"], SE_ratio=m["se"],
                           k=m["k"], Q=m["Q"], P_Q=m["P_Q"], I2=m["I2"])
                for nc in ["n_single", "n_multi"]:
                    vals = [num(pd.Series([row.get(f"{nc}_{c}")])).iloc[0] for c in A]
                    row[nc] = np.nansum(vals)
                rows.append(row)
            wr(pd.DataFrame(rows), f"{out}/score_axes_2026-11-15/{t}_{v}__additivity.tsv")


def grid(C, out, code):
    """meta_grid.py per design x trait; adds trait / panels and restores the variant names (folders are 'safe' names)."""
    U, A = C["UKBB"], C["AoU"]
    jobs = []
    for des in ["grid_univariate", "grid_callsize", "grid_univariate_bmi_allpanels", "joint_sacks34"]:
        for t in ["height", "bmi"]:
            jobs.append((f"{U}/grid/{des}/{t}", f"{A}/grid/{des}/{t}", f"{out}/_grid/{des}_{t}.tsv", t))
    for t in ["height", "bmi"]:   # UKBB drivers run (8 majors + minor_all) against the AoU all_minor run: the minor_all readout
        jobs.append((f"{U}/grid/joint_drivers/{t}", f"{A}/grid/joint_sacks34/{t}", f"{out}/_grid/joint_drivers_{t}.tsv", t))
    for sub in ["all", "noUnionDUP", "noUnionAny", "noAsym13", "noUnionDUP_noAsym13"]:
        for t in ["height", "bmi"]:
            jobs.append((f"{U}/joint/{sub}/{t}", f"{A}/joint/{sub}/{t}", f"{out}/joint/{sub}/{t}/grid_summary.tsv", t))
    res = {}
    for u, a, o, t in jobs:
        if not (os.path.isdir(u) and os.path.isdir(a)):
            print("skip (missing in a cohort):", o); continue
        os.makedirs(os.path.dirname(o), exist_ok=True)
        subprocess.run([sys.executable, "meta_grid.py", "--cohort", f"UKBB={u}", "--cohort", f"AoU={a}", "--out", o], check=True, cwd=code)
        d = pd.read_csv(o, sep="\t") if os.path.getsize(o) > 1 else pd.DataFrame()
        if not len(d):
            print(f"WARNING {o}: no variant could be meta-analysed (missing null_runs.npz in a cohort)"); continue
        names = {}
        for vd in os.listdir(u):
            f = os.path.join(u, vd, "run_info.json")
            if os.path.exists(f):
                info = json.load(open(f)); names[vd] = info["variants"][0]
                d.loc[d.variant == vd, "panels"] = info.get("panels")
        d["variant"] = d.variant.map(lambda x: names.get(x, x)); d.insert(0, "trait", f"{t}_rint")
        if (~d.shared_shifts.astype(bool)).any():
            print(f"WARNING {o}: {(~d.shared_shifts.astype(bool)).sum()} rows without shared shifts (z_emp not replicate-matched)")
        d.to_csv(o, sep="\t", index=False, float_format="%.6g"); res[o] = d
    # figure-ready copies under the UKBB names
    def cat(prefix):
        fs = [f for f in res if os.path.basename(f).startswith(prefix)]
        return pd.concat([res[f] for f in fs]) if fs else None
    gl = [res[f] for f in res if os.path.basename(f) in ("grid_univariate_height.tsv", "grid_univariate_bmi.tsv")]
    G = pd.concat(gl) if gl else None
    if G is not None and len(G):
        G = G.copy(); G["variant"] = G.variant.astype(str)
        wr(G, f"{out}/growthloci_2026-11-15/grid_summary_all.tsv")
    CS = cat("grid_callsize")
    if CS is not None:
        wr(CS, f"{out}/growthloci_2026-11-15/callsize_summary_all.tsv")
    BA = cat("grid_univariate_bmi_allpanels")
    if BA is not None:
        wr(BA, f"{out}/growthloci_bmi_allpanels/grid_summary_all.tsv")
    for t in ["height", "bmi"]:
        f = f"{out}/_grid/joint_sacks34_{t}.tsv"
        if f in res:
            wr(res[f], f"{out}/joint_sacks34_2026-11-15/{t}_{t}_selected.tsv")
        f = f"{out}/_grid/joint_drivers_{t}.tsv"
        if f in res:
            wr(res[f], f"{out}/joint_drivers_2026-11-15/{t}/{t}_selected/grid_summary.tsv")


def gwas(U, A, out, code):
    for t in ["height", "bmi"]:
        o = f"{out}/cnv_gwas/{t}_all"
        subprocess.run([sys.executable, "meta_sumstats.py", "--cohort", f"UKBB={U}", "--cohort", f"AoU={A}", "--trait", t, "--out", o,
                        "--ref_label", "UKBB"], check=True, cwd=code)
        subprocess.run([sys.executable, "plot_meta_miami.py", "--meta_dir", o, "--ref_run_info", f"{U}/taskA_gwas/{t}/run_info.json",
                        "--title", f"{'Height' if t == 'height' else 'BMI'}, UK Biobank + All of Us", "--out", f"{o}/miami_meta_{t}.png"],
                       check=True, cwd=code)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ukbb", required=True, help="UKBB export, stratum folder (…/export/all)")
    ap.add_argument("--aou", required=True, help="AoU export, stratum folder (…/export_aggregate/all)")
    ap.add_argument("--ukbb_gwas", default=None, help="UKBB CNV-GWAS export stratum folder (taskA_gwas, taskA_alleles, locus_effects), if not in --ukbb")
    ap.add_argument("--out", required=True); ap.add_argument("--code", default=CODE)
    ap.add_argument("--parts", default="gwas,burden,additivity,scores,grid")
    a = ap.parse_args(); C = {"UKBB": a.ukbb, "AoU": a.aou}; parts = a.parts.split(",")
    os.makedirs(a.out, exist_ok=True)
    if "gwas" in parts:
        gwas(a.ukbb_gwas or a.ukbb, a.aou, a.out, a.code)
    if "burden" in parts:
        burden(C, a.out)
    if "additivity" in parts:
        additivity(C, a.out)
    if "scores" in parts:
        scores(C, a.out)
    if "grid" in parts:
        grid(C, a.out, a.code)
    json.dump(dict(ukbb=a.ukbb, aou=a.aou, ukbb_gwas=a.ukbb_gwas, parts=parts), open(f"{a.out}/meta_info.json", "w"), indent=1)


if __name__ == "__main__":
    main()
