"""Copy aggregate results to an export folder, enforcing a small-cell rule (All of Us: no count 1..20 may leave the
Workbench, nor anything that lets one be recomputed). Applied to every cohort so UKBB and AoU exports are identical in form.

Primary rule: in every per-probe / per-allele / per-locus table, a row whose carrier count is 1..min_cell loses its count
('<=min_cell') and its effect statistics (BETA/SE/T/P/FREQ...: 'suppressed'), because SE ~ 1/sqrt(n) would reveal the count.
Secondary rule (stratified runs): for each partition of the parent stratum 'all' (e.g. sex = female + male; ancestry =
eur + afr + ...), if for a row key (probe / allele / locus x dosage) any child stratum OR the remainder
(all - sum of children, e.g. participants without an ancestry call) has 1..min_cell carriers, that row is suppressed in
every child stratum of the partition (the parent row follows the primary rule only). Then no small count can be
recovered by subtraction across exported tables.
Per-probe association tables (sumstats_*, lead_*; Jon 2026-11-25: "export all probe level p-values so we can make plots regardless of raw
counts"): a withheld probe row is KEPT with its test statistic T and P (the signed T gives the direction); only the count, BETA, SE and FREQ
are withheld. T and P alone do not give the count back; BETA with SE would (SE ~ 1/sqrt(n)). --drop_small_probe_rows restores the old
behaviour (withheld probe rows left out).
Folders with their own exporter (fig3/: run_fig3.py; gene_assoc*/: export_gene_assoc.py) are skipped. conga_grid_qt.py folders (grid,
joint): a term with n_carriers_term 1..min_cell gets '<=min_cell', and beta/se/t of that (run, term) and its mirror readout are blanked
in observed_terms.tsv and in the null_runs.npz columns, as run_fig3.export_fit does (Milad 2026-11-25).
Person-level files (*genotype*, pheno, .cnv) are never copied. Null-replicate arrays (.npz, and the ctx step's null_z_*.npy; test
statistics only) are copied as float32. ConGA ctx tables (conga_ctx_*.tsv.gz) are copied without their per-set carrier counts.
"""
import argparse
import json
import os
import shutil
from collections import defaultdict
import numpy as np
import pandas as pd

COUNT_COLS = ["NCNVs", "n_carriers", "n", "cases", "ctrl", "carriers_port", "carriers_omar",   # person counts (n_probes is not one)
              "dup_carriers", "n_single", "n_multi", "carriers_nonzero"]                      # Sacks 2026 figure steps (burden, additivity, scores)
STAT_COLS = ["BETA", "SE", "T", "P", "BETA_SD", "raw_mean_diff_unadjusted", "FREQ", "median_call_kb",
             "mean_SD", "beta_SD", "beta", "beta_per_SD", "partial_R2", "slope_single", "SE_single", "slope_multi_minus_single", "SE_diff", "P_interaction"]
NEVER = ("genotype", "pheno", ".cnv", "allele_genotypes", "stratum_pheno")
KEYS = [["CHR", "BP"], ["trait", "allele"], ["allele"], ["trait", "NAME", "CNV"], ["group", "bin"], ["group", "dup_genes"],
        ["dosage", "group", "HC_bin"], ["axis", "model"], ["axis"]]


OWN_EXPORTERS = ("fig3", "gene_assoc", "gene_assoc_rare")   # folders exported by run_fig3.py / export_gene_assoc.py, never walked here


def small_terms(folder, k):
    """(run, readout) cells of a conga_grid_qt.py folder whose term has 1..k carriers (n_carriers_term), plus each one's mirror readout
    (Milad 2026-11-25: n_carriers_term was not masked; same rule as run_fig3.export_fit)."""
    p = os.path.join(folder, "observed_terms.tsv")
    if k <= 0 or not os.path.exists(p):
        return set()
    ob = pd.read_csv(p, sep="\t")
    if "n_carriers_term" not in ob.columns:
        return set()
    n = pd.to_numeric(ob.n_carriers_term, errors="coerce"); s = set()
    for r in ob[n.between(1, k)].itertuples():
        s.add((r.run, r.readout)); s.add((r.run, r.readout.rsplit("|", 1)[0] + "|mirror"))
    return s


def scrub_json(o, k):
    if isinstance(o, dict):
        return {kk: scrub_json(v, k) for kk, v in o.items()}
    if isinstance(o, list):
        return [scrub_json(v, k) for v in o]
    if isinstance(o, (int, np.integer)) and not isinstance(o, bool) and 1 <= o <= k:
        return f"<={k}"
    return o


def count_col(d):
    for c in COUNT_COLS:
        if c in d.columns:
            return c
    return None


def key_cols(d):
    for k in KEYS:
        if all(c in d.columns for c in k):
            return k
    return None


def read_tsv(p):
    return pd.read_csv(p, sep="\t")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True); ap.add_argument("--export", required=True)
    ap.add_argument("--min_cell", type=int, default=20)
    ap.add_argument("--strata", default="", help="comma list of stratum folders under --results (stratified runs)")
    ap.add_argument("--drop_small_probe_rows", action="store_true",
                    help="leave withheld probe rows out (default since 2026-11-25: keep them with T and P only)")
    ap.add_argument("--partitions", default='{"sex": ["female", "male"], "ancestry": ["eur", "afr", "amr", "eas", "sas", "mid"]}',
                    help="JSON: name -> child strata of the parent stratum 'all'")
    a = ap.parse_args()
    k = a.min_cell
    strata = [s for s in a.strata.split(",") if s]
    parts = {n: [c for c in ch if c in strata] for n, ch in json.loads(a.partitions).items()} if strata else {}
    # ---- secondary suppression map: (stratum, rel path within stratum) -> set of row keys to suppress
    extra = defaultdict(set)
    if strata and "all" in strata:
        rels = set()
        for root, dirs, files in os.walk(os.path.join(a.results, "all")):
            dirs[:] = [x for x in dirs if x not in OWN_EXPORTERS]
            for f in files:
                if f.endswith((".tsv", ".tsv.gz")) and not any(x in f for x in NEVER):
                    rels.add(os.path.relpath(os.path.join(root, f), os.path.join(a.results, "all")))
        for rel in sorted(rels):
            par = read_tsv(os.path.join(a.results, "all", rel)); cc = count_col(par); kc = key_cols(par)
            if cc is None or kc is None:
                continue
            pc = par.set_index(kc)[cc].pipe(pd.to_numeric, errors="coerce")
            pc = pc[~pc.index.duplicated()]
            for pname, kids in parts.items():
                if not kids:
                    continue
                tot = pd.Series(0.0, index=pc.index); bad = pd.Series(False, index=pc.index)
                for kid in kids:
                    p = os.path.join(a.results, kid, rel)
                    if not os.path.exists(p):
                        continue
                    d = read_tsv(p)
                    v = d.set_index(kc)[cc].pipe(pd.to_numeric, errors="coerce"); v = v[~v.index.duplicated()].reindex(pc.index).fillna(0)
                    tot += v; bad |= v.between(1, k)
                rem = pc.fillna(0) - tot
                bad |= rem.between(1, k)
                keys = set(bad[bad].index)
                for kid in kids:
                    extra[(kid, rel)] |= keys
    n_sup = 0; n_sec = 0
    for root, dirs, files in os.walk(a.results):
        dirs[:] = [x for x in dirs if x not in OWN_EXPORTERS]           # fig3/ and gene tables have their own exporters
        for f in files:
            if any(x in f for x in NEVER):
                continue
            src = os.path.join(root, f); rel = os.path.relpath(src, a.results); dst = os.path.join(a.export, rel)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            if f.startswith("conga_ctx_") and f.endswith((".tsv", ".tsv.gz")):
                # person-level ConGA (ctx step): per-set carrier counts never leave (they differ between full and removal runs, so a
                # small count could be recovered by subtraction); z, z_emp, P_emp and the null summaries are test statistics (Jon 2026-11-25)
                read_tsv(src).drop(columns=["n_carriers"], errors="ignore").to_csv(dst, sep="\t", index=False)
                continue
            if f == "observed_terms.tsv":                           # conga_grid_qt.py (grid, joint): run_fig3.export_fit's rule
                d = read_tsv(src); sm = small_terms(root, k)
                if sm:
                    hit = np.array([(r, q) in sm for r, q in zip(d.run, d.readout)])
                    d.loc[hit, [c for c in ("beta", "se", "t") if c in d.columns]] = np.nan; n_sup += int(hit.sum())
                if "n_carriers_term" in d.columns and k > 0:
                    v = pd.to_numeric(d.n_carriers_term, errors="coerce"); d["n_carriers_term"] = d.n_carriers_term.astype(object)
                    d.loc[v.between(1, k).to_numpy(), "n_carriers_term"] = f"<={k}"
                d.to_csv(dst, sep="\t", index=False, float_format="%.10g")
                continue
            if f == "null_runs.npz":                                # null columns of withheld (run, term) cells blanked, as export_fit
                z = np.load(src); out = {n: (z[n].astype(np.float32) if z[n].dtype.kind == "f" else z[n]) for n in z.files}
                sm = small_terms(root, k)
                if sm:
                    nf = os.path.join(root, "null_full.npz")
                    if not os.path.exists(nf):
                        raise SystemExit(f"{root}: withheld terms but no null_full.npz to map readout names; export this folder with run_fig3.py")
                    names = {n[len("names__"):]: [str(s) for s in np.load(nf)[n]] for n in np.load(nf).files if n.startswith("names__")}
                    for (rn, ro) in sm:
                        for vn, nm in names.items():
                            key = f"{rn}__{vn}"
                            if key in out and ro in nm and out[key].dtype.kind == "f":
                                out[key][:, nm.index(ro)] = np.nan
                np.savez_compressed(dst, **out)
                continue
            if f.endswith((".tsv", ".tsv.gz")):
                d = read_tsv(src)
                cc = [c for c in COUNT_COLS if c in d.columns]
                small = np.zeros(len(d), bool)
                for c in cc:
                    small |= pd.to_numeric(d[c], errors="coerce").between(1, k).to_numpy()
                n_sup += int(small.sum())
                parts_rel = rel.split(os.sep, 1)
                if len(parts_rel) == 2 and os.path.basename(parts_rel[1]).startswith("lead_"):     # lead rows: use the probe table's key set
                    parts_rel = [parts_rel[0], os.path.join(os.path.dirname(parts_rel[1]),
                                                            os.path.basename(parts_rel[1]).replace("lead_", "sumstats_").replace(".tsv", ".tsv.gz"))]
                if strata and len(parts_rel) == 2 and (parts_rel[0], parts_rel[1]) in extra and key_cols(d):
                    keys = extra[(parts_rel[0], parts_rel[1])]
                    kc = key_cols(d)
                    sec = np.array([tuple(r) in keys if len(kc) > 1 else r[0] in keys for r in d[kc].itertuples(index=False)]) & ~small
                    n_sec += int(sec.sum()); small |= sec
                probe = f.startswith(("sumstats_", "lead_"))
                for c in [c for c in STAT_COLS if c in d.columns and not (probe and not a.drop_small_probe_rows and c in ("T", "P"))]:
                    d[c] = d[c].astype(object); d.loc[small, c] = "suppressed"
                for c in cc:
                    v = pd.to_numeric(d[c], errors="coerce"); d[c] = d[c].astype(object)
                    d.loc[(v.between(1, k)).to_numpy(), c] = f"<={k}"                   # primary
                    d.loc[small & (v > k).to_numpy(), c] = "suppressed"                 # secondary: count > k but withheld
                if f.startswith("sumstats_") and a.drop_small_probe_rows:
                    d = d[~small]                                   # old behaviour: withheld probes left out
                if probe and not a.drop_small_probe_rows:
                    d["stats_withheld"] = small.astype(int)          # 1 = count, BETA, SE, FREQ withheld; T and P shown
                d.to_csv(dst, sep="\t", index=False)
            elif f.endswith(".json"):
                json.dump(scrub_json(json.load(open(src)), k), open(dst, "w"), indent=1)
            elif f == "null_full.npz":
                continue                                            # same arrays as the 'full__' entries of null_runs.npz
            elif f.endswith(".npz"):
                z = np.load(src)
                np.savez_compressed(dst, **{n: (z[n].astype(np.float32) if z[n].dtype.kind == "f" else z[n]) for n in z.files})
            elif f.startswith("null_z_") and f.endswith(".npy"):   # ctx step: score z of every null replicate (test statistics only)
                z = np.load(src)
                if z.dtype.kind != "f":
                    raise SystemExit(f"{src}: expected a float array of null z")
                np.save(dst, z.astype(np.float32))
            elif f.endswith((".png", ".pdf", ".md", ".txt")):
                shutil.copy(src, dst)
    tot = 0; lines = []
    for root, _, files in os.walk(a.export):
        for f in sorted(files):
            p = os.path.join(root, f); sz = os.path.getsize(p); tot += sz; lines.append(f"{sz}\t{os.path.relpath(p, a.export)}")
    with open(os.path.join(a.export, "EXPORT_MANIFEST.txt"), "w") as fh:
        fh.write(f"small-cell rule: counts 1..{k} withheld; rows withheld primary {n_sup}, secondary (cross-stratum) {n_sec}\n"
                 f"total bytes {tot} ({tot / 1e6:.1f} MB) in {len(lines)} files\nbytes\tfile\n" + "\n".join(lines) + "\n")
    print(f"exported to {a.export}; rows suppressed (primary, count 1..{k}): {n_sup}; secondary (cross-stratum): {n_sec}; total {tot / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
