"""Aggregate-only export of the gene-based association tables (gene_assoc_qt.py), with the small-cell rule of export_aggregate.py.

  python export_gene_assoc.py --results <results_dir> --export <new dir> --strata all,eur,afr,... --min_cell 20 [--min_carriers 21]
        [--partitions '{"sex": ["female", "male"], "ancestry": ["eur", "afr", ...]}']

Reads <results>/<stratum>/<set>/<trait>/gene_assoc.tsv.gz, writes <export>/<stratum>/<set>/<trait>/gene_assoc.tsv.gz + run_info.json, for
each table set of --sets that exists: gene_assoc (the cohort's CNV table, i.e. the PGC freeze filters) and the optional, off-by-default
gene_assoc_rare (calls that also pass a position-frequency filter).
max_pos_n (carrier count behind max_pos_freq) is never exported; max_pos_freq is blanked when that count is 1..min_cell.
Rule (rows are DROPPED, never written with masked counts):
  floor      a gene x dosage row needs >= --min_carriers carriers in that stratum: 21 in every cohort, UK Biobank and All of Us alike
             (Jon 2026-11-23: "ok 21"), for the main results and the meta-analysis;
  primary    rows with 1..min_cell carriers are dropped (All of Us: min_cell 20, so >= 21 carriers);
  secondary  for each partition of the parent stratum 'all' (sex; ancestry): if for a gene x dosage any child stratum, or the remainder
             (all - sum of the children, e.g. people without an ancestry call), has 1..min_cell carriers, that row is dropped from EVERY child
             stratum of the partition (the parent row follows the primary rule only), so no small count can be recovered by subtraction.
  min_cell 0 = no small-cell rule (UK Biobank at UCSD): only the floor applies.
The audit at the end re-reads the export and fails if any row has fewer than max(min_carriers, min_cell + 1) carriers.
null model: none (copies analytic OLS tables).
"""
import argparse
import json
import os
import numpy as np
import pandas as pd

KEY = ["entrez", "dosage"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True); ap.add_argument("--export", required=True); ap.add_argument("--strata", default="all")
    ap.add_argument("--min_cell", type=int, default=20); ap.add_argument("--min_carriers", type=int, default=21)
    ap.add_argument("--sets", default="gene_assoc,gene_assoc_rare")
    ap.add_argument("--partitions", default='{"sex": ["female", "male"], "ancestry": ["eur", "afr", "amr", "eas", "sas", "mid"]}')
    a = ap.parse_args(); k = a.min_cell; need = max(a.min_carriers, k + 1)
    if os.path.exists(a.export):
        raise SystemExit(f"{a.export} exists; move it away first (an export is never overwritten)")
    lines = []; tot = 0; summary = []
    jobs = []
    for name in [x for x in a.sets.split(",") if x]:
        strata = [s for s in a.strata.split(",") if os.path.isdir(os.path.join(a.results, s, name))]
        parts = {n: [c for c in ch if c in strata] for n, ch in json.loads(a.partitions).items()}
        jobs += [(name, strata, parts, t) for t in sorted({t for s in strata for t in os.listdir(os.path.join(a.results, s, name))})]
    for name, strata, parts, t in jobs:
        T = {}
        for s in strata:
            f = os.path.join(a.results, s, name, t, "gene_assoc.tsv.gz")
            if os.path.exists(f):
                T[s] = pd.read_csv(f, sep="\t", dtype={"entrez": str})
        drop = {s: set() for s in T}                                    # secondary suppression: keys dropped in the child strata
        if k > 0 and "all" in T:
            pc = T["all"].set_index(KEY).n_carriers
            for pname, kids in parts.items():
                kids = [c for c in kids if c in T]
                if not kids:
                    continue
                tot_k = pd.Series(0.0, index=pc.index); bad = pd.Series(False, index=pc.index)
                for c in kids:
                    v = T[c].set_index(KEY).n_carriers.reindex(pc.index).fillna(0); tot_k += v; bad |= v.between(1, k)
                bad |= (pc - tot_k).between(1, k)
                for c in kids:
                    drop[c] |= set(bad[bad].index)
        for s, d in T.items():
            n = d.n_carriers; keep = n >= need
            sec = np.array([key in drop[s] for key in zip(d.entrez, d.dosage)]) if drop[s] else np.zeros(len(d), bool)
            o = d[keep & ~sec].copy(); dst = os.path.join(a.export, s, name, t); os.makedirs(dst)
            if "max_pos_n" in o.columns:
                o.loc[o.max_pos_n.between(1, k), "max_pos_freq"] = np.nan; o = o.drop(columns="max_pos_n")
            o.to_csv(os.path.join(dst, "gene_assoc.tsv.gz"), sep="\t", index=False, float_format="%.6g")
            info = json.load(open(os.path.join(a.results, s, name, t, "run_info.json")))
            info = {kk: info[kk] for kk in ("cohort", "stratum", "trait", "N", "covariates", "panels", "n_panels", "genes", "test", "segdup_bed", "cnv_table", "cnv_filter") if kk in info}
            info.update(min_cell=k, min_carriers_exported=need, rows=int(len(o)), rows_dropped_below_limit=int((~keep).sum()),
                        rows_dropped_secondary=int((keep & sec).sum()))
            json.dump(info, open(os.path.join(dst, "run_info.json"), "w"), indent=1)
            summary.append((name, s, t, len(d), len(o), int((~keep).sum()), int((keep & sec).sum())))
    bad = []
    for root, _, files in os.walk(a.export):
        for f in sorted(files):
            p = os.path.join(root, f); sz = os.path.getsize(p); tot += sz; lines.append(f"{sz}\t{os.path.relpath(p, a.export)}")
            if f.endswith(".tsv.gz"):
                x = pd.read_csv(p, sep="\t")
                if (x.n_carriers < need).any() or "max_pos_n" in x.columns:
                    bad.append(p)
            if any(x in f for x in ("pheno", "genotype", ".cnv")):
                bad.append(p)
    with open(os.path.join(a.export, "EXPORT_MANIFEST.txt"), "w") as fh:
        fh.write(f"gene-based association export; rows need >= {need} carriers (min_carriers {a.min_carriers}, small-cell rule 1..{k} withheld, rows dropped)\n"
                 f"audit: {'OK' if not bad else 'FAILED'}\ntotal bytes {tot} ({tot / 1e6:.1f} MB) in {len(lines)} files\nbytes\tfile\n" + "\n".join(lines) + "\n")
    print("set\tstratum\ttrait\trows_in\trows_exported\tdropped_below_limit\tdropped_secondary")
    for r in summary:
        print("\t".join(map(str, r)))
    print(f"exported to {a.export}: {tot / 1e6:.1f} MB")
    print("\n".join("AUDIT VIOLATION " + b for b in bad) or f"AUDIT OK: every exported row has >= {need} carriers; no person-level file")
    if bad:
        raise SystemExit("export audit failed: do not download this folder")


if __name__ == "__main__":
    main()
