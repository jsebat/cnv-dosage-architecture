"""Aggregate-only export of the gene-based association tables (gene_assoc_qt.py), with the small-cell rule of export_aggregate.py.

  python export_gene_assoc.py --results <results_dir> --export <new dir> --strata all,eur,afr,... --min_cell 20 [--min_carriers 1]
        [--partitions '{"sex": ["female", "male"], "ancestry": ["eur", "afr", ...]}']

Reads <results>/<stratum>/<set>/<trait>/gene_assoc.tsv.gz, writes <export>/<stratum>/<set>/<trait>/gene_assoc.tsv.gz + run_info.json, for
each table set of --sets that exists: gene_assoc (the cohort's CNV table, i.e. the PGC freeze filters) and the optional, off-by-default
gene_assoc_rare (calls that also pass a position-frequency filter).
Classes: DEL, DUP, and with gene_assoc_qt.py --plof also PLOF (pLoF SNV) and LOF (DEL or pLoF SNV).

Rule (Jon 2026-11-25: "Don't exclude gene based association p-values of any kind. but we should be careful with raw counts data"):
  every gene x class row is exported with its t and P. On a row whose count is withheld (n_carriers = NA, n_carriers_masked = 1) beta, SE,
  freq and max_pos_freq are blanked too: freq x N gives the count back, and so does SE (SE ~ 1/sqrt(carriers) for a rank-normalised trait).
  The signed t carries the direction; meta_gene_assoc.py combines such rows by a sample-size-weighted Z. n_carriers_masked: 1 = count
  1..min_cell (freq_max = (min_cell + 1) / N is then an upper bound on the carrier frequency), 2 = count above min_cell withheld only because
  of a subtraction. Variance explained needs no count: partial R2 = t^2 / (t^2 + df), noise-corrected (t^2 - 1) / df, df ~ N.
  floor      --min_carriers (default 1) only removes rows with fewer carriers than a usable estimate needs; it is not a reporting rule
  primary    counts 1..min_cell are withheld (All of Us: min_cell 20)
  secondary  (strata) for each partition of the parent stratum 'all' (sex; ancestry): if for a gene x class any child stratum, or the
             remainder (all - sum of the children), has 1..min_cell carriers, that count is withheld in EVERY child stratum of the partition
  classes    within a gene and stratum, the shown counts of DEL, PLOF and LOF must not give a count 1..min_cell by subtraction:
             LOF - DEL (pLoF-only carriers), LOF - PLOF (deletion-only carriers) and DEL + PLOF - LOF (carriers of both). Counts are shown in the
             priority DEL, LOF, PLOF; a lower-priority count is withheld if, together with the counts already shown, it would give one.
  min_cell 0 = no small-cell rule (UK Biobank at UCSD): nothing is withheld; n_overlap is dropped in every case.
The audit at the end re-reads the export and fails if any shown count is 1..min_cell or a person-level file is present.
null model: none (copies analytic OLS tables).
"""
import argparse
import json
import os
import numpy as np
import pandas as pd

KEY = ["entrez", "dosage"]


def small(x, k):
    return (x >= 1) & (x <= k)


def class_mask(d, k):
    """Within-gene cross-class rule; returns a boolean Series (True = withhold the count) aligned to d."""
    m = pd.Series(False, index=d.index)
    if k <= 0 or "LOF" not in set(d.dosage):
        return m
    w = d.pivot_table(index="entrez", columns="dosage", values="n_carriers", aggfunc="first")
    for c in ("DEL", "PLOF", "LOF"):
        if c not in w.columns:
            w[c] = 0
    w = w.fillna(0)
    ov = d[d.dosage == "LOF"].set_index("entrez").n_overlap.reindex(w.index).fillna(0) if "n_overlap" in d.columns else pd.Series(0, index=w.index)
    shown = {c: ~small(w[c], k) & (w[c] > 0) for c in ("DEL", "LOF", "PLOF")}
    # LOF vs DEL: LOF - DEL = pLoF-only carriers
    bad = shown["DEL"] & shown["LOF"] & small(w.LOF - w.DEL, k); shown["LOF"] &= ~bad
    # PLOF vs LOF: LOF - PLOF = deletion-only carriers; PLOF vs DEL (+ LOF): DEL + PLOF - LOF = carriers of both
    bad = shown["PLOF"] & ((shown["LOF"] & small(w.LOF - w.PLOF, k)) | (shown["LOF"] & shown["DEL"] & small(ov, k)))
    shown["PLOF"] &= ~bad
    for c in ("LOF", "PLOF"):
        genes = set(w.index[~shown[c] & (w[c] > 0)])
        m |= (d.dosage == c) & d.entrez.isin(genes)
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True); ap.add_argument("--export", required=True); ap.add_argument("--strata", default="all")
    ap.add_argument("--min_cell", type=int, default=20); ap.add_argument("--min_carriers", type=int, default=1)
    ap.add_argument("--sets", default="gene_assoc,gene_assoc_rare")
    ap.add_argument("--partitions", default='{"sex": ["female", "male"], "ancestry": ["eur", "afr", "amr", "eas", "sas", "mid"]}')
    a = ap.parse_args(); k = a.min_cell
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
        sec = {s: set() for s in T}                                      # secondary (cross-stratum): keys whose count is withheld in children
        if k > 0 and "all" in T:
            pc = T["all"].set_index(KEY).n_carriers
            for pname, kids in parts.items():
                kids = [c for c in kids if c in T]
                if not kids:
                    continue
                tot_k = pd.Series(0.0, index=pc.index); bad = pd.Series(False, index=pc.index)
                for c in kids:
                    v = T[c].set_index(KEY).n_carriers.reindex(pc.index).fillna(0); tot_k += v; bad |= small(v, k)
                bad |= small(pc - tot_k, k)
                for c in kids:
                    sec[c] |= set(bad[bad].index)
        for s, d in T.items():
            d = d[d.n_carriers >= a.min_carriers].reset_index(drop=True)
            prim = small(d.n_carriers, k) if k > 0 else pd.Series(False, index=d.index)
            sc = pd.Series([key in sec[s] for key in zip(d.entrez, d.dosage)], index=d.index) if sec[s] else pd.Series(False, index=d.index)
            cl = class_mask(d, k)
            mask = prim | sc | cl
            o = d.copy()
            o["n_carriers"] = o.n_carriers.astype("Int64"); o.loc[mask, "n_carriers"] = pd.NA
            o.loc[mask, ["freq", "beta", "SE"]] = np.nan
            # 1 = count 1..min_cell (frequency below (min_cell + 1) / N, the standard '<=20' label); 2 = count > min_cell but withheld because a
            # subtraction across strata or classes would give a small count
            o["n_carriers_masked"] = np.where(prim, 1, np.where(mask, 2, 0))
            o["freq_max"] = np.where(prim, (k + 1) / o.N, np.nan)
            if "max_pos_n" in o.columns:
                if k > 0:
                    o.loc[small(o.max_pos_n.fillna(0), k) | mask, "max_pos_freq"] = np.nan
                o = o.drop(columns="max_pos_n")
            if "n_overlap" in o.columns:
                o = o.drop(columns="n_overlap")
            dst = os.path.join(a.export, s, name, t); os.makedirs(dst)
            o.to_csv(os.path.join(dst, "gene_assoc.tsv.gz"), sep="\t", index=False, float_format="%.6g")
            info = json.load(open(os.path.join(a.results, s, name, t, "run_info.json")))
            info = {kk: info[kk] for kk in ("cohort", "stratum", "trait", "N", "covariates", "panels", "n_panels", "genes", "test", "segdup_bed",
                                            "cnv_table", "cnv_filter", "plof_table") if kk in info}
            info.update(min_cell=k, min_carriers=a.min_carriers, rows=int(len(o)), counts_withheld_primary=int(prim.sum()),
                        counts_withheld_secondary=int((sc & ~prim).sum()), counts_withheld_cross_class=int((cl & ~prim & ~sc).sum()),
                        rule="all rows exported with t and P; on rows with counts 1..min_cell or counts recoverable by subtraction, n_carriers, freq, beta and SE are withheld")
            json.dump(info, open(os.path.join(dst, "run_info.json"), "w"), indent=1)
            summary.append((name, s, t, len(o), int(mask.sum())))
    bad = []
    for root, _, files in os.walk(a.export):
        for f in sorted(files):
            p = os.path.join(root, f); sz = os.path.getsize(p); tot += sz; lines.append(f"{sz}\t{os.path.relpath(p, a.export)}")
            if f.endswith(".tsv.gz"):
                x = pd.read_csv(p, sep="\t")
                if ((k > 0 and small(x.n_carriers.fillna(0), k).any()) or "max_pos_n" in x.columns or "n_overlap" in x.columns
                        or (k > 0 and (x.n_carriers.isna() & (x.SE.notna() | x.beta.notna() | x.freq.notna())).any())):
                    bad.append(p)
            if any(x in f for x in ("pheno", "genotype", ".cnv", "plof_carriers")):
                bad.append(p)
    with open(os.path.join(a.export, "EXPORT_MANIFEST.txt"), "w") as fh:
        fh.write(f"gene-based association export; every row kept; counts 1..{k} and subtraction-recoverable counts withheld\n"
                 f"audit: {'OK' if not bad else 'FAILED'}\ntotal bytes {tot} ({tot / 1e6:.1f} MB) in {len(lines)} files\nbytes\tfile\n" + "\n".join(lines) + "\n")
    print("set\tstratum\ttrait\trows_exported\tcounts_withheld")
    for r in summary:
        print("\t".join(map(str, r)))
    print(f"exported to {a.export}: {tot / 1e6:.1f} MB")
    print("\n".join("AUDIT VIOLATION " + b for b in bad) or f"AUDIT OK: no shown count in 1..{k}; no person-level file")
    if bad:
        raise SystemExit("export audit failed: do not download this folder")


if __name__ == "__main__":
    main()
