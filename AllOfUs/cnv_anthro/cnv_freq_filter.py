"""Per-position carrier frequency of a CNV call set, and an OPTIONAL call filter on it (off by default).

Jon 2026-11-23, final decision: Fig. 3 uses the PGC freeze filters in every cohort (>= 10 kb; arrays >= 10 probes; calls carried by > 10% of
samples at >= 50% reciprocal overlap removed), exactly as the UK Biobank table was built: "Your 1% rule will eliminate some of the recurrent
CNVs like 15q11.2. It's better to adopt the pgc filters not the other way around." So by default this script only RECORDS the per-position
frequency ("Just record probe level frequency") and tabulates what a position filter would remove; with --max_freq it also writes a filtered
call table (kept for sensitivity analyses; not used by Fig. 3).

Earlier text of this docstring, the rule of the optional filter (Jon 2026-11-23: "The fig 3 dosage architecture analysis would only be comparable if it
uses rare CNVs ... we should do a filtering step for the purposes of fig. 3"; "we filter at the probe level").

Per dosage (DEL, DUP separately) and position x of the autosomes:
    f(x) = people in the dataset with a call of that dosage covering x  /  N        (N = every person of --pheno, carriers or not)
f is piecewise constant between call breakpoints, so it is exact at every base; at a breakpoint it equals the probe frequency NCNVs / N of the
CNV GWAS (cnv_gwas_qt.py, Omar Shanta's PGC pipeline, whose probe filter is frequency > 2%). A person's overlapping calls count once.
Rule: a call is DROPPED when at least --min_frac (0.5) of its length lies where f(x) > --max_freq, i.e. when most of the call is common.
A rare call that only clips a common region is kept; a common call is dropped whatever its exact breakpoints (unlike a reciprocal-overlap
rule, which fragmented WGS breakpoints can evade).
This comes on top of the call-set filters every cohort already has (>= 10 kb; arrays >= 10 probes; the PGC rule: calls carried by > 10% of
samples at >= 50% reciprocal overlap removed).
Outputs (person-level file stays in the cohort's storage):
  --out_cnv                 the filtered call table (same PLINK-CNV columns)
  --out_dir/position_freq.tsv.gz   CHR START END TYPE n_carriers freq, one row per stretch between breakpoints with >= 1 carrier
                                   (aggregate; the export keeps only rows above the small-cell limit)
  --out_dir/run_info.json   N, threshold, calls and people before / after, by dosage; the same counts at the --report thresholds
null model: none (a deterministic filter).
  python cnv_freq_filter.py --cnv X.cnv --pheno pheno.tsv --max_freq 0.01 --out_cnv X_rare.cnv --out_dir DIR [--report 0.005,0.01,0.02,0.05]
"""
import argparse
import json
import os
import numpy as np
import pandas as pd
from common import read_plink_cnv

TYPES = {1: "DEL", 3: "DUP"}


def coverage(c):
    """c: calls of one chromosome and dosage (IID, BP1, BP2). -> breakpoints p (sorted) and carriers n on [p[i], p[i+1])."""
    c = c.sort_values(["IID", "BP1"])
    prev = c.groupby("IID").BP2.cummax().groupby(c.IID).shift()
    grp = (prev.isna() | (c.BP1 > prev)).cumsum()                  # a person's overlapping calls merged: counted once
    m = c.groupby(grp.to_numpy()).agg(BP1=("BP1", "min"), BP2=("BP2", "max"))
    pos = np.concatenate([m.BP1.to_numpy(np.int64), m.BP2.to_numpy(np.int64) + 1]); d = np.r_[np.ones(len(m)), -np.ones(len(m))]
    o = np.argsort(pos, kind="stable"); pos, d = pos[o], d[o]
    p, first = np.unique(pos, return_index=True)
    return p, np.add.reduceat(d, first).cumsum().astype(np.int64)


def common_fraction(p, n, N, thr, b1, b2):
    """Fraction of each call [b1, b2] lying where n / N > thr."""
    bad = (n / N > thr).astype(float); seg = np.diff(p, append=p[-1]); cb = np.r_[0.0, np.cumsum(bad * seg)]
    def L(x):
        i = np.clip(np.searchsorted(p, x, "right") - 1, 0, len(p) - 1)
        return cb[i] + (x - p[i]) * bad[i]
    return (L(b2 + 1) - L(b1)) / (b2 + 1 - b1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cnv", required=True); ap.add_argument("--pheno", required=True, help="phenotype table: its IIDs are the dataset (N)")
    ap.add_argument("--max_freq", type=float, default=None, help="optional: also write --out_cnv without the calls that fail this frequency")
    ap.add_argument("--min_frac", type=float, default=0.5)
    ap.add_argument("--out_cnv", default=""); ap.add_argument("--out_dir", required=True)
    ap.add_argument("--report", default="0.005,0.01,0.02,0.05", help="other thresholds to tabulate in run_info.json")
    a = ap.parse_args(); os.makedirs(a.out_dir, exist_ok=True)
    ids = pd.read_csv(a.pheno, sep="\t", usecols=["IID"], dtype=str).IID.drop_duplicates(); N = len(ids)
    cnv = read_plink_cnv(a.cnv, keep_iids=ids); cnv = cnv[cnv.CHR <= 22].reset_index(drop=True)
    ths = sorted({*([a.max_freq] if a.max_freq is not None else []), *[float(x) for x in a.report.split(",") if x]})
    frac = {t: np.zeros(len(cnv)) for t in ths}; segs = []
    for (ch, ty), c in cnv.groupby(["CHR", "TYPE"]):
        if ty not in TYPES:
            continue
        p, n = coverage(c[["IID", "BP1", "BP2"]])
        for t in ths:
            frac[t][c.index.to_numpy()] = common_fraction(p, n, N, t, c.BP1.to_numpy(np.int64), c.BP2.to_numpy(np.int64))
        k = n[:-1] > 0
        segs.append(pd.DataFrame(dict(CHR=ch, START=p[:-1][k], END=p[1:][k] - 1, TYPE=TYPES[ty], n_carriers=n[:-1][k])))
    S = pd.concat(segs, ignore_index=True); S["freq"] = S.n_carriers / N
    S.to_csv(os.path.join(a.out_dir, "position_freq.tsv.gz"), sep="\t", index=False, float_format="%.6g")
    drop = (frac[a.max_freq] >= a.min_frac) if a.max_freq is not None else np.zeros(len(cnv), bool)
    out = cnv[~drop]; cols = [c for c in ["FID", "IID", "CHR", "BP1", "BP2", "TYPE", "SCORE", "SITES"] if c in out.columns]
    if a.max_freq is not None:
        if not a.out_cnv:
            raise SystemExit("--max_freq needs --out_cnv")
        out[cols].to_csv(a.out_cnv, sep="\t", index=False)
    tname = cnv.TYPE.map(TYPES)

    def tab(d):
        r = {}
        for ty in ("DEL", "DUP"):
            m = (tname == ty).to_numpy(); kept = cnv[m & ~d]
            r[ty] = dict(calls=int(m.sum()), calls_dropped=int((m & d).sum()), carriers=int(cnv[m].IID.nunique()), carriers_after=int(kept.IID.nunique()),
                         people_with_a_dropped_call=int(cnv[m & d].IID.nunique()), calls_per_person_after=round(float(len(kept) / N), 4),
                         kb_dropped_median=float(((cnv[m & d].BP2 - cnv[m & d].BP1 + 1) / 1000).median()) if (m & d).any() else None)
        return r
    info = dict(N=int(N), max_freq=a.max_freq, min_frac=a.min_frac,
                rule=("call dropped when >= min_frac of its length has same-dosage carrier frequency > max_freq" if a.max_freq is not None else
                      "no call dropped: per-position frequency recorded only (the call table already carries the PGC freeze filters)"),
                input=os.path.basename(a.cnv), calls_in=int(len(cnv)), calls_out=int(len(out)), by_dosage=tab(drop),
                max_position_freq={ty: float(S[S.TYPE == ty].freq.max()) for ty in ("DEL", "DUP")},
                mb_above_1pct={ty: round(float((S[(S.TYPE == ty) & (S.freq > 0.01)].eval("END - START + 1")).sum() / 1e6), 3) for ty in ("DEL", "DUP")},
                report={str(t): tab(frac[t] >= a.min_frac) for t in ths})
    json.dump(info, open(os.path.join(a.out_dir, "run_info.json"), "w"), indent=1)
    print(json.dumps({k: info[k] for k in ("N", "max_freq", "calls_in", "calls_out", "by_dosage", "max_position_freq", "mb_above_1pct")}, indent=1))
    for t in ths:
        print(f"a position filter at {t:.3%} would drop: " + "; ".join(f"{ty} calls dropped {v['calls_dropped']} of {v['calls']} ({v['people_with_a_dropped_call']} people)" for ty, v in info["report"][str(t)].items()))


if __name__ == "__main__":
    main()
