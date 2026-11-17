"""Apply the UKBB / PGC Freeze_2024 call filters to an array CNV call set in PLINK-CNV format (AoU array dataset, Jon 2026-11-17).

Filters, identical to the UKBB calls in the PGC freeze (data_10perc, all_platforms_10kb_10probes):
  * autosomes 1-22, DEL (TYPE 1) and DUP (TYPE 3) only
  * length >= --min_size (10 kb) and probes (SITES) >= --min_probes (10)
  * frequency: a call is removed when > --max_freq (10%) of samples carry a same-type call at >= 50% reciprocal overlap
Samples: --samples lists every array sample that passed sample QC (people with no call count in N). --exclude_samples drops people
(e.g. the srWGS samples, so the array and WGS datasets do not share participants); the kept list is written to --write_samples for
aou_pheno_extract.py --keep_samples. Runs inside the Researcher Workbench; writes participant-level files to the VM only.
"""
import argparse
import sys
import pandas as pd
from aou_sv_to_cnv import pgc_frequency_drop


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cnv", required=True, help="PLINK-CNV: FID IID CHR BP1 BP2 TYPE SCORE SITES (GRCh38)")
    ap.add_argument("--samples", required=True, help="one IID per line: all array samples passing sample QC")
    ap.add_argument("--exclude_samples", default="", help="IIDs to drop (e.g. sv_samples.txt from aou_sv_to_cnv.py)")
    ap.add_argument("--out", required=True); ap.add_argument("--write_samples", required=True)
    ap.add_argument("--min_size", type=int, default=10000); ap.add_argument("--min_probes", type=int, default=10)
    ap.add_argument("--max_freq", type=float, default=0.10)
    a = ap.parse_args()
    keep = set(pd.read_csv(a.samples, header=None, dtype=str)[0])
    n_all = len(keep)
    if a.exclude_samples:
        keep -= set(pd.read_csv(a.exclude_samples, header=None, dtype=str)[0])
    c = pd.read_csv(a.cnv, sep=r"\s+", dtype={"FID": str, "IID": str})
    c["CHR"] = pd.to_numeric(c.CHR.astype(str).str.replace("chr", "", regex=False), errors="coerce")
    n0 = len(c)
    c = c[c.IID.isin(keep) & c.CHR.between(1, 22) & c.TYPE.isin([1, 3])]
    n1 = len(c)
    c = c[((c.BP2 - c.BP1) >= a.min_size) & (c.SITES >= a.min_probes)].reset_index(drop=True)
    n2 = len(c)
    c["CHR"] = c.CHR.astype(int)
    site = c.groupby(["CHR", "BP1", "BP2", "TYPE"]).IID.agg(lambda s: set(s)).reset_index()
    recs = [(r.CHR, r.BP1, r.BP2, r.TYPE, r.IID) for r in site.itertuples()]
    drop = pgc_frequency_drop(recs, len(keep), a.max_freq)
    site["drop"] = [j in drop for j in range(len(site))]
    c = c.merge(site[["CHR", "BP1", "BP2", "TYPE", "drop"]], on=["CHR", "BP1", "BP2", "TYPE"])
    c = c[~c["drop"]].drop(columns="drop")
    c.to_csv(a.out, sep="\t", index=False)
    pd.Series(sorted(keep)).to_csv(a.write_samples, index=False, header=False)
    per = c.groupby("IID").size().reindex(sorted(keep), fill_value=0)
    print(f"samples {n_all} -> {len(keep)} after exclusion; calls {n0} -> {n1} (samples, autosomes, DEL/DUP) -> {n2} "
          f"(>= {a.min_size} bp, >= {a.min_probes} probes) -> {len(c)} (frequency <= {a.max_freq:.0%}); "
          f"calls/person median {per.median():.0f}, people with no call {(per == 0).mean():.1%}", file=sys.stderr)


if __name__ == "__main__":
    main()
