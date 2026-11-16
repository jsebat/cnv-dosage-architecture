"""Allele-level carrier frequency at the Table S9 loci across cohorts / technologies (e.g. PGC-array UKBB vs AoU srWGS).
Reads the EXPORTED (small-cell-suppressed) allele_assoc.tsv of each cohort, so it can run outside the Workbench.
Output: one row per allele with frequency per cohort and the WGS/array ratio (direct read-out of the sensitivity
difference between array calls and GATK-SV calls for each allele)."""
import argparse
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--tables", required=True, help="comma list label=path/to/export/taskA_alleles/<trait>/allele_assoc.tsv")
ap.add_argument("--ref", required=True, help="label of the reference cohort for the ratio (e.g. UKBB_array)")
ap.add_argument("--out", required=True)
a = ap.parse_args()
M = None
for item in a.tables.split(","):
    lab, path = item.split("=", 1)
    d = pd.read_csv(path, sep="\t", dtype=str)
    d = d.drop_duplicates("allele")[["allele", "n_carriers", "N", "FREQ"]].rename(
        columns={"n_carriers": f"n_{lab}", "N": f"N_{lab}", "FREQ": f"freq_{lab}"})
    M = d if M is None else M.merge(d, on="allele", how="outer")
for c in [c for c in M.columns if c.startswith("freq_")]:
    M[c] = pd.to_numeric(M[c], errors="coerce")                # 'suppressed' (<= 20 carriers) -> NaN
for c in [c for c in M.columns if c.startswith("freq_") and c != f"freq_{a.ref}"]:
    M[f"ratio_{c[5:]}_vs_{a.ref}"] = M[c] / M[f"freq_{a.ref}"]
M["locus"] = M.allele.str.replace("^X_", "", regex=True)
M.to_csv(a.out, sep="\t", index=False, float_format="%.4g")
print(M.to_string(index=False))
