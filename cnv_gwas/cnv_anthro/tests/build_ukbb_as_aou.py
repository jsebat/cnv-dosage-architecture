"""Plumbing test of the All of Us code path on UKBB data (Expanse only; outputs stay in the group-only data dir).
Writes AoU-shaped inputs for a random subset of UKBB people:
  * a GATK-SV-style VCF (one biallelic DEL/DUP site per distinct call interval, GT 0/1 for carriers);
  * OMOP-style measurement rows (height/BMI concept ids, dates and birth dates chosen to reproduce UKBB age);
  * an ancestry_preds.tsv with the UKBB PCs as pca_features.
Then configs/ukbb_as_aou_test.json runs the unchanged AoU adapters + all analysis steps on them.
"""
import argparse
import gzip
import json
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--pheno", required=True); ap.add_argument("--cnv", required=True); ap.add_argument("--outdir", required=True)
ap.add_argument("--n", type=int, default=5000); ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()
ph = pd.read_csv(a.pheno, sep="\t", dtype={"IID": str})
ph = ph.sample(a.n, random_state=a.seed).reset_index(drop=True)
ph["RID"] = [f"T{i:06d}" for i in range(len(ph))]                   # fresh test IDs, not UKBB IDs
c = pd.read_csv(a.cnv, sep=r"\s+", dtype={"IID": str}, keep_default_na=False)
c = c[c.IID.isin(set(ph.IID)) & (c.CHR <= 22)].merge(ph[["IID", "RID"]], on="IID")
sites = c.groupby(["CHR", "BP1", "BP2", "TYPE"]).RID.apply(set).reset_index().sort_values(["CHR", "BP1"])
S = ph.RID.tolist()
with gzip.open(f"{a.outdir}/sv_test.vcf.gz", "wt") as f:
    f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(S) + "\n")
    for i, r in enumerate(sites.itertuples()):
        t = "DEL" if r.TYPE == 1 else "DUP"
        g = ["0/1" if s in r.RID else "0/0" for s in S]
        f.write(f"chr{r.CHR}\t{r.BP1}\tt{i}\tN\t<{t}>\t.\tPASS\tSVTYPE={t};END={r.BP2}\tGT\t" + "\t".join(g) + "\n")
m = []
for r in ph.itertuples():
    birth = pd.Timestamp("2008-01-01") - pd.Timedelta(days=int(r.age * 365.25 + 100))
    sx = "Male" if r.sex == 1 else "Female"
    if np.isfinite(r.height):
        m.append((r.RID, 3036277, "2008-01-01", r.height, "centimeter", birth, sx))
    if np.isfinite(r.bmi):
        m.append((r.RID, 3038553, "2008-01-01", r.bmi, "kilogram per square meter", birth, sx))
pd.DataFrame(m, columns=["person_id", "cid", "measurement_date", "v", "unit", "birth_datetime", "sex_name"]).to_csv(f"{a.outdir}/measurements_test.csv", index=False)
pcs = ph[[f"PC{i}" for i in range(1, 11)]].to_numpy()
pd.DataFrame({"research_id": ph.RID, "ancestry_pred": "test",
              "pca_features": [json.dumps(list(np.r_[p, np.zeros(6)])) for p in pcs]}).to_csv(f"{a.outdir}/ancestry_test.tsv", sep="\t", index=False)
print(f"people {len(ph)}, sites {len(sites)}, calls {len(c)}")
