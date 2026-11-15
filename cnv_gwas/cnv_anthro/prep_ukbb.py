"""UK Biobank adapter: build the standard phenotype/covariate table used by every analysis module.

Standard table (one row per analysed person, tab-separated):
  IID  sex(1=male,2=female)  age  age2  batch  array  PC1..PC10  height  bmi  height_rint  bmi_rint  [+ cohort extras]
Trait transform (stated choice): rank-based inverse normal transform WITHIN SEX (Blom), applied to the analysed sample.

The All of Us adapter only has to write the same columns; nothing downstream reads UKBB fields.

UKBB inputs:
  --pgc_phe   PGC Freeze_2024 XD .phe (UKBB rows = dataset 'ukbb'; IID = genotype-project IID; C1..C10 = PCs)
  --link      Molly Sacks's genomatrix_UKBB_2025.csv (SNP_Genotype_IID -> f.eid; duplicate_drop; withdrawn_consent)
  --fields    eid + fields 31, 34, 50, 54, 21001, 21003, 21022, 22000, 22001, 22006, 22019, 22020, 22021, 22027 (instance 0)
"""
import argparse
import json
import numpy as np
import pandas as pd
from common import rint_within


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pgc_phe", required=True)
    ap.add_argument("--link", required=True)
    ap.add_argument("--fields", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--pgc_cnv", default=None, help="PGC XD .cnv; UKBB calls are written to --out_cnv")
    ap.add_argument("--out_cnv", default=None)
    a = ap.parse_args()
    log = {}

    phe = pd.read_csv(a.pgc_phe, sep="\t", dtype={"FID": str, "IID": str}, keep_default_na=False)
    phe = phe[phe["dataset"] == "ukbb"].copy()
    log["pgc_ukbb_rows"] = len(phe)
    log["pgc_ukbb_unique_IID"] = phe["IID"].nunique()
    phe = phe.drop_duplicates("IID")
    log["pgc_ukbb_by_AFF"] = phe["AFF"].astype(str).value_counts().to_dict()

    lk = pd.read_csv(a.link, usecols=["SNP_Genotype_IID", "f.eid", "duplicate_drop", "withdrawn_consent", "sex_combined"],
                     dtype=str, keep_default_na=False)
    lk["IID"] = lk["SNP_Genotype_IID"].str.replace(r"\.0$", "", regex=True)
    lk = lk[lk["IID"] != ""]
    d = phe.merge(lk, on="IID", how="inner")
    log["linked_to_eid"] = len(d)
    d = d[(d["duplicate_drop"] != "duplicate_drop") & (d["withdrawn_consent"] != "yes")]
    log["after_duplicate_withdrawn"] = len(d)

    f = pd.read_csv(a.fields, sep="\t", dtype=str, keep_default_na=False, na_values=["NA", ""])
    f = f.rename(columns={"f.eid": "eid", "f.31.0.0": "sex_rep", "f.50.0.0": "height", "f.54.0.0": "centre",
                          "f.21001.0.0": "bmi", "f.21003.0.0": "age", "f.22000.0.0": "batch", "f.22001.0.0": "sex_gen",
                          "f.22019.0.0": "sex_aneuploidy", "f.22027.0.0": "het_missing_outlier", "f.22006.0.0": "white_british",
                          "f.22021.0.0": "kinship"})
    d = d.rename(columns={"f.eid": "eid"}).merge(f, on="eid", how="inner")
    log["with_ukb_fields"] = len(d)
    for c in ["height", "bmi", "age", "batch", "sex_rep", "sex_gen", "sex_aneuploidy", "het_missing_outlier"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")

    # QC exclusions (counted in order)
    ex = {}
    m = d["sex_aneuploidy"] == 1; ex["sex_chromosome_aneuploidy_22019"] = int(m.sum()); d = d[~m]
    m = d["het_missing_outlier"] == 1; ex["het_missing_outlier_22027"] = int(m.sum()); d = d[~m]
    m = d["sex_gen"].notna() & d["sex_rep"].notna() & (d["sex_gen"] != d["sex_rep"]); ex["reported_vs_genetic_sex_mismatch"] = int(m.sum()); d = d[~m]
    m = d["age"].isna() | d["batch"].isna() | d["sex_rep"].isna(); ex["missing_age_batch_or_sex"] = int(m.sum()); d = d[~m]
    m = d[[f"C{i}" for i in range(1, 11)]].apply(pd.to_numeric, errors="coerce").isna().any(axis=1); ex["missing_PCs"] = int(m.sum()); d = d[~m]
    log["exclusions"] = ex
    log["analysis_pool"] = len(d)

    out = pd.DataFrame({"IID": d["IID"], "eid": d["eid"],
                        "sex": np.where(d["sex_rep"] == 1, 1, 2), "age": d["age"], "age2": d["age"] ** 2,
                        "batch": d["batch"].astype(int).astype(str),
                        "array": np.where(d["batch"] < 0, "UKBiLEVE", "UKBAxiom"),
                        "height": d["height"], "bmi": d["bmi"], "AFF": d["AFF"], "ancestry_pgc": d["ancestry"]})
    for i in range(1, 11):
        out[f"PC{i}"] = pd.to_numeric(d[f"C{i}"]).to_numpy()
    for t in ["height", "bmi"]:
        out[f"{t}_rint"] = rint_within(out[t], out["sex"])   # NaN where the trait is missing
        log[f"N_{t}"] = int(out[t].notna().sum())
        log[f"N_{t}_by_sex"] = out.loc[out[t].notna(), "sex"].map({1: "male", 2: "female"}).value_counts().to_dict()
    log["N_by_array"] = out["array"].value_counts().to_dict()
    log["N_batches"] = int(out["batch"].nunique())
    log["N_by_pgc_ancestry"] = out["ancestry_pgc"].value_counts().to_dict()
    log["N_pgc_psychiatric_cases_AFF2"] = int((out["AFF"].astype(str) == "2").sum())
    out.to_csv(a.out, sep="\t", index=False)
    if a.pgc_cnv:
        c = pd.read_csv(a.pgc_cnv, sep=r"\s+", dtype={"FID": str, "IID": str}, keep_default_na=False)
        c = c[c.IID.isin(set(phe.IID))]
        c.to_csv(a.out_cnv, sep="\t", index=False)
        log["ukbb_calls_written"] = int(len(c))
        log["ukbb_calls_by_type"] = c.TYPE.astype(str).value_counts().to_dict()
    json.dump(log, open(a.summary, "w"), indent=1)
    print(json.dumps(log, indent=1))


if __name__ == "__main__":
    main()
