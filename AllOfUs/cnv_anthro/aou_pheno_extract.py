"""All of Us adapter (2/2): height, BMI, age, sex and ancestry PCs from the AoU OMOP CDR -> the standard phenotype
table (same columns prep_ukbb.py writes). Runs inside the Researcher Workbench only (BigQuery; WORKSPACE_CDR).

Choices (stated):
  * height = OMOP measurement concept 3036277 (LOINC 8302-2 'Body height'); BMI = 3038553 (LOINC 39156-5).
    Units normalised by unit concept name (inch -> cm); plausible ranges height 120-220 cm, BMI 12-75 kg/m2.
  * adults only (age >= 18 at measurement). Height and BMI = median of all adult measurements of that trait;
    age = median age of those measurements (age and age^2 of the trait's own measurements are covariates).
  * sex = sex at birth (person.sex_at_birth_concept_id, mapped by concept name Male/Female; others dropped).
  * ancestry PCs = first --n_pcs entries of pca_features in the AoU ancestry release (ancestry_preds.tsv).
  * restricted to people in the CNV source (--keep_samples = SV VCF sample list) so non-sequenced people are
    never counted as non-carriers; optionally to AoU's maximal unrelated set.
  * traits rank-inverse-normal transformed within sex (Blom), as for UKBB.
Individual-level output stays in the workspace bucket / persistent disk; nothing here is exported.
"""
import argparse
import ast
import json
import os
import numpy as np
import pandas as pd
from common import rint_within

SQL = """
SELECT m.person_id, m.measurement_concept_id AS cid, m.measurement_date, m.value_as_number AS v,
       u.concept_name AS unit, p.birth_datetime, sx.concept_name AS sex_name
FROM `{cdr}.measurement` m
JOIN `{cdr}.person` p ON p.person_id = m.person_id
LEFT JOIN `{cdr}.concept` u ON u.concept_id = m.unit_concept_id
LEFT JOIN `{cdr}.concept` sx ON sx.concept_id = p.sex_at_birth_concept_id
WHERE m.measurement_concept_id IN (3036277, 3038553) AND m.value_as_number IS NOT NULL
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("--summary", required=True)
    ap.add_argument("--ancestry", required=True, help="AoU ancestry_preds.tsv (research_id, ancestry_pred, pca_features)")
    ap.add_argument("--keep_samples", required=True, help="IDs present in the CNV source (SV VCF sample list)")
    ap.add_argument("--unrelated", default=None, help="optional AoU maximal-unrelated sample list")
    ap.add_argument("--n_pcs", type=int, default=16)
    ap.add_argument("--batch", default=None, help="array dataset only: TSV without header, IID <tab> genotyping batch -> column array_batch")
    ap.add_argument("--cdr", default=os.environ.get("WORKSPACE_CDR"))
    ap.add_argument("--measurements_csv", default=None, help="testing: pre-extracted rows instead of BigQuery")
    a = ap.parse_args()
    if a.measurements_csv:
        d = pd.read_csv(a.measurements_csv, parse_dates=["measurement_date", "birth_datetime"])
    else:
        from google.cloud import bigquery
        d = bigquery.Client().query(SQL.format(cdr=a.cdr)).to_dataframe()
    log = {"measurement_rows": int(len(d))}
    d["birth_datetime"] = pd.to_datetime(d["birth_datetime"], utc=True).dt.tz_localize(None)
    d["measurement_date"] = pd.to_datetime(d["measurement_date"])
    d["age"] = (d.measurement_date - d.birth_datetime).dt.days / 365.25
    d = d[d.age >= 18].copy(); d["v"] = d["v"].astype(float)
    u = d.unit.fillna("").str.lower()
    d.loc[(d.cid == 3036277) & u.str.contains("inch"), "v"] *= 2.54
    h = d[(d.cid == 3036277) & d.v.between(120, 220)].sort_values("measurement_date")
    b = d[(d.cid == 3038553) & d.v.between(12, 75)].sort_values("measurement_date")
    H = h.groupby("person_id").agg(height=("v", "median"), age_h=("age", "median"), sex_name=("sex_name", "first"))
    B = b.groupby("person_id").agg(bmi=("v", "median"), age_b=("age", "median"), sex_name_b=("sex_name", "first"))
    P = H.join(B, how="outer")
    P["sex_name"] = P.sex_name.fillna(P.sex_name_b)
    P["sex"] = P.sex_name.map({"Male": 1, "Female": 2})
    P = P[P.sex.notna()].reset_index().rename(columns={"person_id": "IID"})
    P["IID"] = P.IID.astype(str)
    keep = set(l.strip() for l in open(a.keep_samples))
    log["with_trait_and_sex"] = int(len(P)); P = P[P.IID.isin(keep)]; log["in_cnv_source"] = int(len(P))
    if a.unrelated:
        un = set(l.strip() for l in open(a.unrelated)); P = P[P.IID.isin(un)]; log["unrelated"] = int(len(P))
    an = pd.read_csv(a.ancestry, sep="\t", dtype={"research_id": str})
    pcs = np.array([ast.literal_eval(x)[: a.n_pcs] for x in an.pca_features])
    A = pd.DataFrame(pcs, columns=[f"PC{i + 1}" for i in range(a.n_pcs)]); A["IID"] = an.research_id; A["ancestry_pred"] = an.ancestry_pred
    P = P.merge(A, on="IID", how="inner"); log["with_PCs"] = int(len(P))
    if a.batch:   # array dataset (Jon 2026-11-17): genotyping batch as a factor covariate, as for UKBB
        B = pd.read_csv(a.batch, sep="\t", header=None, names=["IID", "array_batch"], dtype=str)
        P = P.merge(B, on="IID", how="inner"); log["with_array_batch"] = int(len(P))
    # one covariate age per person: the trait's own age where analysed (height and BMI rows can differ)
    P["age"] = P.age_h.fillna(P.age_b); P["age2"] = P.age ** 2
    out = {}
    for t, ag in [("height", "age_h"), ("bmi", "age_b")]:
        P[f"{t}_rint"] = rint_within(P[t], P.sex)
        out[f"N_{t}"] = int(P[t].notna().sum())
    P["age_bmi"] = P.age_b; P["age_bmi2"] = P.age_b ** 2
    log.update(out)
    log["N_by_ancestry_pred"] = {k: (v if v > 20 else "<=20") for k, v in P.ancestry_pred.value_counts().to_dict().items()}  # small-cell rule
    P.drop(columns=["sex_name", "sex_name_b"]).to_csv(a.out, sep="\t", index=False)
    json.dump(log, open(a.summary, "w"), indent=1)
    print(json.dumps(log, indent=1))


if __name__ == "__main__":
    main()
