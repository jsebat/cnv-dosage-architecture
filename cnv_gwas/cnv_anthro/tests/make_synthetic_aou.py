"""Synthetic GATK-SV VCF + OMOP rows + ancestry file to unit-test the AoU adapters (no real data)."""
import gzip, json, os, numpy as np, pandas as pd
d = os.path.join(os.path.dirname(__file__), "synthetic"); os.makedirs(d, exist_ok=True)
S = [f"S{i}" for i in range(1, 7)]
rows = [  # chrom pos id filter info format genotypes
 ("chr1", 1000000, "PASS", "SVTYPE=DEL;END=1050000", "GT", ["0/1", "0/0", "1/1", "0/0", "./.", "0/0"]),   # DEL 50 kb: S1,S3
 ("chr1", 2000000, "PASS", "SVTYPE=DUP;END=2005000", "GT", ["0/1"] * 6),                                    # 5 kb: dropped (<10 kb)
 ("chr2", 3000000, "MULTIALLELIC", "SVTYPE=CNV;END=3100000", "GT:CN", ["./.:2", "./.:1", "./.:3", "./.:4", "./.:2", "./.:0"]),  # S2,S6 DEL; S3,S4 DUP
 ("chr3", 4000000, "LOW_QUAL", "SVTYPE=DEL;END=4100000", "GT", ["0/1"] * 6),                                 # non-PASS: dropped
 ("chr4", 5000000, "PASS", "SVTYPE=INV;END=5200000", "GT", ["0/1"] * 6),                                     # INV: ignored
 ("chrX", 6000000, "PASS", "SVTYPE=DEL;END=6100000", "GT", ["0/1"] * 6),                                     # chrX: excluded (1-22)
 ("chr5", 7000000, "PASS", "SVTYPE=DUP;END=7020000", "GT", ["0|1", "0/0", "0/0", "0/0", "0/0", "0/0"]),      # DUP: S1
]
with gzip.open(os.path.join(d, "sv.vcf.gz"), "wt") as f:
    f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(S) + "\n")
    for i, (c, p, fl, inf, fm, g) in enumerate(rows):
        f.write(f"{c}\t{p}\tsv{i}\tN\t<X>\t.\t{fl}\t{inf}\t{fm}\t" + "\t".join(g) + "\n")
m = []
for i, s in enumerate(S):
    m.append(dict(person_id=s, cid=3036277, measurement_date="2020-01-01", v=160 + i * 3, unit="centimeter", birth_datetime="1970-01-01", sex_name="Female" if i % 2 else "Male"))
    m.append(dict(person_id=s, cid=3038553, measurement_date="2019-01-01", v=22 + i, unit="kilogram per square meter", birth_datetime="1970-01-01", sex_name="Female" if i % 2 else "Male"))
m.append(dict(person_id="S1", cid=3036277, measurement_date="2021-01-01", v=67, unit="inch (US)", birth_datetime="1970-01-01", sex_name="Male"))
pd.DataFrame(m).to_csv(os.path.join(d, "measurements.csv"), index=False)
an = pd.DataFrame({"research_id": S, "ancestry_pred": ["eur"] * 6, "pca_features": [json.dumps(list(np.round(np.random.default_rng(i).normal(size=16), 3))) for i in range(6)]})
an.to_csv(os.path.join(d, "ancestry_preds.tsv"), sep="\t", index=False)
