# Synthetic cohort (no real data) to test the person-level ConGA port end to end: phenotype table and PLINK .cnv in the AoU pipeline's
# formats, with random deletions/duplications on real gene coordinates and carriers of a few recurrent intervals (loci_sacks34) so
# that the locus-removal runs remove calls. Height is shifted in carriers of duplications or deletions of one context's genes.
# null model: n/a (test data)
"""python make_synthetic_ctx.py --exons <refGene exons> --loci34 <loci_sacks34.tsv> --sets <membership.tsv.gz> --signal_set <set_id>
       --n 4000 --seed 1 --out <dir>
Writes <out>/pheno.tsv (IID, sex, age, age2, age_bmi, age_bmi2, PC1-16, height, bmi, height_rint, bmi_rint) and <out>/calls.cnv."""
import argparse
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import rint_within  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exons", required=True); ap.add_argument("--loci34", required=True); ap.add_argument("--sets", required=True)
    ap.add_argument("--signal_set", required=True); ap.add_argument("--n", type=int, default=4000); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(); rng = np.random.default_rng(a.seed); os.makedirs(a.out, exist_ok=True)
    N = a.n; iid = [f"SYN{a.seed}_{i:05d}" for i in range(N)]
    sex = rng.integers(1, 3, N); age = rng.uniform(40, 70, N); pcs = rng.normal(size=(N, 16))
    ex = pd.read_csv(a.exons, sep="\t", header=None, names=["chr", "s", "e", "tx", "sym", "entrez"])
    ex = ex[ex.chr.isin([f"chr{c}" for c in range(1, 23)])]
    genes = ex.groupby("entrez").agg(chr=("chr", "first"), s=("s", "min"), e=("e", "max")).reset_index()
    genes = genes[(genes.e - genes.s) < 1_500_000]
    sm = pd.read_csv(a.sets, sep="\t", dtype={"entrez": str})
    sig = set(sm[sm.set_id == a.signal_set].entrez)
    calls = []; dose = np.zeros(N)
    for i in range(N):
        for _ in range(rng.poisson(1.2)):                     # random CNVs over a random gene, 20 kb-1 Mb flank
            g = genes.iloc[rng.integers(len(genes))]; typ = int(rng.choice([1, 3]))
            s = max(1, int(g.s) - int(rng.integers(20_000, 1_000_000))); e = int(g.e) + int(rng.integers(20_000, 1_000_000))
            calls.append((iid[i], int(g.chr[3:]), s, e, typ))
            if str(g.entrez) in sig:
                dose[i] += 1                                  # dose-aligned effect: both CNV types lower height
    L = pd.read_csv(a.loci34, sep="\t")                       # columns id, name, dosage, chr, start, end
    L = L[L.chr.astype(str).str.replace("chr", "").isin([str(c) for c in range(1, 23)])].drop_duplicates("id").iloc[::5].head(6)
    for r in L.itertuples():                                  # recurrent-interval carriers (exact interval), both CNV types
        for i in rng.choice(N, 25, replace=False):
            calls.append((iid[i], int(str(r.chr).replace("chr", "")), int(r.start), int(r.end), int(rng.choice([1, 3]))))
    height = 170 + 7 * (sex == 1) - 3.0 * dose + rng.normal(0, 7, N); bmi = 27 + rng.normal(0, 4, N)
    P = pd.DataFrame(dict(IID=iid, sex=sex, age=age, age2=age ** 2, age_bmi=age, age_bmi2=age ** 2, height=height, bmi=bmi))
    for k in range(16):
        P[f"PC{k + 1}"] = pcs[:, k]
    P["height_rint"] = rint_within(height, sex); P["bmi_rint"] = rint_within(bmi, sex)
    P.to_csv(os.path.join(a.out, "pheno.tsv"), sep="\t", index=False)
    C = pd.DataFrame(calls, columns=["IID", "CHR", "BP1", "BP2", "TYPE"]); C.insert(0, "FID", C.IID); C["SCORE"] = 0; C["SITES"] = 10
    C.to_csv(os.path.join(a.out, "calls.cnv"), sep="\t", index=False)
    print(f"N={N} calls={len(C)} signal-set dose carriers={(dose > 0).sum()}")


if __name__ == "__main__":
    main()
