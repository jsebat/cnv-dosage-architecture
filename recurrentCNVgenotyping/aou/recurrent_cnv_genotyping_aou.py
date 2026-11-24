#!/usr/bin/env python3
"""Recurrent CNV genotyping for All of Us -- a port of Molly's recurrentCNVgenotyping/GenotypeRecurrentCNVs_tutorial.ipynb.

Same locus definitions (her public sheet, saved as recurrentCNV_hg38.reciprocal.tsv / .oneway.tsv), same rules:
  * reciprocal loci : a call is assigned if it overlaps the locus by >= 50% of the locus AND >= 50% of the call
                      (bedtools intersect -f .5 -r)
  * one-way loci    : a call is assigned if it covers >= 50% of the locus segment (bedtools intersect -f .5);
                      segments are then combined into allele-specific genotypes by the exact patterns of the
                      notebook (1q21.1 TAR, 15q11-13, 16p11.2, 16p12.2, 16p13.11, 22q11.21 A-D and D-H).
  * genotype codes  : 1 = DEL, 2 = no call, 3 = DUP; when a person has several calls at one locus the first one wins
                      (pivot_table aggfunc='first', calls in input order), as in the notebook.
Fragmented calls are kept as they are (no merging), as in the notebook: e.g. a 22q11.21 pattern '323' gets no A-D genotype.

Differences from the notebook (mechanics only, no rule changes):
  * no bedtools / Google Sheets: the overlap is computed in pandas with the same half-open arithmetic
    (overlap = min(end) - max(start)), and the two locus tabs are read from local TSVs (the Workbench has no
    Google access);
  * input is the pipeline's PLINK-CNV table (FID IID CHR BP1 BP2 TYPE ...; TYPE 1 = DEL, 3 = DUP), so the same
    script runs on AoU srWGS calls (aou_sv_to_cnv.py output) and AoU array calls (filter_plink_cnv.py output);
  * the genotype matrix is restricted to --samples (people analysed); no sample-QC NA step (calls are already QC'd).
Cross-check (crosscheck_molly.py: the notebook run cell by cell with its bedtools v2.28.0 binary, 4,000 simulated people,
~10,700 carrier genotypes): identical except (a) people with both a DEL and a DUP at one locus, where "first call" depends
on bedtools' output order, and (b) the notebook reads each bedtools output with pd.read_csv's default header=0, so the
first hit of each intersect is lost (one person per run); the port keeps it.

Outputs: --out genotype matrix (participant level: keep inside the Workbench) and, with --counts, the per-locus
DEL/DUP carrier counts with 1..min_cell masked (aggregate; safe to export).
"""
import argparse
import os
import re
import numpy as np
import pandas as pd

# ---- combination patterns, copied from the notebook (segment order matters) ----
PATTERNS = {
    # 15q11-13 BP1-BP3: segments 15q11.2_BP1_BP2, 15q11.2-13.1_BP2_BP3
    "15q13": (["15q11.2_BP1_BP2", "15q11.2-13.1_BP2_BP3"], {
        "15q11.2_BP1_BP2": ("12", "32"), "15q11.2-13.1_BP2_BP3": ("21", "23"), "15q11.2-13.1_BP1_BP3": ("11", "33")}),
    # 15q13 BP3-BP5: segments BP3_BP4, BP4_BP5_exclCHRNA7, CHRNA7
    "15q35": (["15q13.1-13.2_BP3_BP4", "15q13.3_BP4_BP5_exclCHRNA7", "15q13.3_CHRNA7"], {
        "15q13.1-13.2_BP3_BP4": ("122", "322"), "15q13.1-13.3_BP3_BP5": ("111", "333"), "15q13.3_BP4_BP5": ("211", "233"),
        "15q13.3_BP4_BP5_exclCHRNA7": ("212", "232"), "15q13.3_CHRNA7": ("221", "223")}),
    "16p11.2": (["16p11.2_BP1_BP2", "16p11.2_BP2_BP3", "16p11.2_BP3_BP4", "16p11.2_BP4_BP5"], {
        "16p11.2_BP1_BP2": ("1222", "3222"), "16p11.2_BP1_BP3": ("1122", "3322"), "16p11.2_BP2_BP3": ("2122", "2322"),
        "16p11.2_BP3_BP4": ("2212", "2232"), "16p11.2_BP4_BP5": ("2221", "2223")}),
    "22q11.21_AD": (["22q11.21_A_B", "22q11.21_B_C", "22q11.21_C_D"], {
        "22q11.21_A_B": ("122", "322"), "22q11.21_B_C": ("212", "232"), "22q11.21_C_D": ("221", "223"),
        "22q11.21_A_C": ("112", "332"), "22q11.21_B_D": ("211", "233"), "22q11.21_A_D": ("111", "333")}),
    "22q11.21_DH": (["22q11.21_D_E", "22q11.21_E_F", "22q11.21_F_G", "22q11.21_G_H"], {
        "22q11.21_D_E": ("1222", "3222"), "22q11.21_E_F": ("2122", "2322"), "22q11.21_F_G": ("2212", "2232"),
        "22q11.21_G_H": ("2221", "2223"), "22q11.21_F_H": ("2211", "2233"), "22q11.21_E_H": ("2111", "2333"),
        "22q11.21_D_H": ("1111", "3333")}),
}
# three-breakpoint regions genotyped with np.where in the notebook (1q21.1 TAR, 16p12.2, 16p13.11)
TRIPLES = [("1q21.1_TAR", "BP1_BP2", "BP2_BP3"), ("16p12.2", "BP1_BP2", "BP2_BP3"), ("16p13.11", "BP1_BP2", "BP2_BP3")]


def read_loci(path):
    d = pd.read_csv(path, sep="\t")
    d = d.rename(columns={"Locus name": "NAME"})[["CHR", "START", "END", "NAME"]]
    d["CHR"] = pd.to_numeric(d.CHR.astype(str).str.replace("chr", "", regex=False), errors="coerce")
    return d.dropna(subset=["CHR"]).astype({"CHR": int, "START": int, "END": int})


def assign(calls, loci, reciprocal):
    """genotype per (IID, locus) for the calls that pass the overlap rule; first call (input order) wins."""
    out = []
    for L in loci.itertuples():
        c = calls[calls.CHR == L.CHR]
        o = np.maximum(0, np.minimum(c.BP2.to_numpy(), L.END) - np.maximum(c.BP1.to_numpy(), L.START))
        m = (o > 0) & (o >= 0.5 * (L.END - L.START))
        if reciprocal:
            m &= o >= 0.5 * (c.BP2 - c.BP1).to_numpy()
        h = c[m]
        out.append(pd.DataFrame({"IID": h.IID.to_numpy(), "locus": L.NAME, "CN": h.TYPE.to_numpy(), "ord": h.ord.to_numpy()}))
    a = pd.concat(out, ignore_index=True).sort_values("ord", kind="mergesort")
    return a.drop_duplicates(["IID", "locus"], keep="first")


def genotype(cnv, samples, rec, one):
    cnv = cnv[cnv.TYPE.isin([1, 3]) & cnv.IID.isin(set(samples))].copy()
    cnv["ord"] = np.arange(len(cnv))
    G = pd.DataFrame(index=pd.Index(samples, name="IID"))
    for loci, recip, suffix in [(rec, True, ""), (one, False, "_oneway")]:
        a = assign(cnv, loci, recip)
        for name in loci.NAME:
            G[name + suffix] = 2
        if len(a):
            w = a.pivot_table(index="IID", columns="locus", values="CN", aggfunc="first")
            for name in w.columns:
                v = w[name].dropna()
                G.loc[v.index, name + suffix] = v.astype(int).to_numpy()
    # allele-specific genotypes from the one-way segments
    for (reg, a1, a2) in TRIPLES:
        s1, s2 = G[f"{reg}_{a1}_oneway"], G[f"{reg}_{a2}_oneway"]
        G[f"{reg}_{a1}"] = np.where(s2 == 2, s1, 2)
        G[f"{reg}_{a2}"] = np.where(s1 == 2, s2, 2)
        G[f"{reg}_BP1_BP3"] = np.where(s1 == s2, s2, 2)
    for key, (segs, alleles) in PATTERNS.items():
        pat = G[[s + "_oneway" for s in segs]].astype(int).astype(str).agg("".join, axis=1)
        for allele, (pdel, pdup) in alleles.items():
            G[allele] = np.select([pat == pdel, pat == pdup], [1, 3], default=2)
    return G


def counts(G, min_cell):
    rows = []
    for col in [c for c in G.columns if not c.endswith("_oneway")]:
        nd, nu = int((G[col] == 1).sum()), int((G[col] == 3).sum())
        m = lambda n: f"<={min_cell}" if 1 <= n <= min_cell else n
        rows.append(dict(locus=col, n_DEL=m(nd), n_DUP=m(nu), N=len(G)))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------------------------------
# Pipeline API (Jon 2026-11-24: "Milad's claude is ignoring our recurrent cnv genotyping code and constantly reinventing
# the wheel"). Every step of AllOfUs/cnv_anthro and cnv_gwas/cnv_anthro that needs recurrent-CNV carriers (Fig. 3 / grid /
# joint locus-removal runs, growth-locus effects, asymmetric-interval strata, additivity exclusions) calls these functions.
# There is no second genotyper: do not add overlap rules for recurrent loci anywhere else.
# null model: n/a (genotyping only).
# ---------------------------------------------------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
SHEET_RECIPROCAL = os.path.join(HERE, "recurrentCNV_hg38.reciprocal.tsv")
SHEET_ONEWAY = os.path.join(HERE, "recurrentCNV_hg38.oneway.tsv")
CODE = {"DEL": 1, "DUP": 3}


def locus_id(name):
    """Locus name in any of the repo's spellings -> one id ('22q11.21_A-D', '22q11.21_A_D' -> '22q11_21_A_D')."""
    return re.sub(r"[.\-]", "_", str(name))


def genotype_calls(cnv):
    """Genotype matrix for every person with a call: IID x allele (columns = locus_id of Molly's allele names),
    codes 1 = DEL, 2 = none, 3 = DUP. cnv: PLINK-CNV DataFrame (IID CHR BP1 BP2 TYPE; TYPE 1/3), in file order."""
    c = cnv[["IID", "CHR", "BP1", "BP2", "TYPE"]].copy()
    c["CHR"] = pd.to_numeric(c.CHR.astype(str).str.replace("chr", "", regex=False), errors="coerce")
    c = c.dropna(subset=["CHR"]).astype({"CHR": int, "BP1": np.int64, "BP2": np.int64, "TYPE": int})
    G = genotype(c, list(pd.unique(c.IID)), read_loci(SHEET_RECIPROCAL), read_loci(SHEET_ONEWAY))
    G = G[[x for x in G.columns if not x.endswith("_oneway")]]
    G.columns = [locus_id(x) for x in G.columns]
    if G.columns.duplicated().any():
        raise ValueError(f"duplicate allele ids: {list(G.columns[G.columns.duplicated()])}")
    return G


def _col(G, name):
    k = locus_id(name)
    if k not in G.columns:
        raise KeyError(f"{name}: not an allele of Molly's recurrent CNV sheet ({len(G.columns)} alleles: {', '.join(G.columns)})")
    return G[k]


def carriers(G, names, dosage="both"):
    """Set of IIDs carrying any of the named alleles; dosage 'DEL', 'DUP' or 'both'."""
    codes = [1, 3] if dosage == "both" else [CODE[dosage]]
    out = set()
    for n in ([names] if isinstance(names, str) else names):
        g = _col(G, n)
        out |= set(g.index[g.isin(codes)])
    return out


def label_calls(cnv, loci, G=None):
    """Per-call locus label (array of loci.id, '' = none). loci: table with id, dosage, chr, start, end[, match]
    (resources/growth_loci/loci_sacks34*.tsv). A call gets recurrent locus L when its person carries L's allele with the
    call's dosage in Molly's genotypes and the call overlaps L (so every piece of a fragmented call is labelled).
    Rows with match == 'gene' (UKBB CNV-GWAS loci, not recurrent CNVs) keep their rule: any overlap of the same dosage.
    The first matching row of the table wins."""
    if G is None:
        G = genotype_calls(cnv)
    lab = np.array([""] * len(cnv), dtype=object)
    iid = cnv.IID.to_numpy(); ch = cnv.CHR.to_numpy(); a = cnv.BP1.to_numpy(); b = cnv.BP2.to_numpy()
    dos = np.where(cnv.TYPE.to_numpy() < 2, "DEL", "DUP")
    for L in loci.itertuples():
        m = (dos == L.dosage) & (ch == int(L.chr)) & (a <= L.end) & (b >= L.start) & (lab == "")
        if getattr(L, "match", "segment") != "gene":
            m &= np.isin(iid, list(carriers(G, L.id, L.dosage)))
        lab[m] = L.id
    return lab


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cnv", required=True, help="PLINK-CNV table: FID IID CHR BP1 BP2 TYPE ... (TYPE 1 = DEL, 3 = DUP)")
    ap.add_argument("--samples", required=True, help="one IID per line, or a TSV with an IID column (people analysed)")
    ap.add_argument("--reciprocal", required=True, help="recurrentCNV_hg38.reciprocal.tsv (Molly's sheet, reciprocal tab)")
    ap.add_argument("--oneway", required=True, help="recurrentCNV_hg38.oneway.tsv (Molly's sheet, one-way tab)")
    ap.add_argument("--out", required=True, help="genotype matrix, IID x locus, codes 1/2/3 (participant level)")
    ap.add_argument("--counts", default=None, help="per-locus DEL/DUP counts, masked (aggregate)")
    ap.add_argument("--min_cell", type=int, default=20)
    a = ap.parse_args()
    s = pd.read_csv(a.samples, sep="\t", dtype=str)
    samples = list(s["IID"] if "IID" in s.columns else pd.read_csv(a.samples, header=None, dtype=str)[0])
    cnv = pd.read_csv(a.cnv, sep=r"\s+", dtype={"FID": str, "IID": str})
    cnv["CHR"] = pd.to_numeric(cnv.CHR.astype(str).str.replace("chr", "", regex=False), errors="coerce")
    G = genotype(cnv.dropna(subset=["CHR"]).astype({"CHR": int}), samples, read_loci(a.reciprocal), read_loci(a.oneway))
    G.to_csv(a.out, sep="\t")
    if a.counts:
        counts(G, a.min_cell).to_csv(a.counts, sep="\t", index=False)
    print(f"genotyped {len(G)} people at {sum(not c.endswith('_oneway') for c in G.columns)} loci -> {a.out}")


if __name__ == "__main__":
    main()
