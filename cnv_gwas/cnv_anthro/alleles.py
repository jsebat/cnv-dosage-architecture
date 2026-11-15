"""Allele-level genotypes at the Table S9 loci -- a port of Omar Shanta's allele calling
(scripts_figures/scripts_trumpet_noSex: 01/02 drivers, make_indiv.py, anno_exons.py, make_indiv_lof.py,
CNV_genotype_combine.py), plus the quantitative-trait association.

Rules (bedtools semantics with PLINK BP1/BP2 used as start/end, overlap = min(end) - max(start)):
  loci_cytoband.txt : reciprocal overlap >= 50% (bedtools -r -f 0.5)
  BP_cytoband.txt   : call covers >= 50% of the locus (bedtools -F 0.5)
  single genes      : DEL contains >= 1 whole exon of the gene; DUP contains >= 1 whole isoform (bedtools -F 1)
  NRXN1, ASTN2 DEL  : loss-of-function DELs only (frameshift by exon phase, or TSS) -- anno_exons.py logic
  compound alleles  : CNV_genotype_combine.py (e.g. 22q11.21 A-B + B-D -> A-D; 15q11.2 BP1-BP2 dropped if BP1-BP3)
A person carries an allele if any of their calls is assigned to it (make_indiv.py).
"""
import argparse
import json
import os
import numpy as np
import pandas as pd
from common import read_plink_cnv, design_matrix, Residualiser

FRAME = {1: 2, 2: 3, 3: 1}


def _ov(s1, e1, s2, e2):
    return np.maximum(0, np.minimum(e1, e2) - np.maximum(s1, s2))


def carriers_by_rule(calls, loci, rule):
    """dict locus NAME -> set of person rows. calls: DataFrame CHR BP1 BP2 row; loci: CHR START END NAME."""
    out = {}
    for L in loci.itertuples():
        c = calls[calls.CHR == L.CHR]
        o = _ov(c.BP1.to_numpy(), c.BP2.to_numpy(), L.START, L.END)
        if rule == "recip":
            m = (o > 0) & (o >= 0.5 * (c.BP2 - c.BP1).to_numpy()) & (o >= 0.5 * (L.END - L.START))
        elif rule == "coverF":
            m = (o > 0) & (o >= 0.5 * (L.END - L.START))
        elif rule == "contains":            # feature (exon / isoform) entirely inside the call (bedtools -F 1)
            m = (o > 0) & (o >= (L.END - L.START))
        out.setdefault(L.NAME, set()).update(c.row.to_numpy()[m])
    return out


def lof_dels(dels, exon_txt, exon_bed):
    """anno_exons.py: classify each DEL that contains >= 1 whole element of the exon table; returns {(row,BP1,BP2): type}."""
    ex = pd.read_csv(exon_txt, sep="\t", header=None, names=["CHR", "BP1", "BP2", "Exon", "Acceptor", "Donor"], na_values=["NA"])
    tab = pd.read_csv(exon_bed, sep="\t", na_values=["NA"])
    hits = []
    for e in ex.itertuples():
        c = dels[dels.CHR == e.CHR]
        o = _ov(c.BP1.to_numpy(), c.BP2.to_numpy(), e.BP1, e.BP2)
        m = (o > 0) & (o >= (e.BP2 - e.BP1))
        for r in c[m].itertuples():
            hits.append((r.IID, r.BP1, r.BP2, r.row, e.BP1, e.Exon))
    if not hits:
        return {}
    h = pd.DataFrame(hits, columns=["IID", "BP1", "BP2", "row", "eBP1", "Exon"])
    h = h.sort_values(["IID", "BP1", "eBP1"], ascending=[True, True, False], kind="mergesort")   # sort -k5,5 -k2,2n -k7,7nr
    out = {}
    for (iid, bp1), g in h.groupby(["IID", "BP1"], sort=False):
        es, ee = g.Exon.iloc[0], g.Exon.iloc[-1]
        si = int(tab.index[tab.Exon == es][0]); ei = int(tab.index[tab.Exon == ee][0])
        don = tab.Donor.iloc[si - 1] if si - 1 >= 0 else np.nan
        acc = tab.Acceptor.iloc[ei + 1] if ei + 1 < len(tab) else np.nan
        if es == "TSS":
            t = "TSS"
        elif np.isnan(don) or np.isnan(acc):
            t = "Frameshift"
        elif FRAME[int(don)] == int(acc):
            t = "InFrame"
        else:
            t = "Frameshift"
        out[(g.row.iloc[0], bp1, g.BP2.iloc[0])] = t
    return out


def call_alleles(cnv, n_rows, res_dir):
    """cnv: calls with column 'row' (person index 0..n_rows-1). Returns DataFrame n_rows x allele columns (0/1)."""
    rd = lambda f: pd.read_csv(os.path.join(res_dir, f), sep="\t")
    all_l = rd("all_cnvs.txt"); cyto = rd("loci_cytoband.txt"); bpc = rd("BP_cytoband.txt")
    sg = rd("gws_loci_singleGene.txt"); sg_ex = rd("gws_loci_singleGene_exons.txt"); sg_iso = rd("gws_loci_singleGene_isoforms.txt")
    G = {}
    for tn, tc in [("DEL", 1), ("DUP", 3)]:
        calls = cnv[cnv.TYPE == tc]
        for name, s in carriers_by_rule(calls, cyto, "recip").items():
            G[f"{name}_{tn}"] = s
        for name, s in carriers_by_rule(calls, bpc, "coverF").items():
            G[f"{name}_{tn}"] = s
        feat = sg_ex if tn == "DEL" else sg_iso
        for name, s in carriers_by_rule(calls, feat, "contains").items():
            G[f"{name}_{tn}"] = s
    dels = cnv[cnv.TYPE == 1]
    for gene in ["nrxn1", "astn2"]:
        a = lof_dels(dels, os.path.join(res_dir, f"{gene}_exons.txt"), os.path.join(res_dir, f"{gene}_exons.bed"))
        b = lof_dels(dels, os.path.join(res_dir, f"{gene}_tss_exons.txt"), os.path.join(res_dir, f"{gene}_tss_exons.bed"))
        merged = dict(a); merged.update(b)                       # mergeExonsTSS: the TSS-table call wins when both exist
        G[f"{gene.upper()}_DEL"] = {k[0] for k, t in merged.items() if t in ("Frameshift", "TSS")}
    D = pd.DataFrame(index=range(n_rows))
    for k, s in G.items():
        v = np.zeros(n_rows, int); v[list(s)] = 1; D[k] = v
    for nm in list(all_l.NAME + "_DUP") + list(all_l.NAME + "_DEL"):
        if nm not in D:
            D[nm] = 0
    L = D  # literal port of CNV_genotype_combine.py below
    b = lambda c: L[c].astype(bool)
    for T in ["DEL", "DUP"]:
        i = b(f"1q21.1_TAR_BP1-BP2_{T}") & b(f"1q21.1_TAR_BP2-BP3_{T}")
        L.loc[i, f"1q21.1_TAR_BP1-BP2_{T}"] = 0; L.loc[i, f"1q21.1_TAR_BP2-BP3_{T}"] = 0; L.loc[i, f"1q21.1_TAR_BP1-BP3_{T}"] = 1
        i = b(f"1q21.1_TAR_BP1-BP3_{T}") & b(f"1q21.1_{T}")
        L[f"1q21.1_collapsed_{T}"] = L[f"1q21.1_{T}"]
        L.loc[i, f"1q21.1_TAR_BP1-BP3_{T}"] = 0; L.loc[i, f"1q21.1_{T}"] = 0
        i = b(f"15q11.2_BP1-BP2_{T}") & b(f"15q11.2-13.1_BP1-BP3_{T}"); L.loc[i, f"15q11.2_BP1-BP2_{T}"] = 0
        i = b(f"15q13.1-13.2_BP3-BP4_{T}") & b(f"15q13.3_BP4-BP5_{T}")
        L.loc[i, f"15q13.1-13.2_BP3-BP4_{T}"] = 0; L.loc[i, f"15q13.3_BP4-BP5_{T}"] = 0; L.loc[i, f"15q13.1-13.3_BP3-BP5_{T}"] = 1
        i = b(f"15q13.3_BP4-BP5_{T}") & b(f"15q13.3_CHRNA7_{T}"); L.loc[i, f"15q13.3_CHRNA7_{T}"] = 0
        i = b(f"15q11.2-13.1_BP1-BP3_{T}") & b(f"15q13.1-13.3_BP3-BP5_{T}")
        L[f"15q11.2-13.1_BP1-BP3_collapsed_{T}"] = L[f"15q11.2-13.1_BP1-BP3_{T}"]
        L.loc[i, f"15q11.2-13.1_BP1-BP3_{T}"] = 0; L.loc[i, f"15q13.1-13.3_BP3-BP5_{T}"] = 0
        i = b(f"16p13.11_BP1-BP2_{T}") & b(f"16p13.11_BP2-BP3_{T}")
        L.loc[i, f"16p13.11_BP1-BP2_{T}"] = 0; L.loc[i, f"16p13.11_BP2-BP3_{T}"] = 0; L.loc[i, f"16p13.11_BP1-BP3_{T}"] = 1
        L[f"16p11.2_BP2-BP3_collapsed_{T}"] = L[f"16p13.11_BP2-BP3_{T}"] | L[f"16p13.11_BP1-BP3_{T}"]
        i = b(f"16p11.2_BP1-BP2_{T}") & b(f"16p11.2_BP2-BP3_{T}")
        L.loc[i, f"16p11.2_BP1-BP2_{T}"] = 0; L.loc[i, f"16p11.2_BP2-BP3_{T}"] = 0; L.loc[i, f"16p11.2_BP1-BP3_{T}"] = 1
        L[f"16p13.11_BP1-BP2_collapsed_{T}"] = L[f"16p13.11_BP1-BP2_{T}"] | L[f"16p13.11_BP1-BP3_{T}"]
        i = b(f"22q11.21_B-C_{T}") & b(f"22q11.21_C-D_{T}")
        L.loc[i, f"22q11.21_B-C_{T}"] = 0; L.loc[i, f"22q11.21_C-D_{T}"] = 0; L.loc[i, f"22q11.21_B-D_{T}"] = 1
        i = b(f"22q11.21_A-B_{T}") & b(f"22q11.21_B-D_{T}")
        L.loc[i, f"22q11.21_A-B_{T}"] = 0; L.loc[i, f"22q11.21_B-D_{T}"] = 0; L.loc[i, f"22q11.21_A-D_{T}"] = 1
        L[f"22q11.21_A-B_collapsed_{T}"] = L[f"22q11.21_A-B_{T}"] | L[f"22q11.21_A-D_{T}"]
        for nm in [f"22q11.21_F-G_{T}", f"22q11.21_G-H_{T}", f"22q11.21_F-H_{T}"]:
            if nm not in L:
                L[nm] = 0
        i = b(f"22q11.21_F-G_{T}") & b(f"22q11.21_G-H_{T}")
        L.loc[i, f"22q11.21_F-G_{T}"] = 0; L.loc[i, f"22q11.21_G-H_{T}"] = 0; L.loc[i, f"22q11.21_F-H_{T}"] = 1
        L[f"USP7-HAPSTR1_{T}"] = L[f"USP7_{T}"] & L[f"HAPSTR1_{T}"]
        L[f"USP7-HAPSTR1_collapsed_{T}"] = L[f"USP7_{T}"] | L[f"HAPSTR1_{T}"]
    L.columns = ["X_" + c.replace("-", ".") for c in L.columns]
    return L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno", required=True); ap.add_argument("--cnv", required=True)
    ap.add_argument("--resources", required=True, help="dir with Omar's locus files (resources/omar_alleles)")
    ap.add_argument("--traits", default="height_rint,bmi_rint"); ap.add_argument("--out", required=True)
    ap.add_argument("--numeric", default="age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10")
    ap.add_argument("--factors", default="sex,batch")
    ap.add_argument("--s9_alleles", default=None, help="Table S9 allele list (NAME, Locus_Type) to report; default all_cnvs.txt x DEL/DUP")
    ap.add_argument("--omar_matrix", default=None, help="optional: Omar's genotype_allele.indiv for the same samples (validation)")
    ap.add_argument("--omar_fid", default="NA", help="FID used in Omar's matchID (UKBB: NA)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    ph = pd.read_csv(a.pheno, sep="\t", dtype={"IID": str})
    cnv = read_plink_cnv(a.cnv, keep_iids=ph["IID"])
    cnv = cnv[cnv.CHR <= 23].copy()
    cnv["row"] = pd.Series(np.arange(len(ph)), index=ph["IID"]).loc[cnv.IID].to_numpy()
    L = call_alleles(cnv, len(ph), a.resources)
    L.insert(0, "IID", ph["IID"].to_numpy())
    L.to_csv(os.path.join(a.out, "allele_genotypes.tsv.gz"), sep="\t", index=False)
    info = {}
    if a.omar_matrix:
        om = pd.read_csv(a.omar_matrix, sep="\t", dtype={"matchID": str})
        om["IID"] = om.matchID.str.split(":", n=1).str[1]
        m = L.merge(om, on="IID", suffixes=("", "_omar"))
        cmp = []
        for c in [c for c in L.columns if c != "IID" and c in om.columns]:
            x = m[c].astype(int); yv = (m[c + "_omar"] > 0).astype(int)
            cmp.append((c, int(x.sum()), int(yv.sum()), int((x != yv).sum())))
        cmp = pd.DataFrame(cmp, columns=["allele", "carriers_port", "carriers_omar", "discordant"])
        cmp.to_csv(os.path.join(a.out, "validation_vs_omar_matrix.tsv"), sep="\t", index=False)
        info["validation_people"] = int(len(m)); info["validation_alleles"] = int(len(cmp))
        info["validation_alleles_discordant"] = int((cmp.discordant > 0).sum()); info["validation_total_discordant_calls"] = int(cmp.discordant.sum())
    rows = []
    num = [c for c in a.numeric.split(",") if c]; fac = [c for c in a.factors.split(",") if c]
    if a.s9_alleles:
        s9 = pd.read_csv(a.s9_alleles, sep="\t")
        alle = ["X_" + f"{r.NAME}_{r.Locus_Type}".replace("-", ".") for r in s9.itertuples()]
    else:
        al = pd.read_csv(os.path.join(a.resources, "all_cnvs.txt"), sep="\t")
        alle = ["X_" + f"{n}_{t}".replace("-", ".") for t in ["DUP", "DEL"] for n in al.NAME]
    from scipy import stats as st
    for tr in a.traits.split(","):
        keep = np.isfinite(ph[tr]).to_numpy() & ph[num + fac].notna().all(axis=1).to_numpy()
        P = ph[keep].reset_index(drop=True); G = L[keep].reset_index(drop=True)
        Z, _ = design_matrix(P, num, fac); R = Residualiser(Z)
        yr = R.resid(P[tr].to_numpy(float)); yy = yr @ yr
        raw = P[tr.replace("_rint", "")].to_numpy(float) if tr.endswith("_rint") else None
        for al in alle:
            if al not in G:
                rows.append((tr, al, 0, np.nan, np.nan, np.nan, np.nan, np.nan, len(P))); continue
            idx = np.flatnonzero(G[al].to_numpy() > 0)
            b, se, t, p = R.test_carriers(idx, yr, yy) if len(idx) >= 2 else (np.nan,) * 4
            rawmean = np.nanmean(raw[idx]) - np.nanmean(raw) if raw is not None and len(idx) else np.nan
            rows.append((tr, al, len(idx), b, se, t, p, rawmean, len(P)))
    out = pd.DataFrame(rows, columns=["trait", "allele", "n_carriers", "BETA_SD", "SE", "T", "P", "raw_mean_diff_unadjusted", "N"])
    out["FREQ"] = out.n_carriers / out.N                     # carrier frequency (array vs WGS sensitivity comparison)
    out["model"] = "OLS trait_rint ~ allele + " + "+".join(num + fac)
    out.to_csv(os.path.join(a.out, "allele_assoc.tsv"), sep="\t", index=False, float_format="%.6g")
    json.dump(info, open(os.path.join(a.out, "allele_info.json"), "w"), indent=1)
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
