"""Rare pLoF SNV/indel carriers per gene in All of Us short-read WGS, for gene_assoc_qt.py --plof (Jon 2026-11-25: "Gene level stats for
LOF DEL and LOF SNV", to show convergence of WGS/WES and CNV signals). Runs INSIDE the Workbench (Hail, Dataproc); the output is
person-level and never leaves the Workbench.

  python plof_carriers_aou.py --vat <VAT tsv(.bgz)> --mt <srWGS exome split Hail MatrixTable> --people <pheno tsv with IID> \
         --out plof_carriers.tsv.gz [--max_af 0.01] [--consequences stop_gained,frameshift_variant,splice_acceptor_variant,splice_donor_variant]

Variant set (from the Variant Annotation Table, VAT): a canonical-transcript consequence in --consequences (high-confidence protein-truncating
classes; the same classes as the exome pLoF masks of Kosmicki et al. 2026 and Akbari et al. 2021, without LOFTEE flags, which the VAT does
not carry) and a cohort alternate-allele frequency (VAT gvs_all_af) below --max_af. Carriers: people in --people with a non-reference
genotype (het or hom-alt) at any such variant of the gene, after the callset's own genotype filters. Output columns: IID, symbol (one row per
person x gene). The gene symbol is matched to the pipeline's exon annotation by gene_assoc_qt.py.
Check before running: the VAT column names used below (vid, gene_symbol, consequence, is_canonical_transcript, gvs_all_af) and the MT/VAT
paths of the current CDR (the Workbench's genomics data page lists them); the script stops if a column is missing.
null model: n/a (variant extraction; no test).
"""
import argparse
import hail as hl

COLS = ["vid", "gene_symbol", "consequence", "is_canonical_transcript", "gvs_all_af"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vat", required=True); ap.add_argument("--mt", required=True); ap.add_argument("--people", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--max_af", type=float, default=0.01)
    ap.add_argument("--consequences", default="stop_gained,frameshift_variant,splice_acceptor_variant,splice_donor_variant")
    a = ap.parse_args()
    hl.init(default_reference="GRCh38", idempotent=True)
    cons = hl.set(a.consequences.split(","))
    v = hl.import_table(a.vat, force_bgz=a.vat.endswith(".bgz"), impute=False, missing="")
    miss = [c for c in COLS if c not in v.row]
    if miss:
        raise SystemExit(f"VAT lacks columns {miss}; set them to this CDR's names")
    v = v.filter((hl.str(v.is_canonical_transcript).lower() == "true")
                 & hl.any(lambda c: cons.contains(c), v.consequence.split(","))
                 & (hl.float64(v.gvs_all_af) < a.max_af) & (v.gene_symbol != ""))
    p = v.vid.split("-")                                                 # VAT vid = chrom-pos-ref-alt (chrom without 'chr')
    v = v.annotate(locus=hl.locus("chr" + p[0], hl.int32(p[1])), alleles=[p[2], p[3]]).key_by("locus", "alleles")
    v = v.select("gene_symbol").distinct()
    people = hl.import_table(a.people, key="IID", types={"IID": hl.tstr})
    mt = hl.read_matrix_table(a.mt)
    mt = mt.filter_cols(hl.is_defined(people[hl.str(mt.s)]))
    n_people = people.count(); n_in = mt.count_cols()                   # people without srWGS would be miscounted as non-carriers
    print(f"people in --people: {n_people}; with short-read WGS genotypes: {n_in} ({n_in / max(n_people, 1):.1%})")
    if n_in < n_people:
        print("WARNING: not every analysed person has srWGS; PLOF/LOF for this dataset would count them as non-carriers. "
              "Use --plof only on a dataset where this is 100%, or restrict the analysis to people with srWGS.")
    mt = mt.semi_join_rows(v)
    mt = mt.annotate_rows(symbol=v[mt.row_key].gene_symbol)
    gt = mt.GT if "GT" in mt.entry else hl.vds.lgt_to_gt(mt.LGT, mt.LA)
    keep = gt.is_non_ref()
    if "FT" in mt.entry:                                                  # genotype filter: PASS (or missing) only
        keep = keep & (hl.is_missing(mt.FT) | (hl.str(mt.FT) == "PASS") | (hl.str(mt.FT) == "true")) if mt.FT.dtype != hl.tbool else keep & ~mt.FT
    mt = mt.filter_entries(keep)
    e = mt.entries().select("symbol")
    e = e.key_by().select(IID=hl.str(e.s), symbol=e.symbol).distinct()
    e.export(a.out)
    n = e.count(); g = len(e.aggregate(hl.agg.collect_as_set(e.symbol)))
    print(f"pLoF carrier x gene pairs: {n}; genes: {g}; variants classes {a.consequences}; max AF {a.max_af}; output stays in the Workbench: {a.out}")


if __name__ == "__main__":
    main()
