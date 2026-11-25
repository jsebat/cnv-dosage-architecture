# Rules for coding agents working in this repo

## Reuse the lab's code. Do not reinvent it.

Before writing any new analysis code, search this repo for an existing implementation and use it. Extend it if needed, and say so in
the commit message. Never write a second version of something that already exists. In particular:

- **Recurrent CNV genotyping = `recurrentCNVgenotyping/` only** (Molly Sacks's `GenotypeRecurrentCNVs_tutorial.ipynb`, and its port
  `recurrentCNVgenotyping/aou/recurrent_cnv_genotyping_aou.py`). Her locus sheet (reciprocal and one-way tabs), her overlap rules and her
  segment patterns define who carries a recurrent CNV. In Python, use `from common import rcg` (in `AllOfUs/cnv_anthro` or
  `cnv_gwas/cnv_anthro`), then:
  - `rcg.genotype_calls(cnv)`: person × allele matrix, coded 1 = DEL, 2 = none, 3 = DUP;
  - `rcg.carriers(G, names, dosage)`: carrier IIDs;
  - `rcg.label_calls(cnv, loci)`: per-call locus label, used for locus-removal (LOO) runs.

  Never write reciprocal-overlap, coverage or bedtools rules for recurrent loci anywhere else (Jon Sebat, 2026-11-24). If a locus you need
  is not in Molly's sheet, ask Jon. Do not add a local rule.

  The one exception is `alleles.py` / `cnv_gwas/replication/replicate_s9_alleles.R`. It is a port of Omar Shanta's PGC Table S9 allele
  calls, used only to replicate that table. Omar's rules stay for the Shanta (PGC) paper (Jon 2026-11-24).
- **CNV → gene rules**: `incidence()` in `conga_axes_qt.py` (DEL: any exon; DUP: a whole transcript).
- **Nulls**: genome shift of the gene annotation (`conga_grid_qt.py`, shared shifts in `resources/shared_shifts`), and the gene-circle
  shift for the gene mirror (`gene_mirror_qt.py`).
- **Meta-analysis**: `cnv_gwas/cnv_anthro/meta_sumstats.py`, `meta_grid.py`, `AllOfUs/meta/`.
- **Pipelines**: `AllOfUs/cnv_anthro/run_pipeline.py` and `run_fig3.py`. Add a step there rather than writing a standalone script.

## Data

- Participant-level data never goes into this public repo. All of Us exports follow the small-cell rule: counts 1–20 are withheld.
- Every script starts with a header comment saying what it does, and a line `# null model: ...` (or `# null model: n/a`).
