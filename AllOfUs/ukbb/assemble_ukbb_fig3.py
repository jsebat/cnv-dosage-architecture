"""UK Biobank Fig. 3 export from the finished run folders (no refit): lays the four UKBB runs out as run_fig3.py expects and writes the config
that points at them. The export itself is then `run_fig3.py --config <work>/ukbb_fig3_config.json --steps export --min_cell 0`.

  python assemble_ukbb_fig3.py --results /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/results --config ../cnv_anthro/configs/ukbb_expanse.json --work <new dir>

Runs used (Jon 2026-11-23): size2t_mirror_2026-11-23, all2t_mirror_2026-11-23, totdose_2026-11-23, gene_mirror_2026-11-22.
min_cell 0: UK Biobank summary statistics stay at UCSD, so no small-cell suppression; the gene mirror keeps its own rule (>= 5 carriers).
null model: none here (files are linked, nothing is computed).
"""
import argparse
import json
import os

RUNS = {"size2t_mirror": "size2t_mirror_2026-11-23", "all2t_mirror": "all2t_mirror_2026-11-23", "totdose": "totdose_2026-11-23",
        "gene_mirror": "gene_mirror_2026-11-22"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True); ap.add_argument("--config", required=True); ap.add_argument("--work", required=True)
    a = ap.parse_args()
    root = os.path.join(a.work, "results", "all", "fig3"); os.makedirs(root)
    for name, src in RUNS.items():
        s = os.path.join(a.results, src); assert os.path.isdir(s), s
        os.symlink(s, os.path.join(root, name))
    C = json.load(open(a.config))
    C["results_dir"] = os.path.join(a.work, "results"); C["fig3"]["export_dir"] = os.path.join(a.work, "export")
    C["cohort"] = "UKBB"; C["note"] = "UK Biobank Fig. 3 export from the finished 2026-11-22/23 runs (assemble_ukbb_fig3.py)"
    json.dump(C, open(os.path.join(a.work, "ukbb_fig3_config.json"), "w"), indent=1)
    json.dump({k: os.path.join(a.results, v) for k, v in RUNS.items()}, open(os.path.join(a.work, "ASSEMBLY_SOURCES.json"), "w"), indent=1)
    print("assembled", root)


if __name__ == "__main__":
    main()
