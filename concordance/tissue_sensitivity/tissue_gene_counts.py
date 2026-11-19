# Genes inside the panel-A CNVs (non-overlapping set without Williams, contiguous TWAS-gene spans as in panel A) that have a TWAS30
# expression model in each tissue (single) and in at least one of the other tissues (all but that tissue). Output: tissue_gene_counts.csv
import os, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(f"{HERE}/tissue_drivers.py").read(); src = src[:src.index("def run(")] + src[src.index("GROUP = "):src.index("rows = []")]
ns = {}; exec(src, ns)
W, s0, Ls, k = ns["W"], ns["s0"], ns["Ls"], ns["kAll"]
idx = np.unique(np.concatenate([np.arange(s, s + L) for s, L in zip(s0[k], Ls[k])]))
H = W.iloc[idx].notna()
out = pd.DataFrame({"n_genes_single": H.sum(), "n_genes_all_but": [int(H.drop(columns=t).any(axis=1).sum()) for t in H.columns]})
out.loc["ALL"] = [int(H.any(axis=1).sum()), np.nan]
out.to_csv(f"{HERE}/tissue_gene_counts.csv"); print(out.sort_values("n_genes_single").to_string())
