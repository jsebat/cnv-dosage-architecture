"""Cross-check recurrent_cnv_genotyping_aou.py against Molly's notebook (run cell by cell with bedtools) on random calls."""
import json, os, subprocess, sys, numpy as np, pandas as pd
os.chdir(sys.argv[1]); rng = np.random.default_rng(7)
loci = pd.concat([pd.read_csv(f, sep='\t') for f in ['recurrentCNV_hg38.reciprocal.tsv', 'recurrentCNV_hg38.oneway.tsv']])
rows = []; N = 4000
for i in range(N):
    for _ in range(rng.integers(0, 4)):
        L = loci.iloc[rng.integers(len(loci))]; span = L.END - L.START
        if rng.random() < .3:   # whole multi-segment events
            L2 = loci[(loci.CHR == L.CHR)].sample(1, random_state=int(rng.integers(1e9))).iloc[0]
            s, e = min(L.START, L2.START), max(L.END, L2.END)
        else:
            s, e = L.START, L.END
        s = int(s + rng.normal(0, .3) * span); e = int(e + rng.normal(0, .3) * span)
        if e - s < 1000: continue
        rows.append((f'chr{L.CHR}', 0, 0, 0, 'x', 'IP,UP', s, e, rng.choice(['DEL', 'DUP']), f'S{i}'))
plt = pd.DataFrame(rows, columns=['chrom', 'locusStart', 'locusEnd', 'locusSize', 'locusType', 'algorithm', 'consStart', 'consEnd', 'consCallType', 'sampleId'])
plt.to_csv('yourData.CNVcalls.txt', sep='\t', index=False)
pd.DataFrame({'sampleID': [f'S{i}' for i in range(N)], 'QC': 'pass'}).to_csv('yourData_manifest.txt', sep='\t', index=False)
# PLINK version of the same calls for the port
pd.DataFrame({'FID': plt.sampleId, 'IID': plt.sampleId, 'CHR': plt.chrom.str[3:], 'BP1': plt.consStart, 'BP2': plt.consEnd,
              'TYPE': plt.consCallType.map({'DEL': 1, 'DUP': 3}), 'SCORE': 0, 'SITES': 0}).to_csv('calls.cnv', sep='\t', index=False)
pd.DataFrame({'IID': [f'S{i}' for i in range(N)]}).to_csv('samples.tsv', sep='\t', index=False)
# Molly's notebook, cell by cell (only the Google URLs and the bedtools download replaced)
nb = json.load(open('GenotypeRecurrentCNVs_tutorial.ipynb')); g = {}
for i, c in enumerate(nb['cells']):
    if c['cell_type'] != 'code' or i in (1,): continue
    src = ''.join(c['source'])
    src = '\n'.join(l for l in src.splitlines() if not l.strip().startswith(("!head", "url_", "sheet_id", "gid")))
    src = src.replace("url_reciprocal", "'recurrentCNV_hg38.reciprocal.tsv'").replace("url_oneway", "'recurrentCNV_hg38.oneway.tsv'")
    if src.startswith('%%bash'):
        env = dict(os.environ, COHORT=g.get('COHORT', 'yourData'), TODAY=g.get('TODAY', 'x'))
        subprocess.run(src.split('\n', 1)[1], shell=True, check=True, env=env, executable='/bin/bash', capture_output=True)
    else:
        src = src.replace("!awk -F'\\t' '{print $4}' recurrentCNV_hg38.reciprocal.bed > loci_names.txt",
                          "subprocess.run(\"awk -F'\\\\t' '{print $4}' recurrentCNV_hg38.reciprocal.bed > loci_names.txt\", shell=True)")
        src = src.replace("!awk -F'\\t' '{print $4 \"_oneway\"}' recurrentCNV_hg38.oneway.bed >> loci_names.txt",
                          "subprocess.run(\"awk -F'\\\\t' '{print $4 \\\"_oneway\\\"}' recurrentCNV_hg38.oneway.bed >> loci_names.txt\", shell=True)")
        g['subprocess'] = subprocess
        exec(src, g)
M = g['genomatrix'].set_index('identifier')
subprocess.run([sys.executable, sys.argv[2], '--cnv', 'calls.cnv', '--samples', 'samples.tsv', '--reciprocal', 'recurrentCNV_hg38.reciprocal.tsv',
                '--oneway', 'recurrentCNV_hg38.oneway.tsv', '--out', 'port.tsv'], check=True)
P = pd.read_csv('port.tsv', sep='\t', index_col=0)
cols = [c for c in P.columns if c in M.columns]
diff = {c: int((P[c].astype(float) != M.loc[P.index, c].astype(float)).sum()) for c in cols}
print('calls', len(plt), '| loci compared', len(cols), '| port-only', sorted(set(P.columns) - set(M.columns))[:5], '| molly-only', sorted(set(M.columns) - set(P.columns) - {'sampleID', 'QC', '15q13', '15q35', '16p11.2', '22q11.21_AD', '22q11.21_DH'})[:8])
print('carriers per locus (molly):', int(sum(((M[c] == 1) | (M[c] == 3)).sum() for c in cols)))
print('mismatching genotypes:', {k: v for k, v in diff.items() if v} or 'none')
# explain each mismatch: does the person have both a DEL and a DUP (or overlapping calls) at that locus?
cnv = pd.read_csv('calls.cnv', sep='\t', dtype={'IID': str})
L = pd.concat([pd.read_csv(f, sep='\t') for f in ['recurrentCNV_hg38.reciprocal.tsv', 'recurrentCNV_hg38.oneway.tsv']]).set_index('Locus name')
why = {'both DEL and DUP at the locus/region': 0, 'other': 0}; other = []
for c in [k for k, v in diff.items() if v]:
    for iid in P.index[(P[c].astype(float) != M.loc[P.index, c].astype(float)).to_numpy()]:
        base = c.replace('_oneway', '')
        reg = L.loc[base] if base in L.index else None
        x = cnv[cnv.IID == iid]
        if reg is not None: x = x[(x.CHR == reg.CHR) & (x.BP2 > reg.START) & (x.BP1 < reg.END)]
        if x.TYPE.nunique() > 1: why['both DEL and DUP at the locus/region'] += 1
        else: why['other'] += 1; other.append((c, iid, P.loc[iid, c], M.loc[iid, c], x[['BP1', 'BP2', 'TYPE']].values.tolist()))
print(why); print(other[:6])
