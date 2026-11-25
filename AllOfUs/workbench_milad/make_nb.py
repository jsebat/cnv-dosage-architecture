"""Build sacks2026_aou_workbench.ipynb: a restartable notebook for the Sacks 2026 All of Us run (HANDOFF_Milad.md Part B)."""
import json, sys, uuid

cells = []
def md(s): cells.append({"cell_type": "markdown", "id": uuid.uuid4().hex[:8], "metadata": {}, "source": s.strip("\n")})
def code(s): cells.append({"cell_type": "code", "id": uuid.uuid4().hex[:8], "metadata": {}, "execution_count": None,
                           "outputs": [], "source": s.strip("\n")})

md(r"""
# Sacks 2026 — All of Us run (HANDOFF_Milad.md, Part B)

**Restartable.** After the app is stopped or recreated, run every cell from the top, in order, up to the "Status" cell. Each cell skips what is already done:
- the repo is cloned only if missing and is then updated;
- the bundle is unpacked only if `resources/` is missing;
- the run folder (`~/cnv_anthro_run`) is restored from the bucket backup if the local disk was wiped;
- the pipeline stages run in the background (`nohup`) and leave a `.done` marker. `launch()` skips finished stages, and every finished stage is backed up to the bucket.

**Keep the notebook itself in the bucket** (`workspace-bucket/cnv_growth_bmi/`) so it survives the app.

**Egress rule.** No participant count from 1 to 20 may leave the Workbench, and that includes pasting into chat. The QC cells below print counts 1–20 as `<=20`. Paste only their output, never raw tables.
""")

md("## 1. Settings (the only cell to edit)")
code(r"""
import os, re, sys, json, glob, hashlib, subprocess, time, shutil
import pandas as pd

REPO_URL   = 'https://github.com/jsebat/cnv-dosage-architecture.git'   # Jon 2026-11-19: code moved here (same history)
BRANCH     = 'main'
PIN_COMMIT = '2e638c1'         # Jon 2026-11-24: one recurrent-CNV genotyper (Molly's, recurrentCNVgenotyping/aou), rerun in section 17; run_fig3 audit fixed. Before: bdfcf29 (Jon 2026-11-23: Fig. 3 joint model + gene-based association (HANDOFF_Milad_Fig3.md). Before: 4b64136 (Jon 2026-11-18: sex strata fix (numeric match). ce1ce67 _nosd fix (additivity rerun done). sv/prep/main(all,eur) c4e09df; gridjoint 05fa65b; burden/main_rest b97150f)

HOME = os.path.expanduser('~')                 # /home/dataproc on the Spark app, /home/jupyter on a plain JupyterLab app
REPO = f'{HOME}/repos/cnv-growth-bmi'
ROOT = f'{REPO}/AllOfUs'
CODE = f'{ROOT}/cnv_anthro'
RUN  = os.path.expanduser('~/cnv_anthro_run')      # = results_dir/export_dir parent in the config (local disk, fast)

BUCKET_MNT = f'{HOME}/workspace/workspace-bucket/cnv_growth_bmi'          # bucket as mounted on the VM
BUCKET_GS  = os.environ['WORKSPACE_BUCKET'].rstrip('/') + '/cnv_growth_bmi'     # same folder, gs:// form
BACKUP_GS  = f'{BUCKET_GS}/run_backup'                                          # mirror of ~/cnv_anthro_run

BUNDLE     = 'cnv_resources_sacks2026_2026-11-16.tar.gz'
BUNDLE_SHA = '27c8b6dc5c66b4cf52c07bc3d195b36cf0658fc1395294f5409215142d9df44b'

# --- CDR v9 genomic data (mounted data collection; found 2026-11-17) ---
V9 = f'{HOME}/workspace/cdrv9/vwb-aou-datasets-controlled-v9/v9/wgs/short_read'
SV_VCF_PATTERN = f'{V9}/structural_variants/vcf/full/AoU_srWGS_SV.v9.chr{{c}}.vcf.gz'   # one GATK-SV VCF per chromosome
ANCESTRY_SRC   = f'{V9}/snpindel/aux/ancestry/ancestry_preds.tsv'
SV_PARALLEL    = 22            # chromosomes converted at once

CFG_TEMPLATE = 'configs/aou_sacks2026.json'          # tracked in git, never edited
CFG          = 'configs/aou_sacks2026_filled.json'   # untracked copy with the two paths filled
MIN_CELL     = 20              # set 39 if the dissemination article asks for >= 40, then rerun export
GRID_JOBS    = 24

GP = os.environ.get('GOOGLE_PROJECT')
def sh(cmd, check=True):
    print('$', cmd); r = subprocess.run(cmd, shell=True, executable='/bin/bash')
    if check and r.returncode: raise RuntimeError(f'exit {r.returncode}: {cmd}')
    return r.returncode
for d in ['data', 'results', 'logs', 'stages']:
    os.makedirs(f'{RUN}/{d}', exist_ok=True)
print('CDR      ', os.environ.get('WORKSPACE_CDR'))
print('bucket   ', BUCKET_GS, '| mounted:', os.path.isdir(os.path.dirname(BUCKET_MNT)))
print('project  ', GP)
""")

md("## 2. Code: clone or update the repo")
code(r"""
if not os.path.isdir(f'{REPO}/.git'):
    os.makedirs(os.path.dirname(REPO), exist_ok=True)
    sh(f'git clone -b {BRANCH} {REPO_URL} {REPO}')
else:
    # local edits (e.g. notebooks opened inside the repo) are stashed, not lost: `git stash list` shows them
    sh(f'cd {REPO} && (git diff --quiet && git diff --cached --quiet || git stash push -m "auto-stash $(date +%F_%T)")')
    # repo moved 2026-11-19 (cnv-growth-bmi:cnv-gwas -> cnv-dosage-architecture:main, same commits): repoint an old clone
    sh(f'cd {REPO} && [ "$(git remote get-url origin)" = "{REPO_URL}" ] || git remote set-url origin {REPO_URL}')
    sh(f'cd {REPO} && git fetch -q origin && git checkout -q -B {BRANCH} origin/{BRANCH}')
if PIN_COMMIT:
    sh(f'cd {REPO} && git checkout -q {PIN_COMMIT}')
COMMIT = subprocess.check_output(['git', '-C', REPO, 'log', '-1', '--format=%h %cs %s']).decode().strip()
print('code at:', COMMIT)
open(f'{RUN}/logs/code_commits.txt', 'a').write(f'{time.strftime("%F %T")}\t{COMMIT}\n')
sh(f'cd {REPO} && git merge-base --is-ancestor dec364d HEAD && echo "OK: at or after dec364d"')
""")

md("## 3. VM check\nNeeds 64–96 vCPUs, ≥128 GB RAM, ≥200 GB free disk. **Set the app's auto-pause / idle timeout to the maximum**: background jobs do not count as activity, and a paused VM kills them.")
code(r"""
sh('nproc; free -g | head -2; df -h $HOME | tail -1; python3 --version')
sh('python3 -c "import numpy, scipy, pandas, matplotlib; print(\'numpy\', numpy.__version__, \'scipy\', scipy.__version__, \'pandas\', pandas.__version__)"')
# BigQuery client (phenotype extraction in the prep steps). On the JupyterLab app (2026-11-20) the preinstalled
# google-cloud-bigquery was newer than google-api-core ("google.api_core has no attribute check_python_version"): upgrade both.
bq = subprocess.run(f'{sys.executable} -c "from google.cloud import bigquery; print(\'bigquery\', bigquery.__version__)"', shell=True, capture_output=True, text=True)
if bq.returncode != 0:
    print('BigQuery import failed:', bq.stderr.strip().splitlines()[-1] if bq.stderr.strip() else '?', '-> upgrading')
    sh(f'{sys.executable} -m pip install -q -U google-api-core google-cloud-bigquery db-dtypes')
    bq = subprocess.run(f'{sys.executable} -c "from google.cloud import bigquery; print(\'bigquery\', bigquery.__version__)"', shell=True, capture_output=True, text=True)
print(bq.stdout.strip() or bq.stderr.strip()[-500:])
""")

md("## 4. Resources bundle (public annotations + UKBB shift lists)")
code(r"""
if not os.path.exists(f'{ROOT}/resources/shared_shifts/SHIFTS_MANIFEST.tsv'):
    local = f'{ROOT}/{BUNDLE}'
    if os.path.exists(f'{BUCKET_MNT}/{BUNDLE}'):
        shutil.copy(f'{BUCKET_MNT}/{BUNDLE}', local)
    else:
        sh(f'gsutil cp {BUCKET_GS}/{BUNDLE} {local}')
    h = hashlib.sha256(open(local, 'rb').read()).hexdigest()
    assert h == BUNDLE_SHA, f'sha256 mismatch: {h}'
    print('sha256 OK')
    sh(f'cd {ROOT} && tar xzf {BUNDLE} && rm -f {BUNDLE}')
print(sorted(os.listdir(f'{ROOT}/resources')))
man = pd.read_csv(f'{ROOT}/resources/shared_shifts/SHIFTS_MANIFEST.tsv', sep='\t', dtype=str)
print(man.to_string(index=False))
assert not man.apply(lambda c: c.str.strip().eq('NO')).any().any(), 'a shift list says NO: stop and tell Jon'
""")

md("## 4b. Restore the run folder from the bucket backup (no-op on a fresh run)")
code(r"""
def backup(run=RUN, gs=BACKUP_GS):
    sh(f'gsutil -m -q rsync -r {run} {gs}')
# every dataset's run folder + the array inputs (same paths as DATASETS in section 7)
RESTORE = [(RUN, BACKUP_GS),
           (os.path.expanduser('~/cnv_anthro_nosegdup'), f'{BUCKET_GS}/run_backup_nosegdup'),
           (os.path.expanduser('~/cnv_anthro_array'), f'{BUCKET_GS}/run_backup_array'),
           (os.path.expanduser('~/cnv_anthro_array/inputs'), f'{BUCKET_GS}/array_inputs')]
for run_, gs_ in RESTORE:
    have_backup = subprocess.run(f'gsutil -q ls {gs_}/', shell=True, capture_output=True).returncode == 0
    # restore ONLY when the local copy is empty/new (no finished stage, no inputs): on a live VM the older backup would overwrite newer local logs
    has_local = bool(glob.glob(f'{run_}/stages/*.done')) or os.path.exists(f'{run_}/BUILD_INFO.txt')
    if have_backup and not has_local:
        os.makedirs(run_, exist_ok=True); sh(f'gsutil -m -q rsync -r {gs_} {run_}')   # adds/updates local files; deletes nothing
        print('restored', run_)
    elif have_backup:
        print('local copy has content, restore skipped:', run_)
    print('  stages done:', sorted(os.path.basename(p)[:-5] for p in glob.glob(f'{run_}/stages/*.done')) or 'none')
""")

md("## 5. Data paths: check that the v9 SV VCFs (chr1–22) and ancestry_preds.tsv are visible")
code(r"""
miss = [c for c in range(1, 23) if not os.path.exists(SV_VCF_PATTERN.format(c=c))]
print('SV VCFs chr1-22:', 'all present' if not miss else f'MISSING {miss}')
print('ancestry_preds.tsv:', os.path.exists(ANCESTRY_SRC))
sh(f'ls -la {os.path.dirname(SV_VCF_PATTERN)} | head -5; du -ch {os.path.dirname(SV_VCF_PATTERN)}/*chr[0-9]*.vcf.gz | tail -1', check=False)
assert not miss and os.path.exists(ANCESTRY_SRC)
""")

md("## 6. Fill the config (writes an untracked copy; the tracked template is not touched)\nAoU v9 ships the SV VCF per chromosome, so the SV → CNV conversion runs as its own stage (`sv`: 22 chromosomes in parallel, same script and filters as the config, then concatenated). The PGC 10% rule is computed within chromosome and over the same sample count, so this equals a single-file run. The config's `prep` then only extracts the phenotypes.")
code(r"""
anc = f'{RUN}/data/ancestry_preds.tsv'
if not os.path.exists(anc):
    if ANCESTRY_SRC.startswith('gs://'):
        sh(f'gsutil -u {GP} cp {ANCESTRY_SRC} {anc}')
    else:
        shutil.copy(ANCESTRY_SRC, anc)
import sys
cfg = json.load(open(f'{CODE}/{CFG_TEMPLATE}'))
p0, p1 = cfg['prep_commands']
assert p0[0] == 'aou_sv_to_cnv.py' and p1[0] == 'aou_pheno_extract.py'
p1[p1.index('--ancestry') + 1] = anc
cfg['prep_commands'] = [p1]                       # SV conversion = stage 'sv' (per chromosome)
SV_OUT = p0[p0.index('--out') + 1].replace('${HOME}', os.path.expanduser('~'))
SV_SAMPLES = p0[p0.index('--write_samples') + 1].replace('${HOME}', os.path.expanduser('~'))
SV_ARGS = ' '.join(f'{p0[i]} {p0[i + 1]}' for i in range(1, len(p0), 2)
                   if p0[i] not in ('--vcf', '--out', '--write_samples', '--chroms'))
SVD = f'{RUN}/data/sv_by_chr'
sv_sh = f'''#!/bin/bash
# AoU v9 SV VCF per chromosome -> one CNV table (aou_sv_to_cnv.py, filters from the config: {SV_ARGS})
set -euo pipefail
mkdir -p {SVD}; cd {CODE}
one() {{
  c=$1
  [ -s {SVD}/chr$c.cnv ] && {{ echo "chr$c done"; return 0; }}
  {sys.executable} aou_sv_to_cnv.py --vcf {SV_VCF_PATTERN.format(c='$c')} --out {SVD}/chr$c.cnv.tmp \\
     --write_samples {SVD}/samples_chr$c.txt --chroms $c-$c {SV_ARGS} > {SVD}/chr$c.log 2>&1
  mv {SVD}/chr$c.cnv.tmp {SVD}/chr$c.cnv; echo "chr$c ok: $(tail -1 {SVD}/chr$c.log)"
}}
export -f one
seq 1 22 | xargs -P {SV_PARALLEL} -I{{}} bash -c 'one {{}}'
# every chromosome must list the same samples in the same order
for c in $(seq 2 22); do cmp -s {SVD}/samples_chr1.txt {SVD}/samples_chr$c.txt || {{ echo "sample list differs: chr$c"; exit 1; }}; done
cp {SVD}/samples_chr1.txt {SV_SAMPLES}
head -1 {SVD}/chr1.cnv > {SV_OUT}.tmp
for c in $(seq 1 22); do tail -n +2 {SVD}/chr$c.cnv >> {SV_OUT}.tmp; done
mv {SV_OUT}.tmp {SV_OUT}
echo "samples $(wc -l < {SV_SAMPLES}); calls $(( $(wc -l < {SV_OUT}) - 1 ))"
'''
open(f'{RUN}/logs/sv_stage.sh', 'w').write(sv_sh)
cfg['min_cell'] = MIN_CELL
assert 'FILL_FROM_DATA_DICTIONARY' not in json.dumps(cfg['prep_commands'])
json.dump(cfg, open(f'{CODE}/{CFG}', 'w'), indent=1)
shutil.copy(f'{CODE}/{CFG}', f'{RUN}/logs/config_used.json')
print('sv stage ', f'{RUN}/logs/sv_stage.sh', '| filters:', SV_ARGS)
print('cnv table', SV_OUT); print('ancestry ', anc, f'({sum(1 for _ in open(anc)) - 1} rows)')
print('min_cell ', cfg['min_cell']); print('config   ', f'{CODE}/{CFG}')
""")

md(r"""
## 7. Run the pipeline in the background

| stage | steps | time |
|---|---|---|
| `sv` | 22 SV VCFs → CNV table (parallel; finished chromosomes are kept on restart) | ~1–3 h |
| `prep` | OMOP → phenotypes (BigQuery) | < 1 h |
| `main` | gwas, plots, alleles, loci, burden (all 9 strata) | hours |
| `additivity` | additivity, scores (stratum all) | minutes |
| `gridjoint` | grid, joint (stratum all), 1,000–2,000 nulls | several hours |
| `export` | small-cell suppression → `export_aggregate/` (run **after QC**, section 9) | minutes |

`launch()` runs the listed stages in order, skips `.done` ones, stops at the first failure, and backs up `~/cnv_anthro_run` to the bucket after each stage. The jobs survive closing the browser or the kernel, but not a stopped VM.
""")
code(r"""
STAGES = {
    'sv':         None,                      # bash logs/sv_stage.sh (section 6)
    'prep':       '--steps prep',
    'main':       '--steps gwas,plots,alleles,loci,burden',   # FAILED at afr on c4e09df; all+eur gwas/alleles/loci valid. Superseded by the next two
    'burden':     '--steps burden',                                                   # b97150f: all 9 strata (WGS single_def changed)
    'main_rest':  '--steps gwas,plots,alleles,loci --strata afr,amr,eas,sas,mid,female,male',
    'sexstrata':  '--steps burden,gwas,plots,alleles,loci --strata female,male',   # 4b64136: female/male were N = 0 (sex 1.0/2.0)
    'additivity': '--steps additivity,scores --strata all',
    'gridjoint':  f'--steps grid,joint --strata all --jobs {GRID_JOBS}',
    'export':     '--steps export',
    'burden_allcalls': 'py:burden_allcalls.py',
    # Fig. 3 of the dosage paper + gene-based association (Jon 2026-11-23, HANDOFF_Milad_Fig3.md; section 16). Always --strata all for fig3.
    'fig3':         '--steps fig3 --strata all --jobs 5',      # 14 fits + 2 gene-mirror + 2 dose-bin runs + in-Workbench summaries
    'fig3_export':  '--steps fig3_export --strata all',        # aggregate-only -> <run>/export_aggregate_fig3, log ends with AUDIT OK
    'genes':        '--steps genes',                           # every stratum, both traits -> results/<stratum>/gene_assoc/<trait>
    'genes_export': '--steps genes_export --strata all',       # aggregate-only -> <run>/export_aggregate_genes, log ends with AUDIT OK
    # 2e638c1 (Jon 2026-11-24): every step with recurrent-CNV carriers now uses Molly's genotyping; rerun into new results folders (section 17).
    # Jon's 'run_pipeline.py --steps loci,grid,joint,additivity,burden', split as before so --jobs only reaches grid/joint (it also sets the
    # number of parallel additivity runs). loci + burden on every stratum, additivity on 'all'. fig3 above = run_fig3.py fit,summary.
    'rcg':          '--steps loci,additivity,burden',
    'rcg_gridjoint': f'--steps grid,joint --strata all --jobs {GRID_JOBS}',
}
# Three datasets, each with its own run folder, filled config and bucket backup (stage markers live in <run>/stages):
#   wgs      = default WGS run (B2)          nosegdup = WGS without the segdup GWAS filter (B6)          array = Molly's array calls (B5)
DATASETS = {
    'wgs':      dict(run=RUN, cfg=CFG, backup=BACKUP_GS, tar='aou_sacks2026',
                     stages=['sv', 'prep', 'burden', 'main_rest', 'sexstrata', 'additivity', 'gridjoint', 'export', 'fig3', 'genes', 'fig3_export', 'genes_export']),   # 'main' (FAILED at afr, c4e09df) superseded by burden/main_rest/sexstrata
    'nosegdup': dict(run=os.path.expanduser('~/cnv_anthro_nosegdup'), cfg='configs/aou_sacks2026_nosegdup_filled.json',
                     backup=f'{BUCKET_GS}/run_backup_nosegdup', tar='aou_sacks2026_nosegdup', stages=['main', 'additivity', 'gridjoint', 'export']),
    'array':    dict(run=os.path.expanduser('~/cnv_anthro_array'), cfg='configs/aou_sacks2026_array_filled.json',
                     backup=f'{BUCKET_GS}/run_backup_array', tar='aou_array_sacks2026', stages=['prep', 'main', 'burden_allcalls', 'additivity', 'gridjoint', 'export', 'fig3', 'genes', 'fig3_export', 'genes_export']),
}
for D in DATASETS.values():
    for d in ['data', 'results', 'logs', 'stages']: os.makedirs(f"{D['run']}/{d}", exist_ok=True)
# drivers: tag '' = main, 'b' = second WGS driver, 'ns' = nosegdup, 'ar' = array. A driver runs one dataset; pid files stay in ~/cnv_anthro_run/logs
def pidf(tag=''): return f'{RUN}/logs/driver{tag}.pid'

def running(tag=''):
    if not os.path.exists(pidf(tag)): return None
    try:
        pid = int(open(pidf(tag)).read()); os.kill(pid, 0); return pid
    except (ValueError, ProcessLookupError, PermissionError):
        return None

def launch(stages, tag='', ds='wgs'):
    assert not running(tag), f'driver{tag} already running (pid {running(tag)}); see status()'
    PIDF = pidf(tag); D = DATASETS[ds]; RUN_, CFG_, BK_ = D['run'], D['cfg'], D['backup']
    assert os.path.exists(f'{CODE}/{CFG_}'), f'filled config missing: {CFG_} (run its fill cell first)'
    import sys
    lines = ['#!/bin/bash', f'echo $$ > {PIDF}', f'cd {CODE}', 'export PYTHONUNBUFFERED=1', f'export NUMEXPR_MAX_THREADS={os.cpu_count()}']
    for st in stages:
        lines += [f'if [ -e {RUN_}/stages/{st}.done ]; then echo "skip {ds}:{st} (done)"; else',
                  f'  rm -f {RUN_}/stages/{st}.failed; echo "=== {ds}:{st} start $(date)"; t0=$(date +%s)',
                  (f'  bash {RUN_}/logs/sv_stage.sh > {RUN_}/logs/{st}.log 2>&1' if STAGES[st] is None else
                   f'  {sys.executable} {RUN_}/logs/{STAGES[st][3:]} {CFG_} > {RUN_}/logs/{st}.log 2>&1' if STAGES[st].startswith('py:') else
                   f'  {sys.executable} run_pipeline.py --config {CFG_} {STAGES[st]} > {RUN_}/logs/{st}.log 2>&1'),
                  f'  rc=$?; echo "=== {ds}:{st} exit $rc after $(( ($(date +%s)-t0)/60 )) min $(date)"',
                  f'  if [ $rc -ne 0 ]; then touch {RUN_}/stages/{st}.failed; exit $rc; fi',
                  f'  echo "{COMMIT}" > {RUN_}/stages/{st}.done',
                  f'  gsutil -m -q rsync -r {RUN_} {BK_} >> {RUN_}/logs/backup.log 2>&1 && echo "backup ok"',
                  'fi']
    open(f'{RUN}/logs/driver{tag}.sh', 'w').write('\n'.join(lines) + '\n')
    if os.path.exists(PIDF): os.remove(PIDF)
    subprocess.Popen(f'nohup setsid bash {RUN}/logs/driver{tag}.sh >> {RUN}/logs/driver{tag}.log 2>&1 &', shell=True)
    time.sleep(2)
    print('started', ds, stages, 'driver', repr(tag), 'pid', running(tag), '\nwatch with status()')

def status(n=15):
    TAGS = ['', 'b', 'ns', 'ar']
    print('drivers:', ' | '.join(f"{t or 'main'}: {running(t) or '-'}" for t in TAGS))
    for ds, D in DATASETS.items():
        marks = {os.path.basename(p).rsplit('.', 1)[0]: p.rsplit('.', 1)[1] for p in glob.glob(f"{D['run']}/stages/*.done") + glob.glob(f"{D['run']}/stages/*.failed")}
        print(f'  {ds:9s}', '  '.join(f"{st}:{'DONE' if marks.get(st) == 'done' else 'FAILED' if marks.get(st) == 'failed' else '-'}" for st in D['stages']))
    sh(' '.join(f"tail -n 3 {RUN}/logs/driver{t}.log 2>/dev/null;" for t in TAGS), check=False)
    logs = sorted(sum((glob.glob(f"{D['run']}/logs/*.log") for D in DATASETS.values()), []), key=os.path.getmtime)
    logs = [l for l in logs if not os.path.basename(l).startswith(('driver', 'backup'))]
    # egress rule: raw logs contain model tables with small counts, so only progress/command/error lines are shown
    keep = re.compile(r'^\[ *[0-9.]+s\]|^=+ |^\+ |^== |Error|Traceback|FAILED|failed strata')
    for l in logs[-3:]:
        lines = [x.rstrip() for x in open(l, errors='replace') if keep.search(x)]
        print(f'\n--- {os.path.relpath(l, os.path.expanduser("~"))} (last {n} progress lines)'); print('\n'.join(x[:300] for x in lines[-n:]))
    sh('uptime; free -g | sed -n 2p', check=False)
""")
code(r"""
# WGS is finished (sv, prep, burden, main_rest, sexstrata, additivity, gridjoint). Do NOT relaunch: 'main' has no .done
# (it failed and was replaced by burden + main_rest), so launching it would rerun the whole WGS main step over the restored results.
# Original first launch, kept for the record:
# launch(['sv', 'prep', 'main', 'additivity', 'gridjoint'])
""")
code(r"""
status()      # rerun any time
""")
code(r"""
# manual backup any time (also runs automatically after each stage)
backup()
""")

md(r"""
## 7b. Pinning a new commit mid-run (e.g. Jon's b97150f)
1. Set `PIN_COMMIT` in section 1, then rerun sections **1, 2, 4 and 6**. Section 6 must be rerun: the filled config is regenerated from the new template. A running stage is not affected when only scripts it does not call changed (b97150f changes `dup_burden.py` and `run_pipeline.py`, not the grid/joint scripts), but check `git diff --stat` first.
2. Run only what Jon asked for. `burden` and `main_rest` write `burden/` and the afr…male folders, never `all/grid` or `all/joint`, so they may run as a second driver while `gridjoint` runs (cell below).
""")
code(r"""
# ce1ce67 (_nosd reference fix): keep the pre-fix additivity results, then rerun additivity+scores (stratum all)
PRE = f'{RUN}/prefix_nosd_20261117'
if not os.path.exists(PRE):
    shutil.copytree(f'{RUN}/results/all/additivity', f'{PRE}/additivity'); shutil.copytree(f'{RUN}/results/all/scores', f'{PRE}/scores')
    print('pre-fix copy ->', PRE)
MARK = f'{RUN}/stages/additivity_rerun_ce1ce67.started'
if not os.path.exists(MARK):                     # one-time: done 2026-11-18, so this cell is a no-op on a rerun
    if os.path.exists(f'{RUN}/stages/additivity.done'): os.remove(f'{RUN}/stages/additivity.done')
    open(MARK, 'w').write('ce1ce67\n'); launch(['additivity'], tag='b')
else:
    print('ce1ce67 additivity rerun already done/started:', os.path.exists(f'{RUN}/stages/additivity.done') and open(f'{RUN}/stages/additivity.done').read().strip())
""")
code(r"""
# sex coding check: female/male strata came out N = 0 (sex written as 1.0/2.0?)
ph = pd.read_csv(f'{RUN}/data/aou_pheno_std.tsv', sep='\t', usecols=['sex'], dtype=str)
print('sex values as written:', mask(ph.sex.value_counts().to_dict()))
""")
code(r"""
# b97150f: burden in all strata, then gwas/plots/alleles/loci in afr..male (second driver, alongside gridjoint)
launch(['burden', 'main_rest'], tag='b')
""")
code(r"""
# 4b64136 (Jon 2026-11-18): make_stratum matches sex numerically, so the existing aou_pheno_std.tsv works (no prep rerun).
# Reruns burden + gwas/plots/alleles/loci for female and male only; other strata untouched. Expect N = 57,965 / 37,058 (before trait filters).
assert COMMIT.startswith('4b64136') or PIN_COMMIT == '4b64136', 'pin 4b64136 first (sections 1, 2, 4, 6)'
launch(['sexstrata'], tag='b')
""")
code(r"""
# check after sexstrata finishes: stratum N and the outputs exist (aggregate only)
R = f'{RUN}/results'
for st in ['female', 'male']:
    sh(f"grep -h '== stratum {st}' {RUN}/logs/sexstrata.log | tail -1", check=False)
    for t in ['height', 'bmi']:
        fs = [f'{R}/{st}/taskA_gwas/{t}', f'{R}/{st}/taskA_alleles/{t}/allele_assoc.tsv', f'{R}/{st}/burden/{t}/run_info.json', f'{R}/{st}/locus_effects/locus_effects_{t}.tsv']
        print(f'  {st} {t}:', ['ok' if os.path.exists(f) else 'MISSING ' + os.path.relpath(f, R) for f in fs])
""")
md("## 8. QC (B3) — outputs are small-cell masked and safe to paste\nSections 8–10 run on one dataset at a time: set `DS` in the next cell (`'wgs'`, `'nosegdup'` or `'array'`) and rerun 8–10. 8.6 is WGS only.")
code(r"""
def mask(o):
    # counts 1..20 -> '<=20' (egress rule); applied to every int in the printed output
    if isinstance(o, dict): return {k: mask(v) for k, v in o.items()}
    if isinstance(o, list): return [mask(v) for v in o]
    if isinstance(o, bool): return o
    if isinstance(o, (int,)) and 1 <= o <= 20: return '<=20'
    if isinstance(o, float) and o.is_integer() and 1 <= o <= 20: return '<=20'
    return o
def flat(d, pre=''):
    out = {}
    for k, v in d.items():
        if isinstance(v, dict): out.update(flat(v, f'{pre}{k}.'))
        elif not isinstance(v, list): out[f'{pre}{k}'] = v
    return out
DS = 'wgs'                       # 'wgs' | 'nosegdup' | 'array'
DRUN = DATASETS[DS]['run']; R = f'{DRUN}/results'; print('QC/export dataset:', DS, '->', R)
def show_info(pattern, keys=r'(^|\.)(N|n_|n$|carrier|union|asym|hit_)', maxk=25):
    import re
    for f in sorted(glob.glob(f'{R}/{pattern}', recursive=True)):
        d = flat(json.load(open(f)))
        d = {k: v for k, v in d.items() if re.search(keys, k)}
        print(os.path.relpath(f, R), json.dumps(mask(dict(list(d.items())[:maxk]))))

# 8.0 phenotype summary from prep (already suppressed by the extractor)
p = f'{R}/aou_pheno_summary.json'
if os.path.exists(p): print(json.dumps(mask(json.load(open(p))), indent=1))
""")
code(r"""
# 8.1 sample sizes: GWAS and burden per stratum, and the carrier subsets
show_info('*/taskA_gwas/*/run_info.json')
show_info('*/burden/*/run_info.json')
show_info('*/joint/strata/union_strata_info.json')
""")
code(r"""
# 8.2 the 'only' additivity runs -> send Jon reference_N and locus_carriers_in_sample (main: no_rare_gene; sensitivity: _no_rare_gene_nosd)
rows = []
for f in sorted(glob.glob(f'{R}/all/additivity/*/only*/run_info.json')):
    d = flat(json.load(open(f))); parts = os.path.relpath(f, R).split('/')
    g = lambda k: next((v for kk, v in d.items() if kk == k or kk.endswith('.' + k)), None)   # keys are nested in run_info
    rows.append(dict(trait=parts[2], run=parts[3], N=g('N'), reference_N=g('reference_N'),
                     locus_carriers_in_sample=g('locus_carriers_in_sample'), only_ref=g('only_ref')))
T_only = pd.DataFrame(rows)
for c_ in ['N', 'reference_N', 'locus_carriers_in_sample']:
    if c_ in T_only: T_only[c_] = T_only[c_].map(lambda v: mask(v) if isinstance(v, (int, float)) and v == v else v)
print(T_only.to_string(index=False))
""")
code(r"""
# 8.2b signed-burden panel (b97150f): send Jon single_dosage_N for all/height (with reference_N from 8.2)
for st in ['all']:
    for t in ['height', 'bmi']:
        f = f'{R}/{st}/burden/{t}/run_info.json'
        if not os.path.exists(f): print('missing', f); continue
        d = flat(json.load(open(f)))
        print(st, t, json.dumps(mask({k: v for k, v in d.items() if re.search(r'single|min_n|^N$', k)})))
""")
code(r"""
# 8.3 hit carriers in the additivity runs: hit_carriers_DEL / _DUP must be > 20
import re
for f in sorted(glob.glob(f'{R}/all/additivity/**/run_info.json', recursive=True)):
    d = flat(json.load(open(f)))
    hc = {k: v for k, v in d.items() if re.search(r'hit_carriers', k)}
    bad = [k for k, v in hc.items() if isinstance(v, (int, float)) and v <= 20]
    print('FLAG ' if bad else 'ok   ', os.path.relpath(f, R), json.dumps(mask(hc)))
""")
code(r"""
# 8.4 recurrent-CNV carrier frequencies (EUR and all): paste this output, it is compared with UKBB outside the Workbench
REC = r'15q11\.2_BP1\.BP2_(DEL|DUP)$|16p11\.2_BP4\.BP5_(DEL|DUP)$|16p11\.2_BP1\.BP4|22q11\.21_A\.D_(DEL|DUP)$|1q21\.1|17q12'
for st in ['eur', 'all']:
    f = f'{R}/{st}/taskA_alleles/height/allele_assoc.tsv'
    if not os.path.exists(f): print('missing', f); continue
    d = pd.read_csv(f, sep='\t').drop_duplicates('allele')
    d = d[d.allele.str.contains(REC)][['allele', 'n_carriers', 'N', 'FREQ']].copy()
    small = d.n_carriers.between(1, 20)
    d['n_carriers'] = d.n_carriers.astype(object); d.loc[small, 'n_carriers'] = '<=20'; d.loc[small, 'FREQ'] = None
    print(f'--- {st}'); print(d.to_string(index=False))
""")
code(r"""
# 8.5 Miami plots (probe-level, aggregate): genomic inflation per stratum
from IPython.display import Image, display
for f in sorted(glob.glob(f'{R}/*/taskA_gwas/*/miami_*.png')):
    print(os.path.relpath(f, R)); display(Image(f, width=900))
""")
code(r"""
# 8.6 1q21.1 TAR DUP breakpoint check (Jon 2026-11-18). All DUP calls in chr1:145.3-146.3 Mb, people in the phenotype table.
# Aggregate only: sites with <= 20 carriers are pooled, counts 1-20 masked, no margins, totals rounded to 10 (no masked cell
# can be recovered by subtraction). Safe to paste.
r10 = lambda n: f'~{int(round(n, -1))}' if n > 20 else mask(int(n))
import numpy as np
SD = f'{CODE}/../resources/segDups.bed'
TAR = {'BP1-BP2': (145627235, 145809100), 'BP2-BP3': (145823900, 146040039), 'BP1-BP3': (145627235, 146040039)}
ids = set(pd.read_csv(f'{RUN}/data/aou_pheno_std.tsv', sep='\t', usecols=['IID'], dtype=str).IID)
cnv = pd.read_csv(f'{RUN}/data/aou_srwgs_sv.cnv', sep='\t', dtype={'IID': str})
c = cnv[(cnv.CHR == 1) & (cnv.TYPE == 3) & (cnv.BP2 > 145_300_000) & (cnv.BP1 < 146_300_000) & cnv.IID.isin(ids)].copy()
print('DUP calls in window:', r10(len(c)), '| carriers:', r10(c.IID.nunique()), '| distinct sites (BP1,BP2):', c.groupby(['BP1', 'BP2']).ngroups)

sd = pd.read_csv(SD, sep='\t', header=None, usecols=[0, 1, 2], names=['chr', 's', 'e'])
sd = sd[sd.chr.astype(str).isin(['chr1', '1'])]
mid = ((c.BP1 + c.BP2) // 2).to_numpy()
c['mid_in_segdup'] = [bool(((sd.s <= m) & (sd.e >= m)).any()) for m in mid]
seg_cov = lambda a, b: (lambda x: x[x[:, 0] < x[:, 1]])(np.c_[np.maximum(sd.s.to_numpy(), a), np.minimum(sd.e.to_numpy(), b)])
def frac_sd(a, b):            # share of the call covered by the union of segdups
    iv = seg_cov(a, b)
    if not len(iv): return 0.0
    iv = iv[np.argsort(iv[:, 0])]; tot, cs, ce = 0, iv[0, 0], iv[0, 1]
    for x, y in iv[1:]:
        if x > ce: tot += ce - cs; cs, ce = x, y
        else: ce = max(ce, y)
    return (tot + ce - cs) / (b - a)
ov = lambda a, b, s, e: max(0, min(b, e) - max(a, s))
def labels(a, b):
    L = b - a
    recip = [k for k, (s, e) in TAR.items() if ov(a, b, s, e) >= 0.5 * L and ov(a, b, s, e) >= 0.5 * (e - s)]
    cover = [k for k, (s, e) in TAR.items() if k != 'BP1-BP3' and ov(a, b, s, e) >= 0.5 * (e - s)]   # alleles.py rule (coverF)
    pipe = 'BP1-BP3' if set(cover) == {'BP1-BP2', 'BP2-BP3'} else (cover[0] if cover else '-')
    return ','.join(recip) or 'none', pipe

# site types from the chr1 SV VCF via its .tbi index (pysam; the VM has no tabix binary). Records are keyed by POS = BP1.
vcf = SV_VCF_PATTERN.format(c=1); st = {}
try:
    try:
        import pysam
    except ImportError:
        subprocess.run(f'{sys.executable} -m pip install -q pysam', shell=True, check=True); import pysam
    tb = pysam.TabixFile(vcf); ctg = 'chr1' if 'chr1' in tb.contigs else '1'
    lo, hi = int(c.BP1.min()) - 1, int(c.BP1.max()) + 1
    for ln in tb.fetch(ctg, lo, hi):
        f = ln.split('\t', 8); info = dict(x.split('=', 1) for x in f[7].split(';') if '=' in x)
        st[(int(f[1]), int(info.get('END', f[1])))] = f"{info.get('SVTYPE')}/{f[6]}"
    print('VCF records read:', len(st))
except Exception as e:
    print('VCF lookup skipped:', repr(e))

S = c.groupby(['BP1', 'BP2']).agg(carriers=('IID', 'nunique'), mid_in_segdup=('mid_in_segdup', 'first')).reset_index()
S['kb'] = (S.BP2 - S.BP1) // 1000
S[['recip50', 'pipeline_label']] = [labels(a, b) for a, b in zip(S.BP1, S.BP2)]
S['segdup_frac'] = [round(frac_sd(a, b), 2) for a, b in zip(S.BP1, S.BP2)]
S['vcf_type'] = [st.get((a, b), '?') for a, b in zip(S.BP1, S.BP2)]

print('\n--- by pipeline label (alleles.py coverF rule, before person-level collapsing) x 50% reciprocal match')
g = c.merge(S[['BP1', 'BP2', 'recip50', 'pipeline_label', 'vcf_type']], on=['BP1', 'BP2'])
print(mask(pd.crosstab([g.pipeline_label, g.recip50], g.vcf_type).to_dict('index')))
print('\n--- midpoint in a segdup, by pipeline label (calls)'); print(mask(pd.crosstab(g.pipeline_label, g.mid_in_segdup).to_dict('index')))

print('\n--- start / end distribution (25 kb bins, Mb; calls)')
for col in ['BP1', 'BP2']:
    b = (c[col] // 25000 * 25000 / 1e6).round(3).value_counts().sort_index()
    print(col, {float(k): int(v) for k, v in b.items() if v > 20}, '| other bins pooled:', r10(int(b[b <= 20].sum())))

print('\n--- sites with >= 21 carriers (others pooled below)')
big = S[S.carriers > 20].sort_values('BP1')
print(big[['BP1', 'BP2', 'kb', 'carriers', 'recip50', 'pipeline_label', 'mid_in_segdup', 'segdup_frac', 'vcf_type']].to_string(index=False))
sm = S[S.carriers <= 20]
print(f'pooled small sites: {len(sm)} sites, carriers', r10(int(sm.carriers.sum())), '| by pipeline label:',
      {k: r10(int(v)) for k, v in sm.groupby('pipeline_label').carriers.sum().items()})
""")

code(r"""
# 8.7 ARRAY ONLY: are recurrent DUPs/DELs split or truncated? (2026-11-20: array DUP carrier frequencies at 22q11.21 A-D,
# 1q21.1, TAR BP1-BP3 and 17q12 are well below UKBB, and TAR BP1-BP2 DUP is above). alleles.py needs ONE call covering >= 50%
# of the locus; consStart/consEnd is where both callers agree, so a call can be shorter than either caller's.
# Per locus and type, people in the phenotype table with any call overlapping >= 10% of the locus, by:
#   best1 = coverage by their single best call, union = coverage by all their calls in the locus merged; ncalls = calls in locus.
# Aggregate only: counts <= 20 masked, others rounded to 10, no totals. Safe to paste.
assert DS == 'array', 'set DS = "array" in section 8 first'
r10 = lambda n: f'~{int(round(n, -1))}' if n > 20 else mask(int(n))
LOCI = {'22q11.21_A-D': (22, 19049819, 21112437), '1q21.1': (1, 147056425, 147922330), 'TAR_BP1-BP3': (1, 145627235, 146040039),
        'TAR_BP1-BP2': (1, 145627235, 145809100), '17q12': (17, 36460090, 37857812), '15q11.2_BP1-BP2': (15, 22778537, 23067754),
        '16p11.2_BP4-BP5 (control)': (16, 29639519, 30189452)}
ids = set(pd.read_csv(f'{DRUN}/data/aou_array_pheno_std.tsv', sep='\t', usecols=['IID'], dtype=str).IID)
cnv = pd.read_csv(f'{DRUN}/data/aou_array.cnv', sep=r'\s+', dtype={'IID': str})
cnv = cnv[cnv.IID.isin(ids)]
def union_len(iv):
    iv = sorted(iv); tot, cs, ce = 0, iv[0][0], iv[0][1]
    for x, y in iv[1:]:
        if x > ce: tot += ce - cs; cs, ce = x, y
        else: ce = max(ce, y)
    return tot + ce - cs
rows = []
for name, (ch, a, b) in LOCI.items():
    L = b - a
    for typ, tname in [(1, 'DEL'), (3, 'DUP')]:
        c = cnv[(cnv.CHR == ch) & (cnv.TYPE == typ) & (cnv.BP2 > a) & (cnv.BP1 < b)].copy()
        if not len(c): rows.append(dict(locus=name, type=tname, people_ge10pct='0')); continue
        c['s'], c['e'] = c.BP1.clip(lower=a), c.BP2.clip(upper=b)
        g = pd.DataFrame([{'best1': (d.e - d.s).max() / L, 'union': union_len(list(zip(d.s, d.e))) / L, 'ncalls': len(d)} for _, d in c.groupby('IID')])
        g = g[g.union >= 0.10]
        if not len(g): rows.append(dict(locus=name, type=tname, people_ge10pct='0')); continue
        rows.append(dict(locus=name, type=tname, people_ge10pct=r10(len(g)),
                         best1_ge50=r10((g.best1 >= .5).sum()), best1_20_50=r10(((g.best1 >= .2) & (g.best1 < .5)).sum()),
                         best1_lt20=r10((g.best1 < .2).sum()), union_ge50=r10((g.union >= .5).sum()),
                         split_2plus_calls=r10((g.ncalls >= 2).sum()), median_best1_cov=round(float(g.best1.median()), 2) if len(g) > 20 else 'masked'))
print(pd.DataFrame(rows).fillna('').to_string(index=False))
# raw PGC call columns (header only): are there per-caller or union coordinates besides consStart/consEnd?
RAW = os.path.expanduser('~/cnv_anthro_array/molly_raw')          # section 12's raw copy (restored from the array backup)
plt_ = sorted(glob.glob(f'{RAW}/4_filtered/*.plt'))[:1]
print('\nraw .plt header:', open(plt_[0]).readline().strip().split('\t') if plt_ else f'no .plt under {RAW}/4_filtered')
""")

md("## 9. Export (after QC; dataset = `DS` from section 8). Check the in-Workbench article *Downloading and Disseminating Genomic Data*: if it asks for ≥ 40, set `MIN_CELL = 39` in section 1, rerun sections 1 and 6, delete `stages/export.done`, and tell Jon.")
code(r"""
launch(['export'], tag={'wgs': '', 'nosegdup': 'ns', 'array': 'ar'}[DS], ds=DS)
""")
code(r"""
status()
E = f'{DRUN}/export_aggregate'
if os.path.exists(f'{E}/EXPORT_MANIFEST.txt'): print(open(f'{E}/EXPORT_MANIFEST.txt').read()[:3000])
""")
code(r"""
# 9.1 audit: no count from 1 to 20 in any count-like column of any exported table
import re
COUNT = re.compile(r'^(n|N)$|^n_|_n$|^N_|carrier|count|^ncnv|^n[A-Z]', re.I)
# not person counts (2026-11-20 audit): probes per call, z statistics of the no-loci runs, mean/SD of per-person feature counts
NOT_COUNT = re.compile(r'^n_probes$|^noloci_|^(mean|sd)_count$', re.I)
bad, nfile = [], 0
for f in glob.glob(f'{E}/**/*.tsv', recursive=True) + glob.glob(f'{E}/**/*.csv', recursive=True):
    nfile += 1
    d = pd.read_csv(f, sep='\t' if f.endswith('.tsv') else ',', dtype=str, low_memory=False)
    for c in [c for c in d.columns if COUNT.search(str(c)) and not NOT_COUNT.search(str(c))]:
        v = pd.to_numeric(d[c], errors='coerce')
        k = int(v.between(1, 20).sum())
        if k: bad.append((os.path.relpath(f, E), c, k))
print(nfile, 'tables checked')
print('\n'.join(f'VIOLATION {f} [{c}] {k} rows' for f, c, k in bad) or 'OK: no 1-20 counts in count columns')
sh(f'find {E} -type f | grep -i -E "pheno|genotype|\\.cnv" || echo "OK: no person-level file names"', check=False)
sh(f'du -sh {E}; find {E} -name null_runs.npz | wc -l')
""")

md("## 10. Package for UCSD (B4; dataset = `DS`)\nAsk support@researchallofus.org for the current egress size limit **before** downloading. Then download the two files from the Jupyter file browser (the dataset's run folder, e.g. `~/cnv_anthro_run/`) and `scp` them to the shared results folder at UCSD.")
code(r"""
T = f'{DATASETS[DS]["tar"]}_export_{time.strftime("%Y%m%d")}.tar.gz'
assert not bad, 'fix the export audit first'
sh(f'cd {DRUN} && tar czf {T} export_aggregate && sha256sum {T} > {T}.sha256 && ls -la {T}* && cat {T}.sha256')
sh(f'gsutil cp {DRUN}/{T} {DRUN}/{T}.sha256 {BUCKET_GS}/')
""")

md("""## 11. Array inputs from Molly's workspace ("PGC CNV project") — look only, nothing is copied yet
Her workspace shows a migrated Researcher-Workbench bucket (`rw-migration-aou-rw-<id>`) and the **C2022Q4R13 (v7)** CDR. Find the array CNV calls, the QC'd sample list and the batch file, then fill the three paths in section 12 (next notebook version).""")
code(r"""
sh('wb workspace list 2>&1 | head -40', check=False)        # note the ID of "PGC CNV project"
""")
code(r"""
MOLLY_WS = 'FILL: workspace id of "PGC CNV project" from the list above'
sh(f'wb resource list --workspace={MOLLY_WS} 2>&1', check=False)
sh(f'wb resource describe --workspace={MOLLY_WS} --name=<migrated bucket resource name> 2>&1', check=False)   # gives the gs:// bucket name
""")
code(r"""
MOLLY_GS = 'gs://FILL-bucket-name-from-describe'
sh(f'gsutil -u {os.environ["GOOGLE_PROJECT"]} ls {MOLLY_GS}/ 2>&1 | head -50', check=False)
# then drill down; file names and sizes only (no head/cat of participant-level files into chat):
# sh(f'gsutil -u {os.environ["GOOGLE_PROJECT"]} ls -lr {MOLLY_GS}/<folder>/ | grep -i -E "cnv|penncnv|quantisnp|sample|batch|qc" | head -50')
""")

md("""## 12. Build the three array inputs from Molly's PGC CNV calls (B5)
Source: `results/4_filtered/GDA-8v1-0/*.plt` (PGC CNV pipeline, consensus calls, 448 batches), sample QC `Freq_Calc == 1` in `results/QC/GDA-8v1-0/GDA-8v1-0.txt`, batch from `results/sample_info/GDA-8v1-0.txt`. Breakpoints = `consStart/consEnd`; only consensus calls (`algorithm` = IP,UP: PennCNV (UP) + QuantiSNP (IP in this project) agree, per Jeff); `4_filtered` = locus merge + sample QC only, so all size/probe/frequency filtering happens in `filter_plink_cnv.py`; TYPE 1 = DEL, 3 = DUP (DELDUP = caller disagreement, dropped); SITES = probes from `probes.txt` inside the call. Everything is written to the VM and mirrored to our own bucket (`array_inputs/`), nothing is printed per person.""")
code(r"""
MOLLY_GS = os.environ.get('MOLLY_GS', 'gs://FILL-bucket-of-the-PGC-CNV-workspace')   # set MOLLY_GS (bucket of Molly's PGC CNV workspace)
P = os.environ['GOOGLE_PROJECT']
ARUN = os.path.expanduser('~/cnv_anthro_array'); AIN = f'{ARUN}/inputs'; RAW = f'{ARUN}/molly_raw'
for d in [AIN, f'{RAW}/4_filtered', f'{ARUN}/data', f'{ARUN}/results']: os.makedirs(d, exist_ok=True)
if not glob.glob(f'{RAW}/4_filtered/*.plt'):
    sh(f'gsutil -u {P} -m -q cp {MOLLY_GS}/results/4_filtered/GDA-8v1-0/*.plt {RAW}/4_filtered/')
for f in ['results/QC/GDA-8v1-0/GDA-8v1-0.txt', 'results/sample_info/GDA-8v1-0.txt', 'inputs/probes.txt']:
    out = f'{RAW}/' + f.replace('/', '__')
    if not os.path.exists(out): sh(f'gsutil -u {P} -q cp {MOLLY_GS}/{f} {out}')
print(len(glob.glob(f'{RAW}/4_filtered/*.plt')), 'call files'); sh(f'du -sh {RAW}')
""")
code(r"""
import numpy as np
if 'mask' not in globals():   # defined in section 8; needed here when section 8 was not run in this kernel
    def mask(o):
        if isinstance(o, dict): return {k: mask(v) for k, v in o.items()}
        if isinstance(o, bool): return o
        if isinstance(o, (int, float)) and float(o).is_integer() and 1 <= o <= 20: return '<=20'
        return o
norm = lambda s: pd.to_numeric(s.astype(str).str.replace('chr', '', regex=False).replace({'X': '23', 'Y': '24'}), errors='coerce')

# samples passing QC + batch
qc = pd.read_csv(f'{RAW}/results__QC__GDA-8v1-0__GDA-8v1-0.txt', sep='\t', usecols=['identifier', 'Freq_Calc'], dtype=str)
si = pd.read_csv(f'{RAW}/results__sample_info__GDA-8v1-0.txt', sep='\t', usecols=['#batch', 'identifier'], dtype=str)
passing = qc.loc[qc.Freq_Calc == '1', 'identifier']
passing.to_csv(f'{AIN}/array_samples_qc_pass.txt', index=False, header=False)
si[['identifier', '#batch']].to_csv(f'{AIN}/array_batch.tsv', sep='\t', index=False, header=False)

# calls
c = pd.concat([pd.read_csv(f, sep='\t', usecols=['chrom', 'consStart', 'consEnd', 'consCallType', 'algorithm', 'sampleId'], dtype={'sampleId': str})
               for f in sorted(glob.glob(f'{RAW}/4_filtered/*.plt'))], ignore_index=True)
n_raw = len(c)
# algorithm x call type (aggregate; masked)
print(mask(pd.crosstab(c.algorithm, c.consCallType).to_dict()))  # no margins: a total minus the shown cells would reveal a masked one
# keep IP,UP consensus only (both callers), whatever the separator/order
both = c.algorithm.astype(str).str.upper().str.split(r'[,;+ ]+').apply(lambda a: {'IP', 'UP'} <= set(a))
n_nonconsensus = int((~both).sum()); c = c[both]
c['CHR'] = norm(c.chrom); c['TYPE'] = c.consCallType.map({'DEL': 1, 'DUP': 3})
c = c.dropna(subset=['CHR', 'TYPE']).astype({'CHR': int, 'TYPE': int, 'consStart': int, 'consEnd': int})

# SITES = probes inside [BP1, BP2]
pr = pd.read_csv(f'{RAW}/inputs__probes.txt', sep='\t', dtype={'Chromosome': str})
pr['CHR'] = norm(pr.Chromosome); pr = pr.dropna(subset=['CHR'])
c['SITES'] = 0
for ch, g in pr.groupby('CHR'):
    pos = np.sort(g.Position.astype(int).values); m = c.CHR == int(ch)
    c.loc[m, 'SITES'] = np.searchsorted(pos, c.loc[m, 'consEnd'].values, 'right') - np.searchsorted(pos, c.loc[m, 'consStart'].values, 'left')

out = pd.DataFrame({'FID': c.sampleId, 'IID': c.sampleId, 'CHR': c.CHR, 'BP1': c.consStart, 'BP2': c.consEnd,
                    'TYPE': c.TYPE, 'SCORE': 0, 'SITES': c.SITES})
out.to_csv(f'{AIN}/aou_array_pgc.cnv', sep='\t', index=False)

# aggregate checks only
print('calls read', n_raw, '| non-IP,UP dropped', mask(n_nonconsensus), '| DEL/DUP kept', len(out),
      '| DELDUP or bad chrom dropped', mask(n_raw - n_nonconsensus - len(out)))
print('samples: QC-pass', len(passing), '| with batch', si.identifier.nunique(), '| batches', si['#batch'].nunique())
print('call IDs not in sample_info:', mask(int((~out.IID.isin(si.identifier)).sum())),
      '| calls from QC-failed samples:', mask(int((~out.IID.isin(passing)).sum())))
print('SITES quantiles (5/50/95%):', out.SITES.quantile([.05, .5, .95]).tolist(), '| share of calls with 0 probes:', round((out.SITES == 0).mean(), 4))
print('share >= 10 kb:', round(((out.BP2 - out.BP1) >= 10000).mean(), 3), '| share >= 10 kb and >= 10 probes:',
      round((((out.BP2 - out.BP1) >= 10000) & (out.SITES >= 10)).mean(), 3))
open(f'{AIN}/BUILD_INFO.txt', 'w').write(f'{time.strftime("%F %T")}\tIP,UP consensus only; DEL/DUP; consStart/consEnd; SITES from probes.txt\n')
sh(f'ls -la {AIN}; gsutil -m -q rsync -r {AIN} {BUCKET_GS}/array_inputs && echo "mirrored to {BUCKET_GS}/array_inputs"')
""")
code(r"""
# preview of Jon's B5 calls-per-person line (same script and filters as the array prep step; temp output only)
sh(f'cd {CODE} && {sys.executable} filter_plink_cnv.py --cnv {AIN}/aou_array_pgc.cnv --samples {AIN}/array_samples_qc_pass.txt '
   f'--exclude_samples {RUN}/data/sv_samples.txt --out {ARUN}/data/preview.cnv --write_samples {ARUN}/data/preview_samples.txt '
   f'--min_size 10000 --min_probes 10 --max_freq 0.10 2>&1 | tail -3')
""")
code(r"""
# 12b. calls-per-person profile under stricter filters (aggregate only: medians and percentages; temp output deleted).
# GDA-8v1-0 has ~1.9M probes vs ~0.82M on the UKBB Axiom array, so 10 probes on GDA span far less DNA than 10 probes on UKBB;
# ~25 GDA probes is roughly the UKBB 10-probe resolution. UKBB reference: median 1 call/person, ~21% with none.
for size, probes, freq in [(10000, 10, 0.10), (10000, 25, 0.10), (20000, 20, 0.10), (50000, 25, 0.10), (10000, 10, 0.01), (20000, 25, 0.01)]:
    print(f'--- >= {size // 1000} kb, >= {probes} probes, freq <= {freq:.0%}', flush=True)
    sh(f'cd {CODE} && {sys.executable} filter_plink_cnv.py --cnv {AIN}/aou_array_pgc.cnv --samples {AIN}/array_samples_qc_pass.txt '
       f'--exclude_samples {RUN}/data/sv_samples.txt --out {ARUN}/data/diag.cnv --write_samples {ARUN}/data/diag_samples.txt '
       f'--min_size {size} --min_probes {probes} --max_freq {freq} 2>&1 | tail -1 | grep -o "calls/person.*"')
for f in ['diag.cnv', 'diag_samples.txt']:
    if os.path.exists(f'{ARUN}/data/{f}'): os.remove(f'{ARUN}/data/{f}')
""")

md("""## 13. WGS without the segdup filter (B6)
Same CNV and phenotype tables as the WGS run (no prep of its own); only `gwas.segdups = null`. Writes to `~/cnv_anthro_nosegdup/`, backup `run_backup_nosegdup/`. Runs on driver `'ns'`, so it can run alongside the array run (watch memory in `status()`; lower `GRID_JOBS` in section 1 if both grids run at once).""")
code(r"""
cfg = json.load(open(f'{CODE}/configs/aou_sacks2026_nosegdup.json'))
assert cfg['prep_commands'] == [] and cfg['gwas']['segdups'] is None
cfg['min_cell'] = MIN_CELL
json.dump(cfg, open(f"{CODE}/{DATASETS['nosegdup']['cfg']}", 'w'), indent=1)
shutil.copy(f"{CODE}/{DATASETS['nosegdup']['cfg']}", f"{DATASETS['nosegdup']['run']}/logs/config_used.json")
assert os.path.exists(os.path.expanduser(cfg['cnv_table'].replace('${HOME}', '~'))) and os.path.exists(os.path.expanduser(cfg['pheno_table'].replace('${HOME}', '~'))), 'WGS sv/prep tables missing'
print('nosegdup config ->', DATASETS['nosegdup']['cfg'], '| results', cfg['results_dir'], '| min_cell', cfg['min_cell'])
""")
code(r"""
# all strata: gwas, plots, alleles, loci, burden; then additivity/scores and grid/joint (stratum all). Export after QC (section 8 with DS = 'nosegdup').
launch(['main', 'additivity', 'gridjoint'], tag='ns', ds='nosegdup')
""")

md("""## 14. AoU array dataset (B5)
Inputs from section 12 (rerun it first: IP,UP consensus calls only). The prep step runs `filter_plink_cnv.py` (drops WGS participants; >= 10 kb, >= 10 probes, PGC 10% rule) and the phenotype extractor (BigQuery) with genotyping batch as a covariate. Writes to `~/cnv_anthro_array/`, backup `run_backup_array/`, driver `'ar'`. Jon 2026-11-20 (cell 12b: almost everyone has a call on this array): the 'only' runs use `no_rare_gene` as the primary reference and `nocall` as the sensitivity (`_nocall`); the burden panel uses `rare_gene` (`burden/`) plus an `all_calls` version (`burden_all_calls/`, its own stage).""")
code(r"""
AIN = os.path.expanduser('~/cnv_anthro_array/inputs')
assert os.path.exists(f'{AIN}/BUILD_INFO.txt') and 'IP,UP' in open(f'{AIN}/BUILD_INFO.txt').read(), 'rerun section 12 (consensus filter) first'
cfg = json.load(open(f'{CODE}/configs/aou_sacks2026_array.json'))
f0, f1 = cfg['prep_commands']
assert f0[0] == 'filter_plink_cnv.py' and f1[0] == 'aou_pheno_extract.py'
f0[f0.index('--cnv') + 1] = f'{AIN}/aou_array_pgc.cnv'
f0[f0.index('--samples') + 1] = f'{AIN}/array_samples_qc_pass.txt'
f1[f1.index('--ancestry') + 1] = f'{RUN}/data/ancestry_preds.tsv'
f1[f1.index('--batch') + 1] = f'{AIN}/array_batch.tsv'
cfg['min_cell'] = MIN_CELL
# Jon 2026-11-20: keep the B5 call filters; median 8 calls/person and 0.1% with no call (UKBB 1 and ~21%), so
# - 'only' runs: primary reference = no rare gene-hitting call (WGS rule); 'nocall' = sensitivity (runs suffixed _nocall)
# - burden panel: primary single_def = rare_gene (-> burden/); all_calls = second version (stage burden_allcalls -> burden_all_calls/)
cfg['figures']['only_ref'] = 'no_rare_gene'
cfg['figures']['only_ref_sensitivity'] = ['nocall']
cfg['figures'].setdefault('rare_freq', 0.01)
cfg['burden']['single_def'] = 'rare_gene'
cfg['burden'].setdefault('rare_freq', 0.01)
# second burden version: same dup_burden.py call as run_pipeline.py (same pheno, covariates, strata), only --single_def all_calls
open(f"{DATASETS['array']['run']}/logs/burden_allcalls.py", 'w').write('''import json, os, subprocess, sys
sys.path.insert(0, os.getcwd())
from run_pipeline import expand
C = expand(json.load(open(sys.argv[1]))); B = C['burden']; root = C['results_dir']; n = 0
for sname, filt in C['strata'].items():
    out = os.path.join(root, sname); pheno = os.path.join(out, 'stratum_pheno.tsv') if filt else C['pheno_table']
    for t, T in C['traits'].items():
        if not os.path.isdir(f'{out}/burden/{t}'):
            print('skip (no primary burden run):', sname, t, flush=True); continue
        fac = ','.join(f for f in T['factors'] if not (f == 'sex' and 'sex' in (filt or {})))
        cmd = [sys.executable, 'dup_burden.py', '--pheno', pheno, '--trait', T['column'], '--cnv', C['cnv_table'], '--exons', B['exons'],
               '--gene_cov', B['gene_cov'], '--clingen', B['clingen'], '--asym_loci', B['asym_loci'], '--numeric', ','.join(T['numeric']),
               '--factors', fac, '--out', f'{out}/burden_all_calls/{t}', '--min_n', str(B.get('min_n', 50)), '--single_def', 'all_calls',
               '--rare_freq', str(B.get('rare_freq', 0.01))]
        print('+', ' '.join(cmd), flush=True); subprocess.run(cmd, check=True); n += 1
print('burden_allcalls runs:', n)
''')
cfg['note'] += " Jon 2026-11-20: 'only'-run reference = no_rare_gene (nocall as sensitivity); burden single_def = rare_gene (+ all_calls in burden_all_calls/)."
assert 'FILL' not in json.dumps(cfg['prep_commands']), [x for x in json.dumps(cfg['prep_commands']).split('"') if 'FILL' in x]   # the note mentions FILL_, so check only the paths
for p_ in [f'{AIN}/aou_array_pgc.cnv', f'{AIN}/array_samples_qc_pass.txt', f'{AIN}/array_batch.tsv', f'{RUN}/data/ancestry_preds.tsv', f'{RUN}/data/sv_samples.txt']:
    assert os.path.exists(p_), p_
json.dump(cfg, open(f"{CODE}/{DATASETS['array']['cfg']}", 'w'), indent=1)
shutil.copy(f"{CODE}/{DATASETS['array']['cfg']}", f"{DATASETS['array']['run']}/logs/config_used.json")
print('array config ->', DATASETS['array']['cfg'], '| min_cell', cfg['min_cell'], '| only_ref', cfg['figures']['only_ref'],
      '+', cfg['figures']['only_ref_sensitivity'], '| burden', cfg['burden']['single_def'], '+ all_calls (burden_all_calls/)')
print('prep 1:', ' '.join(f0)); print('prep 2:', ' '.join(f1))
""")
code(r"""
# prep -> all strata -> additivity/scores -> grid/joint. Export after QC (section 8 with DS = 'array').
launch(['prep', 'main', 'burden_allcalls', 'additivity', 'gridjoint'], tag='ar', ds='array')
""")
code(r"""
# after the array prep stage: calls per person (send Jon this line) and the phenotype summary (both aggregate)
A = DATASETS['array']['run']
sh(f"grep -h -i -E 'calls|person|samples|none' {A}/logs/prep.log | grep -v -E '^\\s*$' | tail -8", check=False)
p_ = f'{A}/results/aou_pheno_summary.json'
if os.path.exists(p_): print(json.dumps(mask(json.load(open(p_))), indent=1) if 'mask' in globals() else 'run section 8 first for mask()')
""")


# ---- section 15: Molly's recurrent CNV genotyping (port), embedded from aou_workbench/ at build time ----
import os as _os; _WB = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..", "recurrentCNVgenotyping", "aou")
_SRC = {f: open(_os.path.join(_WB, f)).read() for f in ['recurrent_cnv_genotyping_aou.py', 'recurrentCNV_hg38.reciprocal.tsv', 'recurrentCNV_hg38.oneway.tsv']}
md("""## 15. Recurrent CNV genotyping with Molly's method (Jon 2026-11-23)
Port of `recurrentCNVgenotyping/GenotypeRecurrentCNVs_tutorial.ipynb` (same locus sheet, same 50% reciprocal / one-way rules, same segment patterns; no merging of fragmented calls). Runs on the WGS and array call tables and the people analysed. The genotype matrices stay in each run's `data/` folder (participant level, never exported). The cell prints masked DEL/DUP counts per locus next to the pipeline's `alleles.py` counts (rounded to 10). Safe to paste.""")
code("""
# write the port and Molly's two locus tables (saved from her public sheet; the Workbench has no Google access)
MG = os.path.expanduser('~/recurrent_cnv_genotyping'); os.makedirs(MG, exist_ok=True)
for f, txt in """ + repr(_SRC) + """.items():
    open(f'{MG}/{f}', 'w').write(txt)
print(sorted(os.listdir(MG)))
""")
code(r"""
r10 = lambda n: f'~{int(round(n, -1))}' if n > MIN_CELL else mask(int(n))
norm = lambda x: re.sub(r'[^a-z0-9]', '', str(x).lower().replace('x_', '', 1))
import re
res = {}
for ds, cnv_, ph_ in [('wgs', f"{DATASETS['wgs']['run']}/data/aou_srwgs_sv.cnv", f"{DATASETS['wgs']['run']}/data/aou_pheno_std.tsv"),
                      ('array', f"{DATASETS['array']['run']}/data/aou_array.cnv", f"{DATASETS['array']['run']}/data/aou_array_pheno_std.tsv")]:
    out_ = f"{DATASETS[ds]['run']}/data/recurrent_genomatrix_molly.tsv"
    sh(f'{sys.executable} {MG}/recurrent_cnv_genotyping_aou.py --cnv {cnv_} --samples {ph_} --reciprocal {MG}/recurrentCNV_hg38.reciprocal.tsv '
       f'--oneway {MG}/recurrentCNV_hg38.oneway.tsv --out {out_} --min_cell {MIN_CELL}')
    G = pd.read_csv(out_, sep='\t', index_col=0)
    G = G[[c for c in G.columns if not c.endswith('_oneway')]]
    al = pd.read_csv(f"{DATASETS[ds]['run']}/results/all/taskA_alleles/height/allele_assoc.tsv", sep='\t')
    al = al.assign(key=al.allele.str.replace(r'_(DEL|DUP)$', '', regex=True).map(norm), typ=al.allele.str[-3:])
    rows = []
    for c in G.columns:
        r = dict(locus=c)
        for typ, code_ in [('DEL', 1), ('DUP', 3)]:
            r[f'{typ}_molly'] = r10(int((G[c] == code_).sum()))
            m = al[(al.key == norm(c)) & (al.typ == typ)]
            r[f'{typ}_alleles_py'] = (r10(int(pd.to_numeric(m.n_carriers, errors='coerce').fillna(0).iloc[0])) if str(m.n_carriers.iloc[0]).isdigit() else m.n_carriers.iloc[0]) if len(m) else '-'
        rows.append(r)
    res[ds] = pd.DataFrame(rows)
    print(f'\n=== {ds}: N = {len(G)} (people analysed); counts rounded to 10, 1-{MIN_CELL} masked; alleles_py = pipeline (all/height)')
    print(res[ds].to_string(index=False))
""")

md(r"""## 16. Fig. 3 joint model and gene-based association (Jon 2026-11-23, `AllOfUs/HANDOFF_Milad_Fig3.md`)
WGS and array only (nosegdup reads the same CNV and phenotype tables as WGS, so it is not needed). Before this section:
set `PIN_COMMIT = 'bdfcf29'` (section 1, done), rerun sections **1, 2, 4 and 6**, rerun section 7 (STAGES now has `fig3`, `fig3_export`,
`genes`, `genes_export`), run section 8 once (for `mask()`), and rerun **only the config cell of section 14** (not its launch cell): both filled
configs then carry the template's new `"fig3"` and `"genes"` blocks. Do not edit `designs/fig3_spec_2026-11-23.json`.

Order: 16.1 check → 16.2 launch `fig3` then `genes` (WGS and array, two drivers) → 16.3 Fig. 3 QC (paste) → 16.4 gene QC (paste) →
16.5 locus-carrier check (eight lines, paste) → 16.6a audit fix → 16.6 exports (`AUDIT OK`) → 16.7 package. Four tarballs; download them at least an hour apart.""")
code(r"""
# 16.1 what changed in the pipeline since the finished stages, and the two filled configs
sh(f'cd {REPO} && git diff --stat 4b64136 bdfcf29 -- AllOfUs/cnv_anthro | tail -20')
for ds in ['wgs', 'array']:
    c_ = json.load(open(f"{CODE}/{DATASETS[ds]['cfg']}"))
    assert c_.get('fig3') and c_.get('genes'), f"{ds}: no fig3/genes block -> rerun section 6 (wgs) / the section 14 config cell (array)"
    shutil.copy(f"{CODE}/{DATASETS[ds]['cfg']}", f"{DATASETS[ds]['run']}/logs/config_used.json")
    spec_ = json.load(open(f"{CODE}/{c_['fig3']['spec']}"))
    print(ds, '| min_cell', c_['min_cell'], '| fig3 strata', c_['fig3']['strata'], '| spec', c_['fig3']['spec'], '| shifts sha1', spec_.get('perms_sha1'))
for t in ['height', 'bmi']:
    f_ = f'{ROOT}/resources/shared_shifts/{t}_gtex_nonbrain_all.npy'      # the config's ../resources, relative to cnv_anthro
    print(t, 'shift list sha1', hashlib.sha1(open(f_, 'rb').read()).hexdigest()[:16], '(UKBB: height 157a247259f02719, bmi be761a32dd575298)')
""")
code(r"""
# 16.2 launch. fig3 = several hours (WGS longer: ~57 calls per person); genes ~1 h. Two drivers ('' = wgs, 'ar' = array) run at once on the
# 251 GB VM; if status() shows memory pressure, stop one, change '--jobs 5' to '--jobs 3' in STAGES['fig3'], rerun section 7 and relaunch
# (finished fits, i.e. those with run_info.json, are skipped).
launch(['fig3', 'genes'], tag='', ds='wgs')
launch(['fig3', 'genes'], tag='ar', ds='array')
""")
code(r"""
# 16.3 Fig. 3 QC (handoff section B): z, P and N only -> safe to paste. Set DS, run once per dataset after its fig3 stage is done.
import numpy as np
DS = 'wgs'
assert 'mask' in globals(), 'run section 8 first (mask)'
F3 = f"{DATASETS[DS]['run']}/results/all/fig3"
print('=== Fig. 3 QC |', DS, '|', COMMIT)
fi = json.load(open(f'{F3}/cnv_filter/run_info.json'))        # call set: nothing is dropped; per-position carrier frequency recorded
print('call set: N', fi['N'], '| calls', fi['calls_in'], '(kept', fi['calls_out'], ') | highest per-position carrier frequency', fi['max_position_freq'])
print('what a position filter WOULD drop (not applied):', json.dumps(mask({t: {d: v['calls_dropped'] for d, v in r.items()} for t, r in fi['report'].items()})))
spec = json.load(open(f'{CODE}/designs/fig3_spec_2026-11-23.json'))
for f in sorted(glob.glob(f'{F3}/*/grid/*/*/run_info.json')):
    i = json.load(open(f)); t = i['trait'].replace('_rint', '')
    print(os.path.relpath(os.path.dirname(f), F3), '| N', i['N'], '| nrep', i['nrep'], '| shifts',
          'OK' if i['perms_sha1'] == spec['perms_sha1'][t] else 'MISMATCH', '| major loci', len(i['major_loci']))
    o = pd.read_csv(f.replace('run_info.json', 'observed_terms.tsv'), sep='\t')
    print('   (run, term) cells with 1-20 carriers, withheld at export:', int(pd.to_numeric(o.n_carriers_term, errors='coerce').between(1, 20).sum()),
          'of', int(o.n_carriers_term.notna().sum()))
for f in sorted(glob.glob(f'{F3}/uni_*/grid/*/univariate/grid_summary.tsv')):      # univariate DEL z / DUP z per class and CNV-size stratum
    d = pd.read_csv(f, sep='\t'); d = d[~d.readout.str.endswith('|mirror')]
    d['z'] = np.where(d.readout.str.contains('n_all'), d.t_analytic, d.z_emp)       # total genes: analytic t (its shift null is degenerate)
    print('---', os.path.relpath(os.path.dirname(f), F3)); print(d[['readout', 'z', 'P_emp', 'loo_held', 'loo_n', 'noloci_z_emp']].round(3).to_string(index=False))
for name in ['size2t_mirror', 'all2t_mirror', 'totdose']:
    for tab in ['mirror4_cross_summary.tsv', 'dose_summary.tsv']:
        p = f'{F3}/{name}/{tab}'
        if os.path.exists(p) and os.path.getsize(p) > 5:
            d = pd.read_csv(p, sep='\t'); d = d[d.stat.isin(['M4', 'armDEL', 'armDUP', 'D4', 'doseOth', 'doseTot', 'mirTot'])]
            print('---', name, tab); print(d[['trait', 'size', 'stratum', 'stat', 'z', 'P', 'loo_held', 'loo_n', 'noloci_z']].round(3).to_string(index=False))
for t in ['height', 'bmi']:
    d = pd.read_csv(f'{F3}/gene_mirror/{t}/gene_mirror_summary.tsv', sep='\t')
    print('--- gene mirror', t); print(d[['size', 'stratum', 'stat', 'obs', 'z', 'P', 'n_genes', 'loo_held', 'n_loo', 'noloci_z']].round(3).to_string(index=False))
for t in ['height', 'bmi']:      # dose-response bins: already masked when written (bins with 1-20 people have no n and no estimate)
    i = json.load(open(f'{F3}/dose_bins/{t}/run_info.json')); d = pd.read_csv(f'{F3}/dose_bins/{t}/fig3_dose_bins.tsv', sep='\t')
    d = d[(d.run == 'full') & (d.side != 'none')]; est = d.beta.notna()
    print('--- dose bins', t, '| N', i['N'], '| single-dosage or non-carriers', i['full_N_single'], '| bins with an estimate', int(est.sum()), 'of', len(d),
          '| withheld', int((d.n.astype(str) == '<=20').sum()), '| empty', int((d.n.astype(str) == '0').sum()))
    print(d[est & (d.partition == 'uni')].groupby(['size', 'cls', 'side']).size().unstack('side').to_string())
""")
code(r"""
# 16.4 gene-based association QC (handoff section E): N, numbers of genes and test statistics only -> safe to paste. Never open gene_assoc.tsv.gz itself.
DS = 'wgs'
R = f"{DATASETS[DS]['run']}/results"
print('=== gene QC |', DS)
for f in sorted(glob.glob(f'{R}/*/gene_assoc*/*/gene_assoc.tsv.gz')):
    d = pd.read_csv(f, sep='\t'); st, t = f.split('/')[-4] + ' ' + f.split('/')[-3], f.split('/')[-2]; d = d[d.n_carriers >= 21]
    for dos in ['DEL', 'DUP']:
        x = d[d.dosage == dos]
        print(st, t, dos, '| N', int(x.N.iloc[0]) if len(x) else 0, '| genes with >= 21 carriers', len(x),
              '| lambda', round(float((x.t ** 2).median() / 0.4549), 2) if len(x) else None, '| genes with P < 2.5e-6:', int((x.P < 2.5e-6).sum()))
""")
code(r"""
# 16.5 locus carriers at 1q21.1 and 22q11.21 (handoff section F): people, 1-20 masked -> the eight lines for Jon. Do not change the locus rule.
sys.path.insert(0, CODE)
from common import read_plink_cnv
from conga_axes_qt import assign_loci
X = lambda v: os.path.expandvars(v)
loci = pd.read_csv(f'{CODE}/../resources/growth_loci/loci_sacks34.tsv', sep='\t')
for DS in ['wgs', 'array']:
    cfg_ = json.load(open(f"{DATASETS[DS]['run']}/logs/config_used.json"))
    ids = pd.read_csv(X(cfg_['pheno_table']), sep='\t', usecols=['IID'], dtype=str).IID
    c = read_plink_cnv(X(cfg_['cnv_table']), keep_iids=ids); c = c[c.CHR <= 22].reset_index(drop=True)
    c['locus'] = assign_loci(c, loci)
    for name, ch, a, b, fam in [('1q21.1 distal', 1, 147056425, 147922330, ['1q21_1']),
                                ('22q11.21 A-D', 22, 19049819, 21112437, ['22q11_21_A_D', '22q11_21_A_B', '22q11_21_B_D', '22q11_21_B_C', '22q11_21_C_D'])]:
        for typ, tn in [(1, 'DEL'), (3, 'DUP')]:
            d = c[(c.CHR == ch) & (c.TYPE == typ) & (c.BP2 > a) & (c.BP1 < b)].copy()
            d['cov'] = (d.BP2.clip(upper=b) - d.BP1.clip(lower=a)) / (b - a)
            assigned = set(d[d.locus.isin(fam)].IID); cover50 = set(d[d['cov'] >= .5].IID); any10 = set(d[d['cov'] >= .1].IID)
            print(DS, '|', name, tn, '| assigned by the 50% reciprocal rule:', mask(len(assigned)), '| one call covers >= 50% of the locus:', mask(len(cover50)),
                  '| of those not assigned:', mask(len(cover50 - assigned)), '| any call covering >= 10%, not assigned:', mask(len(any10 - assigned)))
""")
code(r"""
# 16.6a fix the export audit of run_fig3.py at bdfcf29 (2026-11-24): in audit(), the fig3_dose_bins.tsv branch assigned the table to `d`,
# the os.walk folder variable, so the next file of dose_bins/<trait>/ crashed (TypeError: ... not DataFrame) after the export was written
# but before the audit and EXPORT_MANIFEST.txt. Rename only; no rule changes. No-op once fixed upstream. Rerunning section 2 stashes it.
p_ = f'{CODE}/run_fig3.py'; s_ = open(p_).read()
old_ = ('d = pd.read_csv(p, sep="\t"); n = pd.to_numeric(d.n, errors="coerce")\n'
        '                if n.between(1, k).any() or (d.run.ne("full") & n.notna() & d.side.ne("none")).any():')
new_ = ('db_ = pd.read_csv(p, sep="\t"); n = pd.to_numeric(db_.n, errors="coerce")\n'
        '                if n.between(1, k).any() or (db_.run.ne("full") & n.notna() & db_.side.ne("none")).any():')
if old_ in s_:
    open(p_, 'w').write(s_.replace(old_, new_)); print('patched run_fig3.py audit()')
else:
    print('nothing to patch' if new_ in s_ or 'd = pd.read_csv(p, sep="\t"); n = pd.to_numeric(d.n' not in s_ else 'PATTERN CHANGED: check run_fig3.py audit()')
sh(f'cd {REPO} && git diff --stat', check=False)
# an export that crashed before its audit leaves an unaudited folder: move it aside (never download it)
for ds in ['wgs', 'array']:
    D = DATASETS[ds]['run']
    if os.path.exists(f'{D}/export_aggregate_fig3') and not glob.glob(f'{D}/export_aggregate_fig3/*/fig3/EXPORT_MANIFEST.txt'):
        sh(f'mv {D}/export_aggregate_fig3 {D}/export_aggregate_fig3_noaudit_{time.strftime("%Y%m%d_%H%M")}')
""")
code(r"""
# 16.6 exports, after 16.3 and 16.4 look fine. Each export refuses to overwrite an existing folder (move it away to redo).
launch(['fig3_export', 'genes_export'], tag='', ds='wgs')
launch(['fig3_export', 'genes_export'], tag='ar', ds='array')
""")
code(r"""
# 16.7 audit lines and packaging (set DS). Then download each .tar.gz + .sha256 (one tarball per hour at most) and scp to
# the shared results folder at UCSD
DS = 'wgs'
DRUN = DATASETS[DS]['run']
for kind, man in [('fig3', 'export_aggregate_fig3/all/fig3/EXPORT_MANIFEST.txt'), ('genes', 'export_aggregate_genes/EXPORT_MANIFEST.txt')]:
    last = open(f'{DRUN}/logs/{kind}_export.log').read().strip().splitlines()[-3:]
    assert any('AUDIT OK' in l for l in last), f'{kind}_export.log does not end with AUDIT OK: {last}'
    T = f"{DATASETS[DS]['tar']}_{kind}_export_{time.strftime('%Y%m%d')}.tar.gz"
    sh(f"head -3 {DRUN}/{man}")
    sh(f'cd {DRUN} && tar czf {T} export_aggregate_{kind} && sha256sum {T} > {T}.sha256 && ls -la {T}* && gsutil cp {T} {T}.sha256 {BUCKET_GS}/')
""")

md(r"""## 17. Rerun with one recurrent-CNV genotyper (Jon 2026-11-24, commit `2e638c1`, top of `AllOfUs/HANDOFF_Milad_Fig3.md`)
Every step that needs recurrent-CNV carriers (locus-removal runs, locus effects, interval strata) now takes them from Molly's genotyping
(`recurrentCNVgenotyping/aou`, imported by `common.py`; keep the whole repo checkout). Not affected: probe GWAS, gene-based association,
alleles, frequency filters. The fit steps skip finished folders, so everything goes to a new results folder inside each run folder
(`<run>/rcg_20261124/`); the CNV and phenotype tables are reused as they are.

Before this section: `PIN_COMMIT = '2e638c1'` (section 1, done), rerun sections **1, 2, 4, 6 and 7**, section 8 once (for `mask()`) and the
**config cell of section 14** (not its launch cell). Section 2 stashes the 16.6a patch (it is upstream now).

Rerunning section 7 resets `DATASETS`: run 17.1 again after it (safe to repeat).

Order: 17.1 configs → 17.2 launch (`rcg`, `rcg_gridjoint`, `fig3`, arrays also `burden_allcalls`) → 16.3 QC with `DS = 'wgs_rcg'` / `'array_rcg'` (paste)
→ 17.3 exports (`export` + `fig3_export`) → 17.4 audit and package. Four tarballs; download them at least an hour apart.""")
code(r"""
# 17.1 new configs: same filled configs, only results_dir / export_dir moved to <run>/rcg_20261124
assert COMMIT.startswith('2e638c1'), f'code at {COMMIT}: pin 2e638c1 and rerun sections 1, 2, 4, 6, 7 and the section 14 config cell'
sh(f'cd {REPO} && git status --short | head; git stash list | head -3', check=False)
sh(f'cd {CODE} && {sys.executable} -c "import common; print(\'genotyper:\', common.rcg.__file__)"')
RCG = 'rcg_20261124'
for ds in ['wgs', 'array']:
    D = DATASETS[ds]; c_ = json.load(open(f"{CODE}/{D['cfg']}"))
    assert c_.get('fig3'), f"{ds}: filled config has no fig3 block -> rerun section 6 (wgs) / the section 14 config cell (array)"
    run_ = f"{D['run']}/{RCG}"
    for d in ['results', 'logs', 'stages']: os.makedirs(f'{run_}/{d}', exist_ok=True)
    c_['results_dir'], c_['export_dir'] = f'{run_}/results', f'{run_}/export_aggregate'
    cfg_ = D['cfg'].replace('_filled.json', f'_{RCG}_filled.json')
    json.dump(c_, open(f'{CODE}/{cfg_}', 'w'), indent=1); shutil.copy(f'{CODE}/{cfg_}', f'{run_}/logs/config_used.json')
    if ds == 'array':                         # second burden version (all_calls), same script as section 14
        shutil.copy(f"{D['run']}/logs/burden_allcalls.py", f'{run_}/logs/burden_allcalls.py')
    DATASETS[f'{ds}_rcg'] = dict(run=run_, cfg=cfg_, backup=f"{D['backup']}/{RCG}", tar=f"{D['tar']}_rcg",
                                 stages=['rcg', 'rcg_gridjoint', 'fig3'] + (['burden_allcalls'] if ds == 'array' else []) + ['export', 'fig3_export'])
    print(ds, '->', cfg_, '| results', c_['results_dir'], '| cnv', c_['cnv_table'], '| min_cell', c_['min_cell'])
""")
code(r"""
# 17.2 launch both datasets (two drivers). rcg = loci + burden on 9 strata and additivity on 'all'; rcg_gridjoint = grid + joint on 'all'
# (several hours, like gridjoint before); fig3 = run_fig3.py fit,summary (several hours). If status() shows memory pressure, stop one driver and relaunch later
# (finished stages and finished fits are skipped).
launch(['rcg', 'rcg_gridjoint', 'fig3'], tag='', ds='wgs_rcg')
launch(['rcg', 'rcg_gridjoint', 'fig3', 'burden_allcalls'], tag='ar', ds='array_rcg')
""")
code(r"""
status()      # the two new datasets show up as wgs_rcg and array_rcg
""")
code(r"""
# 17.3 exports, after 16.3 (DS = 'wgs_rcg' / 'array_rcg') looks fine. export_aggregate.py copies everything under results/ and has no rule
# for fig3/ (the first run exported before fig3 existed; here fig3 ran first and the export failed inside results/all/fig3). Fig. 3 has its
# own audited export, so results/all/fig3 is held aside while the main export runs (< 1 min, in this cell) and is always put back.
# Then fig3_export runs in the background as before (log ends with AUDIT OK). Each export refuses to overwrite an existing folder.
def export_rcg(ds, tag):
    D = DATASETS[ds]; RUN_ = D['run']; A = f'{RUN_}/results/all'
    need = [s for s in D['stages'] if s not in ('export', 'fig3_export')]
    assert all(os.path.exists(f'{RUN_}/stages/{s}.done') for s in need), f'{ds}: wait for {need}'
    assert not running(tag), f'driver{tag} still running'
    assert not os.path.exists(f'{RUN_}/export_aggregate'), f'{RUN_}/export_aggregate exists: move it away first'
    if not os.path.exists(f'{RUN_}/stages/export.done'):
        os.rename(f'{A}/fig3', f'{RUN_}/fig3_hold')
        try:
            rc = sh(f"cd {CODE} && {sys.executable} run_pipeline.py --config {D['cfg']} --steps export > {RUN_}/logs/export.log 2>&1", check=False)
        finally:
            os.rename(f'{RUN_}/fig3_hold', f'{A}/fig3')
        if rc:
            x = open(f'{RUN_}/logs/export.log', errors='replace').read().splitlines()
            k = next((n for n, l in enumerate(x) if l.startswith('Traceback')), 0)     # first traceback = export_aggregate.py's own
            print('\n'.join(l[:250] for l in x[k:k + 30] if '\t' not in l)); raise RuntimeError(f'{ds}: export failed')
        if os.path.exists(f'{RUN_}/stages/export.failed'): os.remove(f'{RUN_}/stages/export.failed')
        open(f'{RUN_}/stages/export.done', 'w').write(COMMIT + '\n')
        sh(f'ls {RUN_}/export_aggregate; du -sh {RUN_}/export_aggregate; find {RUN_}/export_aggregate -path "*fig3*" | head -2', check=False)
    launch(['fig3_export'], tag=tag, ds=ds)

export_rcg('wgs_rcg', '')
# export_rcg('array_rcg', 'ar')      # after array_rcg fig3 and burden_allcalls are DONE
""")
code(r"""
# 17.3b small grid cells in the main export. Since 965703c conga_grid_qt.py writes n_carriers_term into observed_terms.tsv (grid_univariate,
# grid_callsize, joint); export_aggregate.py's COUNT_COLS does not list it, so cells with 1-20 carriers went out unmasked (the first run's grid,
# at 05fa65b, had no such column). Same rule as run_fig3.export_fit, applied to the exported copy only: n_carriers_term 1..20 -> '<=20',
# beta/se/t of that (run, readout) and of its |mirror readout -> blank. Shifted-genome nulls (null_runs.npz) are left as they are.
def mask_grid_terms(ds, k=MIN_CELL):
    E_ = f"{DATASETS[ds]['run']}/export_aggregate"; nf = nc = 0
    for f in glob.glob(f'{E_}/**/observed_terms.tsv', recursive=True):
        ob = pd.read_csv(f, sep='\t', dtype={'n_carriers_term': object})
        if 'n_carriers_term' not in ob.columns: continue
        n = pd.to_numeric(ob.n_carriers_term, errors='coerce'); sm = n.between(1, k)
        if not sm.any(): continue
        small = {(r, q) for r, q in zip(ob.run[sm], ob.readout[sm])}
        small |= {(r, q.rsplit('|', 1)[0] + '|mirror') for r, q in list(small)}
        hit = np.array([(r, q) in small for r, q in zip(ob.run, ob.readout)])
        for c in [c for c in ['beta', 'se', 't'] if c in ob.columns]: ob.loc[hit, c] = np.nan
        ob.loc[sm.to_numpy(), 'n_carriers_term'] = f'<={k}'
        ob.to_csv(f, sep='\t', index=False, float_format='%.10g'); nf += 1; nc += int(sm.sum())
    open(f'{E_}/EXPORT_MANIFEST.txt', 'a').write(f'notebook 17.3b: observed_terms n_carriers_term 1..{k} masked, beta/se/t blanked: {nc} cells in {nf} tables\n')
    print(ds, '| tables masked', nf, '| cells', nc)
mask_grid_terms('wgs_rcg')
# mask_grid_terms('array_rcg')     # after export_rcg('array_rcg', 'ar')
""")
code(r"""
# 17.4 audit and packaging (set DS). Then download each .tar.gz + .sha256 (one tarball per hour at most) and scp to
# the shared results folder at UCSD
import re
DS = 'wgs_rcg'                          # 'wgs_rcg' | 'array_rcg'
DRUN = DATASETS[DS]['run']
last = open(f'{DRUN}/logs/fig3_export.log').read().strip().splitlines()[-3:]
assert any('AUDIT OK' in l for l in last), f'fig3_export.log does not end with AUDIT OK: {last}'
sh(f'head -3 {DRUN}/export_aggregate_fig3/all/fig3/EXPORT_MANIFEST.txt')
# same audit as 9.1 on the main export: no count from 1 to 20 in any count-like column
E = f'{DRUN}/export_aggregate'
COUNT = re.compile(r'^(n|N)$|^n_|_n$|^N_|carrier|count|^ncnv|^n[A-Z]', re.I)
NOT_COUNT = re.compile(r'^n_probes$|^noloci_|^(mean|sd)_count$', re.I)
bad, nfile = [], 0
for f in glob.glob(f'{E}/**/*.tsv', recursive=True) + glob.glob(f'{E}/**/*.csv', recursive=True):
    nfile += 1
    d = pd.read_csv(f, sep='\t' if f.endswith('.tsv') else ',', dtype=str, low_memory=False)
    for c in [c for c in d.columns if COUNT.search(str(c)) and not NOT_COUNT.search(str(c))]:
        k = int(pd.to_numeric(d[c], errors='coerce').between(1, 20).sum())
        if k: bad.append((os.path.relpath(f, E), c, k))
print(nfile, 'tables checked'); print('\n'.join(f'VIOLATION {f} [{c}] {k} rows' for f, c, k in bad) or 'OK: no 1-20 counts in count columns')
sh(f'find {E} -type f | grep -i -E "pheno|genotype|\\.cnv" || echo "OK: no person-level file names"', check=False)
assert not bad, 'fix the export audit first'
for kind, folder in [('export', 'export_aggregate'), ('fig3_export', 'export_aggregate_fig3')]:
    T = f"{DATASETS[DS]['tar']}_{kind}_{time.strftime('%Y%m%d')}.tar.gz"
    sh(f'cd {DRUN} && tar czf {T} {folder} && sha256sum {T} > {T}.sha256 && ls -la {T}* && gsutil cp {T} {T}.sha256 {BUCKET_GS}/')
""")

nb = {"cells": cells, "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                                   "language_info": {"name": "python"}}, "nbformat": 4, "nbformat_minor": 5}
json.dump(nb, open(sys.argv[1], "w"), indent=1)
print("wrote", sys.argv[1], len(cells), "cells")
