#!/usr/bin/env python
"""Model figure (Jon 2026-11-15): how the dose response was dissected. Univariate screen of TWAS and constraint measures, selection,
joint model, first in the 34 recurrent CNV intervals (Sacks et al.) and then genome-wide (UK Biobank), and the schematic of the
resulting architecture of height and BMI as the final panel (style of the dosage report's axis schematic: arrows in DEL x DUP space,
thick filled = held, thin open = not).
Inputs (aggregate tables only):
  locus grid      two_axis_CNV_dose_response/growth_newloci_2026-11-15/ref34/results/eight_traits/growth_conga/growth_grid_conga.csv
  locus joint     .../growth_conga_axes_summary.csv (ConGA fit, unsigned TWAS count), growth_conga_presence_vs_count.csv,
                  v2_2026-11-08/results/eight_traits/growth_conga/yaxis_monogenic.csv (duplication-axis single hit)
  genome-wide     UKBB_anthro_CNV/results_aggregate/growthloci_2026-11-15/grid_summary_all.tsv (univariate, with recurrent-loci removal)
                  UKBB_anthro_CNV/results_aggregate/joint_sacks34_2026-11-15/{height_height,bmi_bmi}_selected.tsv (joint)
Usage: python make_fig_models.py  -> Fig_models.{pdf,png}
AllOfUs copy (2026-11-16): GENOMEWIDE_ROOT / GENOMEWIDE_LABEL / FIG_OUT / SEBATLAB_DROPBOX as in make_fig2_fig3_reorg.py (panels E-G and the
single-hit tests from the meta-analysis when GENOMEWIDE_ROOT is the meta_sacks2026.py output).
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

OUTDIR = Path(os.environ.get('FIG_OUT', Path(__file__).resolve().parent))
JC = Path(os.environ.get('SEBATLAB_DROPBOX', Path.home() / 'sebatlab Dropbox')) / 'Jon/Claude'
LOC = JC / 'two_axis_CNV_dose_response/growth_newloci_2026-11-15/ref34/results/eight_traits/growth_conga'
V2 = JC / 'two_axis_CNV_dose_response/v2_2026-11-08/results/eight_traits/growth_conga'
RA = Path(os.environ.get('GENOMEWIDE_ROOT', JC / 'UKBB_anthro_CNV/results_aggregate'))
SRC = os.environ.get('GENOMEWIDE_LABEL', 'UK Biobank')
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial'], 'font.size': 8,
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold',
    'axes.titlesize': 8.4, 'axes.labelsize': 8, 'axes.linewidth': .7, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'legend.fontsize': 7, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})
INK = '#1a1a1a'; GREY = '#8a8a8a'; LIGHT = '#dcdcdc'
HCOL, BCOL = INK, '#c47b1a'
SIGC, UNSC, CONC = '#9a9a9a', '#e08214', '#1b7837'          # signed TWAS, unsigned TWAS, constraint
MIRC, HITC = '#e08214', '#b03a2e'
RHO = 0.351
_fp = lambda p: f'{p:.2g}' if p >= 1e-4 else f'{p:.0e}'


class Fig:                                                       # same layout helper as make_fig2_fig3_reorg.py
    def __init__(self, w, hgt, row_h, gaps, top=.885, bot=.10):
        self.fig = plt.figure(figsize=(w, hgt)); self.W, self.H = w, hgt
        gaps = [g / hgt for g in gaps]; unit = (top - bot - sum(gaps)) / sum(row_h)
        self.rows = []; y = top
        for i, rh in enumerate(row_h):
            b = y - unit * rh; self.rows.append((y, b)); y = b - (gaps[i] if i < len(gaps) else 0)
        self.AX = {}; self.sections = []

    def section(self, row, x0, x1, letters, title, wspace=.35, wr=None):
        t, b = self.rows[row]
        gs = self.fig.add_gridspec(1, len(letters), left=x0, right=x1, top=t, bottom=b, wspace=wspace, width_ratios=wr)
        axs = []
        for i, L in enumerate(letters):
            ax = self.fig.add_subplot(gs[0, i]); self.AX[L] = ax; axs.append(ax)
        self.sections.append((letters, title, row, x0))
        return axs

    def finish(self, suptitle, note, path):
        fig = self.fig; fig.canvas.draw(); rend = fig.canvas.get_renderer(); F = fig.transFigure.inverted()
        top = {}
        for L, ax in self.AX.items():
            p = ax.get_position()
            xl = min(p.x0, F.transform(ax.yaxis.get_tightbbox(rend).p0)[0]) if ax.yaxis.get_visible() else p.x0
            t = fig.text(xl - .012, p.y1 + .016, L, size=11, weight='bold', va='bottom', ha='left')
            yt = t.get_window_extent(rend).transformed(F).y1
            if ax.title.get_text():
                yt = max(yt, ax.title.get_window_extent(rend).transformed(F).y1)
            top[L] = yt
        row_top = {}
        for letters, _, r, _ in self.sections:
            row_top[r] = max(row_top.get(r, 0), max(top[L] for L in letters))
        for n, (letters, text, r, x0) in enumerate(self.sections):
            fig.text(x0 - .03, row_top[r] + .010, f'{n + 1}.  {text}', size=8.8, weight='bold', color=INK, ha='left', va='bottom')
        fig.text(.52, .05, note, ha='center', va='top', size=6.3, color=GREY, linespacing=1.6)
        fig.suptitle(suptitle, size=13, weight='bold', y=.975)
        fig.savefig(path.with_suffix('.pdf')); fig.savefig(path.with_suffix('.png'), dpi=200)
        print('wrote', path.with_suffix('.pdf'))


# ------------------------------------------------------------------ data
L = pd.read_csv(LOC / 'growth_grid_conga.csv'); L = L[L.gene_count == 'adj'].copy()
L['trait'] = L.trait.str.lower(); L['measure'] = L.measure.astype(str)
G = pd.read_csv(RA / 'growthloci_2026-11-15/grid_summary_all.tsv', sep='\t')
G['measure'] = G.variant.str.split('|').str[0]; G['adj'] = G.variant.str.split('|').str[1]
G = G[(G.adj == 'adj') & (G.readout.str.split('|').str[0] == G.measure)].copy()
G['readout'] = G.readout.str.split('|').str[-1]; G['trait'] = G.trait.str.replace('_rint', '')
G['measure'] = G.measure.str.replace('^g:', '', regex=True)


def family(m):
    return 'con' if m in ('none_H', 'n95', 'pres:n95', 'pres:n95HI') else ('uns' if ':uns:' in m else 'sig')


SEL = {'height': ['HC:sum:uns:scaled', 'none_H'], 'bmi': ['HC:c4:uns:mean', 'none_H', 'pres:n95']}   # genome-wide selected measures
COLS = [('height', 'DEL'), ('height', 'DUP'), ('height', 'mirror'), ('bmi', 'DEL'), ('bmi', 'DUP'), ('bmi', 'mirror')]


def screen(ax, D, title, held_frac=.9, ylim=(-7.6, 7.6)):
    rng = np.random.default_rng(5)
    for i, (t, ro) in enumerate(COLS):
        d = D[(D.trait == t) & (D.readout == ro)]
        for fam, col, off in [('sig', SIGC, -.24), ('uns', UNSC, 0), ('con', CONC, .24)]:
            q = d[d.measure.map(family) == fam]
            held = (q.loo_held >= held_frac * q.loo_n).to_numpy()
            x = i + off + rng.uniform(-.07, .07, len(q))
            ax.scatter(x[held], q.z_emp[held], s=13, color=col, edgecolor='white', lw=.3, zorder=4)
            ax.scatter(x[~held], q.z_emp[~held], s=11, facecolor='white', edgecolor=col, lw=.8, zorder=3)
            sel = q.measure.isin(SEL[t]).to_numpy() & ((fam != 'con') | (ro != 'mirror'))
            ax.scatter(x[sel], q.z_emp[sel], s=58, facecolor='none', edgecolor=INK, lw=1.0, zorder=5)
    for v in (-1.96, 1.96):
        ax.axhline(v, color=GREY, lw=.5, ls=':', zorder=1)
    ax.axhline(0, color=LIGHT, lw=.7, zorder=0); ax.axvline(2.5, color=INK, lw=.6)
    ax.set_xticks(range(6)); ax.set_xticklabels(['DEL', 'DUP', 'mirror'] * 2)
    ax.set_xlim(-.55, 5.55); ax.set_ylim(*ylim)
    for x0, t in [(1, 'Height'), (4, 'BMI')]:
        ax.text(x0, ylim[1] * .97, t, ha='center', va='top', size=8, weight='bold', color=HCOL if t == 'Height' else BCOL)
    ax.set_ylabel('z (genome-shift null)'); ax.set_title(title, size=8.4, pad=5, weight='bold')


def screen_key(ax, loc='lower left', ring='ring: selected for the joint model'):
    h = [Line2D([], [], marker='o', ls='', color=SIGC, ms=4.5, label='signed TWAS measures (20)'),
         Line2D([], [], marker='o', ls='', color=UNSC, ms=4.5, label='unsigned TWAS measures (20)'),
         Line2D([], [], marker='o', ls='', color=CONC, ms=4.5, label='constraint measures (4)'),
         Line2D([], [], marker='o', ls='', mfc='white', mec=GREY, ms=4.5, label='open: not held in removals'),
         Line2D([], [], marker='o', ls='', mfc='none', mec=INK, ms=7.5, label=ring)]
    ax.legend(handles=h, frameon=False, fontsize=6.2, loc=loc, handletextpad=.3, labelspacing=.25, borderaxespad=.2)


def forest(ax, rows, title, xlim=(-7.5, 7.5), right_text=None):
    """rows: list of (label, [(z, held_bool, colour, marker, offset, kind)])"""
    y = 0; yt = []; yl = []
    for lab, pts in rows:
        for z, held, col, mk, off, kind in pts:
            if z is None or np.isnan(z):
                continue
            if kind == 'open':
                ax.scatter(z, y + off, s=26, marker=mk, facecolor='white', edgecolor=col, lw=1.0, zorder=4)
            else:
                ax.scatter(z, y + off, s=30, marker=mk, color=col if held else 'white', edgecolor=col, lw=1.0, zorder=4)
        yt.append(y); yl.append(lab); y -= 1
    for v in (-1.96, 1.96):
        ax.axvline(v, color=GREY, lw=.5, ls=':')
    ax.axvline(0, color=INK, lw=.6)
    ax.set_yticks(yt); ax.set_yticklabels(yl, size=6.8); ax.set_ylim(y + .4, .6); ax.set_xlim(*xlim)
    ax.set_xlabel('z (joint model)'); ax.set_title(title, size=8.4, pad=5, weight='bold')
    ax.tick_params(axis='y', length=0)


# ------------------------------------------------------------------ locus-level joint and selection numbers
A = pd.read_csv(LOC / 'growth_conga_axes_summary.csv'); A['trait'] = A.disorder.str.lower()
ax_ = lambda t, r: A[(A.trait == t) & (A.readout == r)].iloc[0]
PV = pd.read_csv(LOC / 'growth_conga_presence_vs_count.csv'); PV['trait'] = PV.trait.str.lower()
pv = lambda t, d, f: PV[(PV.trait == t) & (PV.dosage == d) & (PV.form == f)].iloc[0]
MY = pd.read_csv(V2 / 'yaxis_monogenic.csv'); MY['trait'] = MY.trait.str.lower()
hit = {t: MY[(MY.trait == t) & (MY.set == 'top5_or_hi') & (MY.model == 'alone')].iloc[0] for t in ['height', 'bmi']}
J = {'height': pd.read_csv(RA / 'joint_sacks34_2026-11-15/height_height_selected.tsv', sep='\t'),
     'bmi': pd.read_csv(RA / 'joint_sacks34_2026-11-15/bmi_bmi_selected.tsv', sep='\t')}
jz = lambda t, ro: J[t][J[t].readout == ro].iloc[0]
# Jon's LOO rule (2026-11-15): genome-wide = matched null per driver + one matched null dropping all minor units (drivers kept). The 24 single-unit
# removals above already include the 8 drivers; the minor_all null comes from the drivers-mode rerun (UKBB session, Expanse job 55744860, same seed/inputs)
for t in J:
    D_ = pd.read_csv(RA / f'joint_drivers_2026-11-15/{t}/{t}_selected/grid_summary.tsv', sep='\t')[['readout', 'z_emp', 'minor_all_z_emp']]
    J[t] = J[t].merge(D_.rename(columns={'z_emp': 'z_emp_drv'}), on='readout', how='left')
    assert np.allclose(J[t].z_emp, J[t].z_emp_drv), t                   # same full-sample fit in both runs


def heldA(r):
    return r.loo_held >= .9 * r.loo_n


# ================================================================== figure
F = Fig(15.0, 12.6, [1.0, 1.0, 1.08], [1.05, 1.05], top=.895, bot=.11)
LX = (.06, .975)
a = F.section(0, *LX, 'ABCD', 'Recurrent CNVs (34 intervals): univariate screen, selection of the constraint and single-hit terms, joint model', wspace=.5,
              wr=[1.3, .72, .95, .95])
screen(a[0], L, 'Univariate screen: 44 measures × 3 readouts', ylim=(-5.2, 5.2)); screen_key(a[0], 'lower left', 'ring: measure selected genome-wide (F, G)')
# B: constraint form
# B (Jon 2026-11-15): height duplication offset ~ mirror line + one gene term; presence vs count for each gene set, and the union.
#   two_axis_CNV_dose_response/v2_2026-11-08/code/growth_yaxis_presence_vs_count.R (REML meta-regression, ML AIC)
ax = a[1]
PC = pd.read_csv(V2 / 'yaxis_presence_vs_count.csv'); PC = PC[PC.trait == 'Height']
SET_COL = {'con': CONC, 'top5': '#e0533d', 'hi': '#f5a03a', 'union': '#8b1a10'}
ROWS_B = [('con', 'count', 'constrained genes, count'), ('top5', 'count', 'top-5% s$_{het}$, count'), ('top5', 'presence', 'top-5% s$_{het}$, present'),
          ('hi', 'presence', 'ClinGen HI, present'), ('union', 'count', 'union, count'), ('union', 'presence', 'union, present')]
for i, (st, fm, lab) in enumerate(ROWS_B):
    r = PC[(PC.set == st) & (PC.form == fm)].iloc[0]; gain = -r.dAIC_vs_mirror_only; held = r.loo_interval_held >= .9 * r.loci
    ax.barh(-i, gain, .62, color=SET_COL[st] if fm == 'presence' else 'white', edgecolor=SET_COL[st], lw=1.3, hatch=None if fm == 'presence' else '////', zorder=3)
    ax.text(gain + .5, -i, f'z {r.z:+.1f}' + ('' if held else ' (not held)'), va='center', size=6.4, color=INK)
ax.set_yticks([-i for i in range(len(ROWS_B))]); ax.set_yticklabels([r[2] for r in ROWS_B], size=6.6); ax.tick_params(axis='y', length=0)
ax.set_xlim(0, 30); ax.set_ylim(-len(ROWS_B) + .4, .6); ax.axvline(0, color=INK, lw=.6)
ax.set_xlabel('AIC improvement over the mirror line')
hh = PC[(PC.set == 'union') & (PC.form == 'presence + count')].iloc[0]
ax.text(29.5, .45, f'union, both terms\nin one model:\npresent z {hh.z_presence_joint:+.1f}\ncount z {hh.z_count_joint:+.1f}', ha='right', va='top', size=6.2, color=INK)
ax.set_title('Height DUP offset: presence vs count', size=8.4, pad=5, weight='bold')
# C: locus joint model
rows = []
for lab, key in [('mirror (unsigned\nTWAS count)', 'mirrorU'), ('TWAS count, DEL', 'countU_DEL'), ('TWAS count, DUP', 'countU_DUP'),
                 ('constraint count, DEL', 'conU_DEL'), ('constraint count, DUP', 'conU_DUP')]:
    rows.append((lab, [(ax_(t, key).z_emp, heldA(ax_(t, key)), c, 'o', o, 'f') for t, c, o in [('height', HCOL, .13), ('bmi', BCOL, -.13)]]))
rows.append(('single hit, DUP\n(top-5% s$_{het}$ or HI gene)', [(hit[t].z_diff, hit[t].loo_held >= 31, c, 'D', o, 'f') for t, c, o in [('height', HCOL, .13), ('bmi', BCOL, -.13)]]))
forest(a[3], rows, 'Joint model, recurrent CNVs', xlim=(-5.5, 5.5))
# C: deletion vs duplication effect of each interval, by whether it carries a top-5% s_het gene and/or a ClinGen HI gene (report fig. 3)
ax = a[2]
LH = pd.read_csv(LOC / 'growth_conga_loci.csv'); LH = LH[LH.trait == 'Height'].merge(
    pd.read_csv(V2 / 'yaxis_monogenic_intervals.csv').query("trait == 'Height'")[['NAME', 'top5', 'hi_clingen', 'top5_or_hi']], on='NAME', how='left')
NEITHER = '#2b6cb0'; SHC = '#e0533d'; HIC = '#f5a03a'; BOTHC = '#8b1a10'   # Jon 2026-11-15: no rings; s_het red, HI orange, both dark red
for sel, col in [((LH.top5 == 0) & (LH.hi_clingen == 0), NEITHER), ((LH.top5 == 1) & (LH.hi_clingen == 0), SHC),
                 ((LH.top5 == 0) & (LH.hi_clingen == 1), HIC), ((LH.top5 == 1) & (LH.hi_clingen == 1), BOTHC)]:
    d = LH[sel]
    ax.errorbar(d.b_DEL, d.b_DUP, xerr=d.se_DEL, yerr=d.se_DUP, fmt='none', ecolor=col, alpha=.35, lw=.6, zorder=2)
    ax.scatter(d.b_DEL, d.b_DUP, s=28, color=col, edgecolor='white', lw=.4, zorder=4)
assert len(LH) == 34 and LH.top5_or_hi.sum() == 11
xx = np.array([-1.3, .6]); ax.plot(xx, -RHO * xx, color=GREY, lw=.9, ls=(0, (4, 3)), zorder=1)
ax.axhline(0, color=LIGHT, lw=.7, zorder=0); ax.axvline(0, color=LIGHT, lw=.7, zorder=0)
LAB = {'16p11.2_BP4_BP5': ('16p11.2 BP4–5', 0, .065, 'center'), '7q11.23_WBS': ('7q11.23 (WBS)', 0, .07, 'center'),
       '22q11.21_A_B': ('22q11.2 A–B', 0, -.065, 'center'), '22q11.21_A-D': ('22q11.2 A–D', 0, .055, 'center'), '17q12': ('17q12', -.05, -.03, 'right'),
       '15q13.1-13.2_BP3_BP4': ('15q13 BP3–4', .04, .045, 'left'), '2q11.2': ('2q11.2', .05, -.04, 'left'), '17p12': ('17p12', -.06, -.06, 'right'),
       '16p11.2_BP1-BP3': ('16p11.2 BP1–3', 0, -.12, 'center'), '7q11.23_distal': ('7q11.23 distal', .06, .06, 'left'), '16p11.2_BP2_BP3': ('16p11.2 BP2–3', -.05, 0, 'right')}
for nm, (t, dx, dy, ha) in LAB.items():
    r = LH[LH.NAME == nm].iloc[0]; ax.text(r.b_DEL + dx, r.b_DUP + dy, t, size=5.9, ha=ha, va='center', color=INK)
ax.set_xlim(-1.3, .6); ax.set_ylim(-.62, .62)
ax.set_xlabel('Deletion effect on height (SD)'); ax.set_ylabel('Duplication effect on height (SD)')
h = hit['height']
ax.set_title(f'Single hit: {int(h.intervals_with)}/{int(h.intervals)} intervals, DUP {h["diff"]:+.2f} SD', size=8.4, pad=5, weight='bold')
ax.legend(handles=[Line2D([], [], marker='o', ls='', color=SHC, ms=4.5, label='top-5% s$_{het}$ gene'),
                   Line2D([], [], marker='o', ls='', color=HIC, ms=4.5, label='ClinGen HI gene (score 2–3)'),
                   Line2D([], [], marker='o', ls='', color=BOTHC, ms=4.5, label='both'),
                   Line2D([], [], marker='o', ls='', color='#2b6cb0', ms=4.5, label='neither (no single-hit gene)')],
          frameon=False, fontsize=6.0, loc='upper right', handletextpad=.2, borderaxespad=.1)
a[3].legend(handles=[Line2D([], [], marker='o', ls='', color=HCOL, ms=4.5, label='height'), Line2D([], [], marker='o', ls='', color=BCOL, ms=4.5, label='BMI'),
                     Line2D([], [], marker='o', ls='', mfc='white', mec=GREY, ms=4.5, label='open: not held')],
            frameon=False, fontsize=6.2, loc='lower right', handletextpad=.2)

b = F.section(1, *LX, 'EFG', f'Genome-wide ({SRC}): univariate screen, then joint models with and without the recurrent intervals', wspace=.42,
              wr=[1.45, .92, .92])
screen(b[0], G, 'Univariate screen: 44 measures × 3 readouts'); screen_key(b[0], 'lower left')
for axx, t, col, title, terms in [
        (b[1], 'height', HCOL, 'Joint model, height', [('mirror (HC genes, unsigned\nTWAS sum)', 'HC:sum:uns:scaled|mirror'), ('  … DEL arm', 'HC:sum:uns:scaled|DEL'),
                                                       ('  … DUP arm', 'HC:sum:uns:scaled|DUP'), ('constraint count, DEL', 'none_H|DEL'),
                                                       ('constraint count, DUP', 'none_H|DUP'), ('genes, DEL', 'n_all|DEL'), ('genes, DUP', 'n_all|DUP')]),
        (b[2], 'bmi', BCOL, 'Joint model, BMI', [('mirror (HC genes, unsigned\nTWAS count |z| > 4)', 'HC:c4:uns:mean|mirror'), ('  … DEL arm', 'HC:c4:uns:mean|DEL'),
                                                 ('  … DUP arm', 'HC:c4:uns:mean|DUP'), ('top-5% s$_{het}$ gene\npresent, DEL', 'pres:n95|DEL'),
                                                 ('constraint count, DUP', 'none_H|DUP'), ('genes, DEL', 'n_all|DEL'), ('genes, DUP', 'n_all|DUP')])]:
    rows = []
    # Jon 2026-11-16: the gene count (mean shift) is not rotated by the shift null and absorbs the scrambled annotation terms under it (UKBB BMI DEL:
    # t +4.7 vs null mean +7.4 -> z_emp -4.3, opposite to the raw effect). Its rows show the analytic t (against 0), per run, from per_run_z.tsv.
    PR = pd.read_csv(RA / f'joint_drivers_2026-11-15/{t}/{t}_selected/per_run_z.tsv', sep='\t')
    for lab, ro in terms:
        r = J[t][J[t].readout.str.replace('g:', '', regex=False) == ro].iloc[0]
        if ro.startswith('n_all|'):
            p = PR[PR.readout.str.replace('g:', '', regex=False) == ro].set_index('run').t
            lo = p[p.index.str.startswith('loo_')]; held = ((lo.abs() > 1.96) & (np.sign(lo) == np.sign(p['full']))).mean() >= .9
            rows.append((lab, [(p['full'], held, col, 'o', .2, 'f'), (p['minor_all'], True, col, '^', 0, 'f'), (p['no_loci'], True, col, 'x', -.2, 'f')]))
            continue
        rows.append((lab, [(r.z_emp, r.loo_held >= .9 * r.loo_n, col, 'o', .2, 'f'), (r.minor_all_z_emp, True, col, '^', 0, 'f'),
                           (r.noloci_z_emp, True, col, 'x', -.2, 'f')]))
    forest(axx, rows, title)
    axx.legend(handles=[Line2D([], [], marker='o', ls='', color=col, ms=4.5, label='all CNVs'), Line2D([], [], marker='o', ls='', mfc='white', mec=col, ms=4.5, label='all CNVs, not held'),
                        Line2D([], [], marker='^', ls='', color=col, ms=4.5, label='minor recurrent units removed\n(8 drivers kept)'),
                        Line2D([], [], marker='x', ls='', color=col, ms=5, label='34 recurrent intervals removed')], frameon=False, fontsize=6.0, loc='lower left' if t == 'height' else 'lower right', handletextpad=.2)


# ------------------------------------------------------------------ schematic (final solution)
def arrow(ax, x0, y0, x1, y1, col, held, circle=False):
    lw = 3.0 if held else 1.2
    if circle:
        ax.plot([x0, x1], [y0, y1], color=col, lw=lw, zorder=4, solid_capstyle='butt')
        ax.plot(x1, y1, 'o', ms=9 if held else 7, color=col if held else 'white', mec=col, mew=1.2, zorder=5)
    else:
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>' if held else '->', mutation_scale=14 if held else 11, color=col,
                                     lw=lw, shrinkA=0, shrinkB=0, zorder=4, fill=held))


def plane(ax, title):
    lim = 7.2
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
    ax.axhline(0, color=LIGHT, lw=.7, zorder=0); ax.axvline(0, color=LIGHT, lw=.7, zorder=0)
    xx = np.array([-lim, lim]); ax.plot(xx, -RHO * xx, color=GREY, lw=.9, ls=(0, (4, 3)), zorder=1)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel('Deletion effect'); ax.set_ylabel('Duplication effect')
    ax.set_title(title, size=8.4, pad=5, weight='bold')
    for (x, y, t, ha, va) in [(-lim * .97, lim * .97, 'DUP raises,\nDEL lowers', 'left', 'top'), (lim * .97, lim * .97, 'both raise', 'right', 'top'),
                              (lim * .97, -lim * .97, 'DEL raises,\nDUP lowers', 'right', 'bottom'), (-lim * .97, -lim * .97, 'both lower', 'left', 'bottom')]:
        ax.text(x, y, t, ha=ha, va=va, size=6.2, color=GREY)


def txt(ax, x, y, s, col, ha='left', va='center'):
    ax.text(x, y, s, ha=ha, va=va, size=6.7, color=col, weight='bold', linespacing=1.2)


u = np.array([-1, RHO]) / np.hypot(1, RHO)                     # unit vector along the mirror, DUP raises
# genome-wide single hits (Jon 2026-11-15): presence of a top-5% s_het or ClinGen HI gene, fitted with the constrained-gene count
# (additivity_constraint.py, hit models; analytic z); thick only when P < 0.05 with and without the 34 recurrent intervals
HM = RA / 'additivity_hit3_2026-11-15'


def gw_hit(t, term, hdef='top5_hi'):
    out = {}
    for ex in ['all', 'noSacks34']:
        m = pd.read_csv(HM / t / f'{hdef}_{ex}' / 'additivity_models.tsv', sep='\t')
        r = m[(m.model == 'hit: HC count + hit presence') & (m.term == term)].iloc[0]; out[ex] = r
    za, zn = float(out['all']['T']), float(out['noSacks34']['T'])
    return za, bool((out['all']['P'] < .05) and (out['noSacks34']['P'] < .05) and np.sign(za) == np.sign(zn)), zn


c = F.section(2, .12, .66, 'HI', 'Final solution: the gene dosage architecture of height and BMI', wspace=.25)
# Jon 2026-11-15: one panel per trait, held components only, short labels. Arrow length = z of the established test:
#   mirror: recurrent intervals (height, D) / genome-wide joint (BMI, G); constraint: genome-wide joint (F, G); single hit: union presence (B)
# Jon 2026-11-15: the full mirror = TWAS + rho. Long orange arrow = full mirror (z of the DEL-DUP locus correlation in the bivariate joint model
# without the TWAS term, single hit in the model); short yellow arrow overlaid = the part TWAS explains (twas_share x full length; share = 1 - cov_with/cov_without).
# two_axis_CNV_dose_response/v2_2026-11-08/code/growth_mirror_decomp.R (unsigned TWAS |z|>3 count, union-presence hit)
DEC = pd.read_csv(V2 / 'joint_bivariate_2026-11-15/mirror_decomposition_union_presence.csv'); DEC = DEC[DEC.twas == 'twU']
TWC = '#f2c12e'


def mirror_pair(ax, t, lab_off, ha, va):
    d = DEC[DEC.trait == t].iloc[0]; L = abs(d.z_total); s = max(d.twas_share, 0)
    arrow(ax, 0, 0, *(L * u), MIRC, True)
    ax.add_patch(FancyArrowPatch((0, 0), tuple(s * L * u), arrowstyle='-|>', mutation_scale=11, color=TWC, lw=2.2, shrinkA=0, shrinkB=0, zorder=6))
    txt(ax, *(L * u + [.3, .45]), 'mirror', MIRC, 'center', 'bottom')
    if va == 'below':
        txt(ax, *(s * L * u + [-.3, -.45]), f'TWAS {100 * s:.0f}%', '#b8860b', 'center', 'top')
    else:
        txt(ax, *(s * L * u + [.25, .35]), f'TWAS\n{100 * s:.0f}%', '#b8860b', 'center', 'bottom')


# Jon 2026-11-15 ("Locus level only"): H and I are drawn from ONE model per trait, the 34-interval bivariate joint model
#   (growth_joint_bivariate.R; mirror decomposition growth_mirror_decomp.R). Only components that hold are drawn: height = mirror + single hit DUP
#   (hit_DUP z_emp, 16/16 families, 34/34 intervals); constraint is not significant at locus level for either trait (it is genome-wide, F and G); BMI = mirror only.
BV = pd.read_csv(V2 / 'joint_bivariate_2026-11-15/joint_bivariate_summary.csv'); BV = BV[BV.twas.str.startswith('unsigned')]
bv = lambda t, ro: BV[(BV.trait == t) & (BV.readout == ro)].iloc[0]
ax = c[0]; plane(ax, 'Height')
mirror_pair(ax, 'Height', [-.2, .35], 'right', 'bottom')
z = bv('Height', 'hit_DUP').z_emp; arrow(ax, .7, 0, .7, z, HITC, True, circle=True); txt(ax, 1.1, z + .2, 'single hit', HITC)
ax = c[1]; plane(ax, 'BMI')
mirror_pair(ax, 'BMI', [-.2, .35], 'right', 'below')
# legend
ax = F.fig.add_axes([.71, F.rows[2][1], .2, F.rows[2][0] - F.rows[2][1]]); ax.axis('off'); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_autoscale_on(False)
items = [(MIRC, False, 'Mirror: deletion and\nduplication effects opposite\n(total effect); yellow =\nTWAS portion'),
         (HITC, True, 'Single hit: presence of\none top-5% s$_{het}$ or\nhaploinsufficient gene')]
for i, (col, circ, s) in enumerate(items):
    y = .93 - i * .22
    if circ:
        ax.plot([.0, .16], [y, y], color=col, lw=3); ax.plot(.16, y, 'o', ms=8, color=col)
    else:
        ax.add_patch(FancyArrowPatch((0, y), (.18, y), arrowstyle='-|>', mutation_scale=13, color=col, lw=3, transform=ax.transAxes))
    ax.text(.24, y, s, va='center', size=6.8, color=INK, transform=ax.transAxes, linespacing=1.2)
ax.text(0, .42, 'One model per trait: the\n34-interval joint model.\nArrow length = z of the\ncomponent; dashed line =\nmirror line. Constraint is\nnot significant here (it is\ngenome-wide: F, G).', size=6.4, color=GREY, transform=ax.transAxes, va='top')

NOTE = ('A–D: 34 recurrent CNV intervals with published height and BMI effects (Sacks et al.). E–G: ' + SRC + (' genome-wide CNV calls (All of Us to be added). ' if SRC == 'UK Biobank' else ' genome-wide CNV calls (meta-analysis). ') + (
        'Measures: TWAS genes (all or highly constrained, HC), summed z or counts at |z| > 1–4, signed or unsigned, mean or tissue-scaled; constraint = count or presence of constrained genes.\n'
        'Readouts: deletion arm, duplication arm and mirror ((z$_{DUP}$ − z$_{DEL}$)/√2). z against burden-matched circular shifts of the gene annotation (A: 1,000; E–G: 1,000–2,000); '
        'held = |z| > 1.96 with the same sign in ≥ 90% of single-locus removals (34 intervals in A–D; 24 recurrent units, overlapping intervals merged, in E–G).\n'
        'F, G: matched nulls for each of the 24 single-unit removals, for removal of all 16 minor units together with the 8 driver units kept (triangles) '
        'and for removal of all recurrent intervals (crosses).\n'
        'B: duplication effect ~ deletion effect (mirror line) + one gene term, REML meta-regression, ML AIC; union = top-5% s$_{het}$ or ClinGen HI (11 intervals); hatched = count, solid = presence; '
        'no term fits the BMI duplication offset except the constrained-gene count. D: analytic z. F, G: the gene-count rows (mean shift) show the analytic z, because the shift null does not rotate the gene count.\n'
        'C: bars ±1 SE; dashed line = mirror. H, I: components that hold; mirror (orange) = total mirror effect, z of the deletion–duplication locus correlation in the recurrent-interval joint model '
        '(single hit in the model); yellow = TWAS portion, the share removed by adding the unsigned TWAS count (height 58%, BMI 17%); single hit = duplication single-hit term of the same model; '
        'constraint terms of that model are not significant and are not drawn.'))
F.finish('Dissecting the dose response: univariate screen, selection and joint models', NOTE, OUTDIR / 'Fig_models')
