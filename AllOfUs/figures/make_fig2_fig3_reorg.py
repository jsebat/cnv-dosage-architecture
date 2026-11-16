#!/usr/bin/env python
"""Figures 2 and 3 reorganised (Jon 2026-11-15). Every panel of the manuscript Figure 2 (Figure_height_main_v3_ad, A-O) is kept and
drawn by its own code, copied unchanged from Molly/ContinuousVsSingleHit/perp_offset_bimodality/code/make_height_figure_v3.py (same
inputs, same rcParams, colours and section style). New UK Biobank panels are drawn in the same style from aggregate tables.

Story: the height dose response of CNVs is a mixture of three components: mirror, constraint (deletion axis) and a monogenic single
hit (duplication axis).
  Figure 2  Mirror and constraint. Genome-wide CNV burden is non-linear, while recurrent CNVs show strong mirror effects (Fig. 1F).
            The mirror is the additive effect of many genes (gene count, TWAS) and is carried by large recurrent CNVs; constrained
            genes add a deletion-axis component; mirror and constraint scores add up across CNVs, for height and BMI.
  Figure 3  Asymmetric CNVs are monogenic. Bimodal offsets, class differences, a single constrained gene on the duplication axis
            (presence, not count; more sub-additive), single-hit simulations, and the two components (mirror + single hit) explain
            the U-shaped genome-wide burden.
Old -> new panel map: old A,B -> 3A,B; old C-F -> 3C-F; old G,H -> 2C,D; old I,J -> 2E,F; old K,L -> 3J,K; old N,O -> 3L,M;
old M -> 3N.

AllOfUs copy (2026-11-16): the genome-wide panels (2A-B, 2G-L, 3H-I, 3O-P) read GENOMEWIDE_ROOT, by default the UK Biobank aggregates.
Point it at the output of meta/meta_sacks2026.py to draw the same panels from the UK Biobank + All of Us meta-analysis:
  GENOMEWIDE_ROOT=<meta dir> GENOMEWIDE_LABEL="UK Biobank + All of Us" FIG_OUT=<dir> python make_fig2_fig3_reorg.py
The locus panels read the lab Dropbox (SEBATLAB_DROPBOX, default ~/sebatlab Dropbox). Writes Fig2_reorg.{pdf,png}, Fig3_reorg.{pdf,png}.
"""
from __future__ import annotations
import os, pickle, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

OUTDIR = Path(os.environ.get('FIG_OUT', Path(__file__).resolve().parent))
DB = Path(os.environ.get('SEBATLAB_DROPBOX', Path.home() / 'sebatlab Dropbox'))
MOLLY = DB / 'Manuscripts/PGC_joint_submission/Molly'
HERE = MOLLY / 'ContinuousVsSingleHit/perp_offset_bimodality'
CODE = HERE / 'code'
sys.path.insert(0, str(CODE))
IN = MOLLY / 'ContinuousVsSingleHit/CNV_height_BMI_corrected_directional_models_v2'
OUT = HERE / 'outputs' / 'me/all_loci'; OUT_CRIT = OUT; LOCI = IN
CONC = MOLLY / 'CNV_ASD_growth_analysis' / 'concordance'
TWAS_SET = 'ad'
JC = DB / 'Jon/Claude'
RA = Path(os.environ.get('GENOMEWIDE_ROOT', JC / 'UKBB_anthro_CNV/results_aggregate'))
V2 = JC / 'two_axis_CNV_dose_response/v2_2026-11-08/results/eight_traits/growth_conga'
ALT = 'equal_var'

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial'], 'font.size': 8,
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
    'mathtext.bf': 'Arial:bold',
    'axes.titlesize': 8.4, 'axes.labelsize': 8, 'axes.linewidth': .7,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
    'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42,
})
INK = '#1a1a1a'; GREY = '#8a8a8a'; LIGHT = '#dcdcdc'
BLUE = '#1f5c8b'; RED = '#b03a2e'; GOLD = '#c47b1a'; PURPLE = '#6a3d9a'
HCOL, BCOL = INK, GOLD          # height, BMI in the genome-wide panels

# ------------------------------------------------------------------ inputs of the original figure (v3 lines 111-137)
h = pd.read_csv(LOCI / 'height' / 'inputs' / 'sacks_height_all_loci_exact.csv')
po = pd.read_csv(OUT / 'perp_offsets_per_locus.csv')
po = po[po.trait == 'Height'].tail(len(h))
h = h.merge(po[['CNV_ID', 'perp_offset']], on='CNV_ID')
h['cls'] = h.Dose_response_class.str.lower()
MIR = h[h.cls == 'mirror']; ASY = h[h.cls == 'asymmetric']
GRP = [('mirror', BLUE, 'o'), ('asymmetric', RED, 's')]
res = pickle.loads((OUT / '_res.pkl').read_bytes())
R = [r for r in res if r['mirror_subset'] and r['name'] == 'Height'][0]
tab = pd.read_csv(OUT / 'bimodality_test.csv')
ROW = tab[(tab.trait == 'Height') & (tab.mirror_fit == 'mirror-class only') & (tab.alternative == ALT)].iloc[0]
sc = pd.read_csv(OUT_CRIT / 'model_criteria_scores.csv')
BIMODAL = ROW.p_bootstrap_nullB < .05
mace = pd.read_csv(IN / 'code' / 'mace_burden_reference.csv')
mace = mace[mace.trait == 'height'].set_index('burden_model')
x = R['x']; f = R['fits'][ALT]
cls = np.asarray([str(s).lower() for s in R['cls']])
nullB, o = f['nullB'], f['obs']
_fp = lambda p: f'{p:.2g}' if p >= 1e-4 else f'{p:.0e}'


def dens(x_, mu, sd):
    return np.exp(-0.5 * ((x_ - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))


# ------------------------------------------------------------------ layout helpers (v3 style: numbered sections, measured headings)
class Fig:
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
            xl = min(p.x0, F.transform(ax.yaxis.get_tightbbox(rend).p0)[0]) if ax.axison else p.x0
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
        fig.text(.52, .055, note, ha='center', va='top', size=6.3, color=GREY, linespacing=1.6)
        fig.suptitle(suptitle, size=13, weight='bold', y=.975)
        fig.savefig(path.with_suffix('.pdf')); fig.savefig(path.with_suffix('.png'), dpi=200)
        print('wrote', path.with_suffix('.pdf'))


# ================================================================== original panel code (v3, unchanged apart from ax handles)
def panel_class_diff(axs):                                   # v3 lines 173-198 (old C-F)
    rng = np.random.default_rng(3)
    PANELS = [('CNV_length_Mb', 'CNV length (Mb)', np.log, True), ('Gene_count', 'Gene count', np.log1p, True),
              ('Deletion_beta', 'Deletion effect (SD)', None, False), ('Duplication_beta', 'Duplication effect (SD)', None, False)]
    for i, (v, lab, tr, from_zero) in enumerate(PANELS):
        ax = axs[i]
        a = ASY[v].to_numpy(float); m = MIR[v].to_numpy(float)
        aa, mm = (tr(a), tr(m)) if tr is not None else (a, m)
        p = stats.ttest_ind(aa, mm, equal_var=False).pvalue
        for j, (vals, cc) in enumerate([(m, BLUE), (a, RED)]):
            ax.scatter(np.full(len(vals), j) + rng.normal(0, .085, len(vals)), vals, s=9, color=cc, alpha=.8, lw=0)
            ax.plot([j - .30, j + .30], [np.median(vals)] * 2, color=INK, lw=1.5, zorder=4)
        span = max(a.max(), m.max()) - min(a.min(), m.min())
        if from_zero:
            ax.set_ylim(0, max(a.max(), m.max()) * 1.12)
        else:
            ax.set_ylim(min(a.min(), m.min()) - .07 * span, max(a.max(), m.max()) + .12 * span)
            ax.axhline(0, color=LIGHT, lw=.6, zorder=1)
        ax.set_xlim(-.62, 1.62); ax.set_xticks([0, 1]); ax.set_xticklabels(['mir', 'asy'])
        ax.set_ylabel(lab, size=7.6)
        ptxt = f'$p$ = {p:.3f}' if p >= 1e-4 else f'$p$ = {p:.0e}'
        ax.set_title(ptxt, size=8, pad=5, color=RED if p < .05 else INK, weight='bold' if p < .01 else 'normal')


def panel_bimodal(axA, axB):                                  # v3 lines 200-250 (old A-B)
    ax = axA
    lo, hi = x.min() - .10, x.max() + .10
    ax.hist(x, bins=np.linspace(lo, hi, 15), density=True, color=BLUE, alpha=.20, edgecolor='white', lw=.5)
    g = np.linspace(lo, hi, 500)
    S = R['me']['S'][:, None]; mu1, t1 = R['me']['one']; pi_, m1, m2, t2 = R['me']['two']
    one = dens(g[None, :], mu1, np.sqrt(t1 ** 2 + S ** 2)).mean(0)
    two = (pi_ * dens(g[None, :], m1, np.sqrt(t2 ** 2 + S ** 2)) + (1 - pi_) * dens(g[None, :], m2, np.sqrt(t2 ** 2 + S ** 2))).mean(0)
    ax.plot(g, one, color=GREY, lw=1.5, ls='--', label='1 component')
    ax.plot(g, two, color=BLUE, lw=2.0, label='2 components')
    ymax = ax.get_ylim()[1]
    for lab, cc in [('mirror', BLUE), ('asymmetric', RED)]:
        s = cls == lab
        ax.plot(x[s], np.full(int(s.sum()), -ymax * .06), '|', color=cc, ms=7, mew=1.3, clip_on=False, label=f'{lab} ({int(s.sum())})')
    ax.set_ylim(-ymax * .11, ymax * 1.05)
    ax.set_xlabel('Perpendicular offset from mirror line (SD)'); ax.set_ylabel('Density')
    ax.set_title('Offsets are bimodal' if BIMODAL else 'Offsets: bimodality not significant', size=8.4, pad=5, weight='bold')
    ax.legend(frameon=False, loc='upper left', handletextpad=.5, labelspacing=.22, fontsize=6.8)
    ax = axB
    top = max(np.quantile(nullB, .995), o) * 1.10
    ax.hist(nullB, bins=np.linspace(0, top, 40), density=True, color=GREY, alpha=.45, edgecolor='white', lw=.4)
    ax.axvline(o, color=BLUE, lw=2.2); ax.set_xlim(0, top)
    ax.set_xlabel('Likelihood ratio, 2 vs 1 component'); ax.set_ylabel('Density')
    ax.set_title(f'vs unimodal null:  $p$ = {ROW.p_bootstrap_nullB:.3f}', size=8.4, pad=5, weight='bold', color=RED if BIMODAL else INK)
    _r = o / top
    _x, _y, _ha = (.97, .92, 'right') if not .66 < _r < .99 else (_r - .02, .62, 'right')
    ax.text(_x, _y, f'observed {o:.1f}\n{len(nullB):,} null draws', transform=ax.transAxes, ha=_ha, va='top', size=7, color=INK)


def panel_genecount(axs):                                     # v3 lines 252-282 (old G-H)
    for c, (eff, elab) in enumerate([('Deletion_beta', 'Deletion effect (SD)'), ('Duplication_beta', 'Duplication effect (SD)')]):
        ax = axs[c]; v = 'Gene_count'
        for cls_, col, mk in GRP:
            z = h[h.cls == cls_]
            ax.scatter(z[v], z[eff], s=20, marker=mk, color=col, alpha=.85, edgecolor='white', lw=.4)
            sl, ic = stats.theilslopes(z[eff], z[v])[:2]
            xx = np.linspace(z[v].min(), z[v].max(), 40)
            rr, pp = stats.spearmanr(z[eff], z[v])
            ax.plot(xx, ic + sl * xx, color=col, lw=1.9 if pp < .05 else 1.0, ls='-' if pp < .05 else (0, (4, 2)), zorder=4)
        sl, ic = stats.theilslopes(h[eff], h[v])[:2]
        xx = np.linspace(h[v].min(), h[v].max(), 40)
        ra, pa = stats.spearmanr(h[eff], h[v])
        ax.plot(xx, ic + sl * xx, color=INK, lw=1.4, ls=(0, (1, 1.6)), zorder=2)
        ax.axhline(0, color=LIGHT, lw=.7)
        ax.set_xlabel('Gene count'); ax.set_ylabel(elab)
        ents = [('all ', INK, ra, pa)]
        for cls_, col, _ in GRP:
            z = h[h.cls == cls_]
            rr, pp = stats.spearmanr(z[eff], z[v])
            ents.append((cls_[:4], col, rr, pp))
        y0 = .30 if eff == 'Deletion_beta' else .99
        for k, (lbl, col, rr, pp) in enumerate(ents):
            ax.text(.03, y0 - k * .092, f'{lbl}  $\\rho$={rr:+.2f}  $p$={pp:.3f}', transform=ax.transAxes, ha='left', va='top',
                    size=6.9, color=col, weight='bold' if pp < .05 else 'normal')
        ax.set_title('deletion' if c == 0 else 'duplication', size=8.4, pad=5)


def panel_twas(axI, axJ):                                     # v3 lines 423-486 (old I-J)
    _SUB = ['16p11.2_BP2_BP3', '1q21.1_TAR_BP1_BP2', '1q21.1_TAR_BP2_BP3', '16p13.11_BP1_BP2', '16p13.11_BP2_BP3']
    TW_DROP = _SUB + ['22q11.21_A_B', '22q11.21_B_C', '22q11.21_C_D', '22q11.21_B-D']
    TW_SETNAME = 'primary (spanning intervals; 22q A-D)'
    tw = pd.read_csv(CONC / 'twas_dosage_panel.csv'); tw = tw[~tw.CNV_ID.isin(TW_DROP)].copy()
    tst = pd.read_csv(CONC / 'twas_nonoverlapping.csv'); tst = tst[tst.cnv_set == TW_SETNAME].set_index('test')
    WB = (tw.CNV_ID == '7q11.23_WBS').to_numpy()
    ax = axI
    for cls_, col, mk in GRP:
        s = (tw.cls == cls_).to_numpy() & ~WB
        ax.scatter(tw.summed_twas_z[s], tw.dosage_effect[s], s=20, marker=mk, color=col, alpha=.85, edgecolor='white', lw=.4, zorder=4)
    ax.scatter(tw.summed_twas_z[WB], tw.dosage_effect[WB], s=22, marker='s', facecolor='white', edgecolor=RED, lw=1.0, zorder=4)
    _sl, _ic = stats.theilslopes(tw.dosage_effect[~WB], tw.summed_twas_z[~WB])[:2]
    _xx = np.linspace(tw.summed_twas_z[~WB].min(), tw.summed_twas_z[~WB].max(), 40)
    ax.plot(_xx, _ic + _sl * _xx, color=INK, lw=1.6, zorder=3)
    ax.axhline(0, color=LIGHT, lw=.7); ax.axvline(0, color=LIGHT, lw=.7)
    _LAB = {'16p11.2_BP4_BP5': ('16p11.2 BP4–BP5', 0, .06, 'center'), '3q29': ('3q29', 0, -.08, 'center'),
            '1q21.1': ('1q21.1', -.5, 0, 'right'), '22q11.21_A-D': ('22q11.2 A–D', -.45, .06, 'right'),
            '22q11.21_B-D': ('22q11.2 B–D', -.45, .06, 'right'), '7q11.23_WBS': ('7q11.23 (Williams)', .55, 0, 'left'),
            '16p11.2_BP1-BP3': ('16p11.2 BP1–BP3', 0, .07, 'center')}
    for _cid, (_t, _dx, _dy, _ha) in _LAB.items():
        _r = tw[tw.CNV_ID == _cid]
        if len(_r):
            _r = _r.iloc[0]
            ax.text(_r.summed_twas_z + _dx, _r.dosage_effect + _dy, _t, size=6.6, ha=_ha, va='center', color=BLUE if _r.cls == 'mirror' else RED)
    _ra = tst.loc['dosage rho, all without Williams']
    ents = [('all ', INK, _ra.estimate, _ra.p_size_matched),
            ('mirr', BLUE, tst.loc['dosage rho, mirror'].estimate, tst.loc['dosage rho, mirror'].p_size_matched),
            ('asym', RED, tst.loc['dosage rho, asymmetric without Williams'].estimate, tst.loc['dosage rho, asymmetric without Williams'].p_size_matched)]
    for k, (lbl, col, rr, pp) in enumerate(ents):
        ax.text(.03, .99 - k * .085, f'{lbl}  $\\rho$={rr:+.2f}  $p$={_fp(pp)}', transform=ax.transAxes, ha='left', va='top',
                size=6.9, color=col, weight='bold' if pp < .05 else 'normal')
    ax.text(.03, .99 - 3 * .085, f'separate lines for the two classes: $p$={_fp(tst.loc[[i for i in tst.index if i.startswith("separate")][0]].p_welch)}',
            transform=ax.transAxes, ha='left', va='top', size=6.6, color=GREY)
    ax.set_ylim(min(-.5, tw.dosage_effect.min() - .08), max(1.45, tw.dosage_effect.max() + .78))
    ax.set_xlabel('Summed TWAS z of genes in the CNV\n(common-variant predicted expression; + = higher expression, taller)')
    ax.set_ylabel('CNV dosage effect on height (SD)\n(DUP − DEL) / 2')
    ax.set_title('Predicted expression (TWAS) predicts dosage effects', size=8.4, pad=5, weight='bold')
    ax = axJ
    _pts = pd.read_csv(CONC / 'twas_class_strip_points_all.csv'); _pts = _pts[~_pts.CNV_ID.isin(TW_DROP)]
    _ss = pd.read_csv(CONC / 'twas_class_strip_stats_sets.csv'); _ss = _ss[_ss.twas_set == TWAS_SET].set_index('group')
    _rng2 = np.random.default_rng(3)
    for j, (cls_, col) in enumerate([('mirror', BLUE), ('asymmetric', RED)]):
        v = _pts.size_matched_z[_pts.cls == cls_].to_numpy(float)
        ax.scatter(np.full(len(v), j) + _rng2.normal(0, .085, len(v)), v, s=9, color=col, alpha=.8, lw=0)
        ax.plot([j - .30, j + .30], [np.median(v)] * 2, color=INK, lw=1.5, zorder=4)
        ax.text(j, 1.02, f'$p$ = {_fp(_ss.loc[cls_, "p_two_sided"])}', ha='center', va='bottom', size=6.9, color=INK)
    ax.axhline(0, color=LIGHT, lw=.6, zorder=1)
    ax.set_xlim(-.62, 1.62); ax.set_xticks([0, 1]); ax.set_xticklabels(['mir', 'asy']); ax.set_ylim(-1.9, 1.35)
    ax.set_ylabel('Summed TWAS z vs equal-gene-count\nrandom windows (SD)', size=7.6)
    ax.set_title(f'asy vs mir: $p$ = {_fp(_ss.loc["asymmetric minus mirror", "p_two_sided"])}', size=8, pad=5, color=INK)


MODEL_STYLE = [('mirror', 'Mirror alone', '#4d4d4d', (0, (4, 2))), ('open_continuous', '+ continuous', GOLD, '-'),
               ('open_single', '+ single hit', PURPLE, '-')]


def panel_models(axM, axN, axO):                              # v3 lines 284-396 (old M-O)
    import score_models_on_criteria as SC
    ax = axM
    mb = float(mace.loc['duplication'].beta); mlo = float(mace.loc['duplication'].ci_low_derived)
    mhi = float(mace.loc['duplication'].ci_high_derived)
    XMAX = 4.0; xb = np.linspace(0, XMAX, 50)
    ax.fill_between(xb, mlo * xb, mhi * xb, color=GREY, alpha=.22, lw=0, zorder=1)
    ax.plot(xb, mb * xb, color=INK, lw=1.8, zorder=3); ax.axhline(0, color=LIGHT, lw=.8, zorder=0)
    _scm = sc.set_index('key'); hd = []
    for key, lab, col, ls in MODEL_STYLE:
        s = float(_scm.loc[key, 'burden_DUP_per_Mb'])
        ax.plot(xb, s * xb, color=col, lw=2.0, ls=ls, zorder=4)
        hd.append(Line2D([], [], color=col, lw=2.0, ls=ls, label=f'{lab}  {s:+.3f}'))
    hd.append(Line2D([], [], color=INK, lw=1.8, label=f'observed  {mb:+.3f}'))
    hd.append(Patch(color=GREY, alpha=.22, lw=0, label='observed 95% CI'))
    ax.legend(handles=hd, frameon=False, loc='upper left', fontsize=6.1, handlelength=2.0, labelspacing=.22, title='slope, cm/Mb',
              title_fontsize=6.1, alignment='left')
    ax.set_xlim(0, XMAX); ax.set_ylim(-.62, 1.15)
    ax.set_xlabel('Duplication burden (Mb)'); ax.set_ylabel('Predicted change in height (cm)')
    ax.set_title('DUP burden and short stature', size=8.4, pad=5, weight='bold')
    _cmp = pd.read_csv(LOCI / 'height' / 'outputs' / 'height_axis_mixture_model_comparison.csv').set_index('model')
    _base = pd.read_csv(LOCI / 'height' / 'outputs' / 'height_mirror_parameters.csv').iloc[0].to_dict()
    _rho = _base['rho']; _k = np.sqrt(1 + _rho * _rho)
    assert abs(_rho - R['rho']) < 1e-6
    _G = h.Gene_count.to_numpy(float)
    _draws = np.load(OUT_CRIT / 'model_criteria_draws.npz')
    _rng = np.random.default_rng(SC.SEED)
    SIMX, E = {}, {}
    for key, lab, col, ls in MODEL_STYLE:
        pars = _cmp.loc[SC.LABEL_IN_CSV[key]].to_dict()
        Ds, Us = SC.simulate(key, pars, _G, h.Deletion_SE.to_numpy(float), h.Duplication_SE.to_numpy(float), _base, -1,
                             len(_draws[key + '_lrt']), _rng)
        SIMX[key] = ((_rho * Ds + Us) / _k).ravel()
        pv = (nullB[None, :] >= _draws[key + '_lrt'][:, None]).mean(1)
        E[key] = float(np.median(pv))
    ax = axN
    lo, hi = -.62, .32; g = np.linspace(lo, hi, 400)
    ax.hist(x, bins=np.linspace(x.min() - .10, x.max() + .10, 15), density=True, color=BLUE, alpha=.16, edgecolor='white', lw=.5, zorder=1)
    for key, lab, col, ls in MODEL_STYLE:
        ax.plot(g, stats.gaussian_kde(SIMX[key], bw_method=.12)(g), color=col, lw=2.0, ls=ls, zorder=4)
    ymax = ax.get_ylim()[1]
    for lab, cc in [('mirror', BLUE), ('asymmetric', RED)]:
        s = cls == lab
        ax.plot(x[s], np.full(int(s.sum()), -ymax * .06), '|', color=cc, ms=7, mew=1.3, clip_on=False)
    ax.set_ylim(-ymax * .11, ymax * 1.05); ax.set_xlim(lo, hi)
    ax.set_xlabel('Perpendicular offset (SD)'); ax.set_ylabel('Density')
    ax.set_title('Expected offset distribution', size=8.4, pad=5, weight='bold')
    hd = [Line2D([], [], color=BLUE, lw=6, alpha=.25, label='observed')]
    hd += [Line2D([], [], color=col, lw=2.0, ls=ls, label=lab) for key, lab, col, ls in MODEL_STYLE]
    ax.legend(handles=hd, frameon=False, loc='upper left', fontsize=6.1, handlelength=2.0, labelspacing=.22)
    ax = axO
    _all = np.concatenate([nullB] + [_draws[k + '_lrt'] for k in SIMX])
    top = max(np.quantile(_all, .99), o) * 1.05; gl = np.linspace(0, top, 400)
    ax.hist(nullB, bins=np.linspace(0, top, 36), density=True, color=GREY, alpha=.35, edgecolor='white', lw=.3, zorder=1)
    for key, lab, col, ls in MODEL_STYLE:
        lr = _draws[key + '_lrt']; kd = stats.gaussian_kde(np.concatenate([lr, -lr]), bw_method=.10)
        ax.plot(gl, 2 * kd(gl), color=col, lw=2.0, ls=ls, zorder=4)
    ax.axvline(o, color=BLUE, lw=2.0, zorder=3.5); ax.set_xlim(0, top)
    ax.set_xlabel('Likelihood ratio, 2 vs 1 component'); ax.set_ylabel('Density')
    ax.set_title('Expected bimodality', size=8.4, pad=5, weight='bold')
    hd = [Line2D([], [], color=BLUE, lw=2.0, label='observed'), Patch(color=GREY, alpha=.35, lw=0, label='unimodal null')]
    hd += [Line2D([], [], color=col, lw=2.0, ls=ls, label=f'{lab}: median $p$ = {E[key]:.3f}') for key, lab, col, ls in MODEL_STYLE]
    ax.set_ylim(0, ax.get_ylim()[1] * 1.12)
    leg = ax.legend(handles=hd, frameon=True, loc='upper right', fontsize=6.0, handlelength=2.0, labelspacing=.22, facecolor='white',
                    edgecolor='none', framealpha=1, borderpad=.3)
    leg.set_zorder(10)


def panel_clouds(axs):                                        # v3 lines 398-421 (old K-L)
    import make_supp_model_clouds as MC
    CL = MC.cloud_data()
    for _c, (_key, _name) in enumerate(MC.SHOW):
        ax = axs[_c]
        MC.draw_cloud(ax, _key, *CL[_key]); MC.observed_key(ax)
        _row = MC.cmp_.loc[MC.SC.LABEL_IN_CSV[_key]]
        _sub = f'$\\Delta$AIC {_row.delta_conditional_AIC:.1f}'
        if _key == 'open_single':
            _sub += f';  P(hit | 23 genes) = {1 - (1 - _row.q) ** 23:.2f}'
        ax.set_title(f'{_name}\n{_sub}', size=8.6, pad=5, weight='bold', linespacing=1.5)


# ================================================================== new panels (UK Biobank aggregates; same style)
SRC = os.environ.get('GENOMEWIDE_LABEL', 'UK Biobank')
XL = ['21+', '11–20', '6–10', '3–5', '1–2', '0', '1–2', '3–5', '6–10', '11–20', '21+']
SB = {t: pd.read_csv(RA / f'dup_burden_2026-11-15/{t}__signed_burden_bins.tsv', sep='\t') for t in ['height', 'bmi']}
SM = {t: pd.read_csv(RA / f'dup_burden_2026-11-15/{t}__signed_burden_models.tsv', sep='\t') for t in ['height', 'bmi']}
TL = {'height': 'Height', 'bmi': 'BMI'}


def burden_axes(ax, trait):
    ax.axhline(0, color=LIGHT, lw=.7, zorder=0); ax.axvline(0, color=LIGHT, lw=.7, zorder=0)
    ax.set_xticks(range(-5, 6)); ax.set_xticklabels(XL, size=5.6); ax.set_xlim(-5.6, 5.6)
    ax.set_xlabel('genes deleted            genes duplicated'); ax.set_ylabel(f'{TL[trait]} (SD), adjusted')


def burden_line(ax, trait, grp, col, ls, mk, lab, dx=0.0, dup_only=False):
    d = SB[trait][SB[trait].group == grp].sort_values('order')
    if dup_only:
        d = d[d.order >= 0]
    ax.errorbar(d.order + dx, d.mean_SD, yerr=1.96 * d.SE, color=col, ls=ls, marker=mk, ms=3.6, lw=1.3, capsize=1.6, label=lab,
                mec='white', mew=.4, zorder=4)


def p_u(trait, grp='all'):
    m = SM[trait]; r = m[(m.group == grp) & (m.term == 's^2 (U-shape)')].iloc[0]
    return r.beta_SD, r.P


def panel_burden(ax, trait):
    b, p = p_u(trait)
    burden_line(ax, trait, 'all', HCOL if trait == 'height' else BCOL, '-', 'o', None)
    burden_axes(ax, trait)
    ax.set_title(f'{TL[trait]}: genome-wide burden', size=8.4, pad=5, weight='bold')
    ax.text(.5, .04 if b < 0 else .96, f'quadratic term: $p$ = {_fp(p)}', transform=ax.transAxes, ha='center',
            va='bottom' if b < 0 else 'top', size=6.9, color=INK)


# grid (genome-wide mirror and constraint, with and without the 34 recurrent intervals)
G = pd.read_csv(RA / 'growthloci_2026-11-15/grid_summary_all.tsv', sep='\t')
G['measure'] = G.variant.str.split('|').str[0]; G['adj'] = G.variant.str.split('|').str[1]
G['ro'] = G.readout.str.split('|').str[-1]; G['trait'] = G.trait.str.replace('_rint', '')
CS = pd.read_csv(RA / 'growthloci_2026-11-15/callsize_summary_all.tsv', sep='\t'); CS['trait'] = CS.trait.str.replace('_rint', '')


def panel_grid_mirror(ax):
    d0 = G[(G.adj == 'adj') & (G.ro == 'mirror') & G.measure.str.startswith('g:') & (G.readout.str.split('|').str[0] == G.measure)]
    for t, col, mk in [('height', HCOL, 'o'), ('bmi', BCOL, 's')]:
        d = d0[d0.trait == t]
        ax.scatter(d.z_emp, d.noloci_z_emp, s=13, marker=mk, color=col, alpha=.8, edgecolor='white', lw=.3, label=f'{TL[t]} ({len(d)} measures)', zorder=4)
    lim = 7.5
    ax.plot([-lim, lim], [-lim, lim], color=LIGHT, lw=.7, zorder=1)
    for v in (-1.96, 1.96):
        ax.axhline(v, color=GREY, lw=.5, ls=':'); ax.axvline(v, color=GREY, lw=.5, ls=':')
    ax.axhline(0, color=LIGHT, lw=.6); ax.axvline(0, color=LIGHT, lw=.6)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel('Mirror z, all CNVs'); ax.set_ylabel('Mirror z, 34 recurrent intervals removed')
    ax.set_title('Genome-wide mirror is carried by recurrent CNVs', size=8.4, pad=5, weight='bold')
    ax.legend(frameon=False, loc='upper left', fontsize=6.6, handletextpad=.3)


def panel_callsize(ax):
    ms = [('g:HC:c3:uns:mean', 'bmi', 'BMI\nunsigned count\nHC genes |z|>3'), ('g:all:sum:uns:scaled', 'bmi', 'BMI\nunsigned\nsummed z'),
          ('g:HC:c3:sig:mean', 'height', 'Height\nsigned count\nHC genes |z|>3')]
    for i, (mm, t, _) in enumerate(ms):
        for j, (cl_, col) in enumerate([('L5', INK), ('S5', '#b8b8b8')]):
            r = CS[(CS.readout == f'{cl_}~{mm}|mirror') & (CS.trait == t)].iloc[0]
            ax.bar(i + (j - .5) * .36, r.z_emp, .34, color=col, label=('calls > 5 genes' if cl_ == 'L5' else 'calls ≤ 5 genes') if i == 0 else None)
    ax.axhline(0, color=INK, lw=.6)
    for v in (-1.96, 1.96):
        ax.axhline(v, color=GREY, lw=.5, ls=':')
    ax.set_xticks(range(3)); ax.set_xticklabels([m[2] for m in ms], size=6.4); ax.set_ylabel('Mirror z (genome-shift null)')
    ax.set_ylim(-5.4, 3.4); ax.set_title('Mirror comes from large CNVs', size=8.4, pad=5, weight='bold')
    ax.legend(frameon=False, fontsize=6.6, loc='upper right')


ORDER_HC = ['0', '1', '2', '3-4', '5-9', '10+']


def panel_hc_dose(ax, trait):
    b = pd.read_csv(RA / f'additivity_2026-11-15/{trait}__dose_response_bins.tsv', sep='\t')
    for dos, col, mk, dx in [('DEL', RED, 'o', -.09), ('DUP', BLUE, 's', .09)]:
        g = b[(b.dosage == dos) & (b.group == 'single CNV')].set_index('HC_bin').reindex(ORDER_HC)
        ax.errorbar(np.arange(6) + dx, g.mean_SD, yerr=1.96 * g.SE, color=col, marker=mk, ms=3.6, lw=1.3, capsize=1.6, mec='white', mew=.4,
                    label='deletion' if dos == 'DEL' else 'duplication', zorder=4)
    ax.axhline(0, color=LIGHT, lw=.7, zorder=0)
    ax.set_xticks(range(6)); ax.set_xticklabels([s.replace('-', '–') for s in ORDER_HC])
    ax.set_xlabel('Constrained genes in the CNV (s$_{het}$ > 0.029)'); ax.set_ylabel(f'{TL[trait]} (SD), adjusted')
    ax.set_title(f'{TL[trait]}', size=8.4, pad=5, weight='bold')
    ax.legend(frameon=False, fontsize=6.6, loc='lower left' if trait == 'height' else 'upper left')


def panel_constraint_grid(ax):
    d0 = G[(G.measure == 'none_H') & (G.adj == 'adj') & (G.readout.str.startswith('none_H|'))]
    d0 = d0.drop_duplicates(['trait', 'readout'])
    items = [('height', 'DEL'), ('height', 'DUP'), ('bmi', 'DEL'), ('bmi', 'DUP')]
    for i, (t, dos) in enumerate(items):
        r = d0[(d0.trait == t) & (d0.ro == dos)].iloc[0]
        col = HCOL if t == 'height' else BCOL
        ax.bar(i - .18, r.z_emp, .34, color=col, label='all CNVs' if i == 0 else None)
        ax.bar(i + .18, r.noloci_z_emp, .34, color=col, alpha=.45, hatch='////', edgecolor='white', lw=0,
               label='34 recurrent intervals removed' if i == 0 else None)
    ax.axhline(0, color=INK, lw=.6)
    for v in (-1.96, 1.96):
        ax.axhline(v, color=GREY, lw=.5, ls=':')
    ax.set_xticks(range(4)); ax.set_xticklabels(['Height\nDEL', 'Height\nDUP', 'BMI\nDEL', 'BMI\nDUP'], size=6.6)
    ax.set_ylabel('Constrained-gene count z (genome-shift null)'); ax.set_ylim(-5.6, 4.6)
    ax.set_title('Constraint does not depend on recurrent CNVs', size=8.4, pad=5, weight='bold')
    ax.legend(frameon=False, fontsize=6.4, loc='upper left')


def panel_additivity(ax):
    rows = []
    for t, f_ in [('height', 'height_height_selected'), ('bmi', 'bmi_bmi_selected')]:
        a = pd.read_csv(RA / f'score_axes_2026-11-15/{f_}__additivity.tsv', sep='\t').set_index('axis')
        for axis in ['mirror', 'constraint']:
            r = a.loc[axis]; s = r.slope_single
            rows.append((t, axis, (s + r.slope_multi_minus_single) / s, r.SE_diff / s, r.P_interaction))
    y = 0; yt = []; yl = []
    for t, axis, est, se, p in rows:
        col = HCOL if t == 'height' else BCOL
        ax.errorbar(est, y, xerr=1.96 * se, fmt='o' if axis == 'mirror' else 's', color=col, ms=4.5, capsize=1.8, lw=1.3, mec='white', mew=.4)
        ax.text(1.95, y, f'$p$ = {_fp(p)}', va='center', ha='right', size=6.7, color=RED if p < .05 else INK)
        yt.append(y); yl.append(f'{TL[t]}, {axis}'); y -= 1
    ax.axvline(1, color=INK, lw=.7); ax.axvline(0, color=LIGHT, lw=.6)
    ax.set_yticks(yt); ax.set_yticklabels(yl, size=7); ax.set_ylim(y + .4, .6); ax.set_xlim(-.1, 2.0)
    ax.set_xlabel('Score effect in carriers of 2+ CNVs\nrelative to carriers of one CNV (1 = additive)')
    ax.set_title('Effects add up across CNVs', size=8.4, pad=5, weight='bold')


CLS = pd.read_csv(MOLLY / 'pathway_analysis/Final/cnv_genes_annotated.csv')[['CNV_ID', 'Dose_response_class']].drop_duplicates().set_index('CNV_ID').Dose_response_class


def panel_hit_offset(ax):
    Y = pd.read_csv(V2 / 'yaxis_monogenic_intervals.csv'); Y = Y[Y.trait == 'Height'].copy(); Y['cls'] = Y.NAME.map(CLS).str.lower()
    rng = np.random.default_rng(3)
    for j, flag in enumerate([0, 1]):
        d = Y[Y.top5_or_hi == flag]
        for c_, col, mk in GRP:
            dd = d[d.cls == c_]
            ax.scatter(np.full(len(dd), j) + rng.normal(0, .085, len(dd)), dd.off, s=14, marker=mk, color=col, alpha=.85, edgecolor='white', lw=.3, zorder=4)
        ax.plot([j - .30, j + .30], [np.median(d.off)] * 2, color=INK, lw=1.5, zorder=5)
    MY = pd.read_csv(V2 / 'yaxis_monogenic.csv'); r = MY[(MY.trait == 'Height') & (MY.set == 'top5_or_hi') & (MY.model == 'alone')].iloc[0]
    ax.axhline(0, color=LIGHT, lw=.7, zorder=0)
    ax.set_xlim(-.62, 1.62); ax.set_xticks([0, 1]); ax.set_xticklabels(['no', 'yes'])
    ax.set_xlabel('Interval contains a top-5% s$_{het}$\nor haploinsufficient gene'); ax.set_ylabel('DUP offset from the mirror line (SD)')
    ax.set_title(f'difference {r["diff"]:+.2f} SD, $p$ = {_fp(2 * stats.norm.sf(abs(r.z_diff)))}', size=8, pad=5, color=RED, weight='bold')


def panel_presence(ax):
    for nm, col, mk, lab, dx in [('onlyasym', RED, 's', 'asymmetric DUPs', -.09), ('onlymirror', BLUE, 'o', 'mirror DUPs', .09)]:
        b = pd.read_csv(RA / f'additivity_2026-11-15/height__{nm}_only_DUP__dose_response_bins.tsv', sep='\t')
        g = b[(b.dosage == 'DUP') & (b.group == 'single CNV')].set_index('HC_bin').reindex(ORDER_HC); ok = (g.n >= 20).to_numpy()
        ax.errorbar(np.arange(6)[ok] + dx, g.mean_SD[ok], yerr=1.96 * g.SE[ok], color=col, marker=mk, ms=3.6, lw=1.3, capsize=1.6,
                    mec='white', mew=.4, label=lab, zorder=4)
    ax.axhline(0, color=LIGHT, lw=.7, zorder=0)
    ax.set_xticks(range(6)); ax.set_xticklabels([s.replace('-', '–') for s in ORDER_HC])
    ax.set_xlabel('Constrained genes duplicated'); ax.set_ylabel('Height (SD), adjusted')
    ax.set_title('Presence, not number, of a constrained gene', size=8.4, pad=5, weight='bold')
    ax.legend(frameon=False, fontsize=6.6, loc='lower left')


def panel_subadd(ax):
    CM = pd.read_csv(RA / 'additivity_2026-11-15/subadditivity_asym_vs_mirror_2026-11-15.tsv', sep='\t')
    CM = CM[(CM.trait == 'height') & (CM.model == 'linear') & (CM.carriers == 'both dosages')]
    y = 0; yt = []; yl = []
    for row in CM.itertuples():
        for grp, s, col, mk, off in [('asymmetric', row.asym, RED, 's', .14), ('mirror', row.mirror, BLUE, 'o', -.14)]:
            est = float(s.split(' ')[0]); se = float(s.split('SE ')[1].split(',')[0])
            ax.errorbar(est, y + off, xerr=1.96 * se, fmt=mk, color=col, ms=4.2, capsize=1.6, lw=1.2, mec='white', mew=.4)
        ax.text(.148, y, f'difference\n$p$ = {_fp(row.P_diff)}', ha='right', va='center', size=6.5, color=INK)
        yt.append(y); yl.append('constrained genes\n' + ('deleted' if row.term.endswith('DEL') else 'duplicated')); y -= 1
    ax.axvline(0, color=INK, lw=.7)
    ax.set_yticks(yt); ax.set_yticklabels(yl, size=7); ax.set_ylim(y + .45, .55); ax.set_xlim(-.04, .15)
    ax.set_xlabel('2+ CNVs × constrained-gene count\n(SD per gene; > 0 = less than additive)')
    ax.set_title('Asymmetric CNVs are less additive', size=8.4, pad=5, weight='bold')
    ax.legend(handles=[Line2D([], [], marker='s', ls='', color=RED, label='asymmetric'), Line2D([], [], marker='o', ls='', color=BLUE, label='mirror')],
              frameon=False, fontsize=6.6, loc='center right')


def panel_u_resolved(ax):
    burden_line(ax, 'height', 'all', INK, '-', 'o', 'all people')
    burden_line(ax, 'height', 'asymmetric-interval carriers removed', RED, (0, (4, 2)), 's', 'asymmetric intervals removed', dx=.12, dup_only=True)
    burden_line(ax, 'height', 'HC-gene DUP carriers removed', BLUE, '-', 'D', 'DUPs of constrained genes removed', dx=.24, dup_only=True)
    burden_axes(ax, 'height')
    _, p = p_u('height', 'HC-gene DUP carriers removed')
    ax.set_title('Removing single hits removes the DUP arm', size=8.4, pad=5, weight='bold')
    ax.legend(frameon=False, fontsize=6.3, loc='lower center', handlelength=2.2)


def panel_dup_slopes(ax):
    grps = [('all', 'none'), ('asymmetric-interval carriers removed', 'asymmetric\nintervals'),
            ('hit-gene DUP carriers removed', 'top-5% s$_{het}$\nor HI gene'), ('HC-gene DUP carriers removed', 'any\nconstrained')]
    for t, col, mk, off in [('height', HCOL, 'o', -.13), ('bmi', BCOL, 's', .13)]:
        for i, (g_, _) in enumerate(grps):
            r = SM[t][(SM[t].group == g_) & (SM[t].term == 'DUP genes (per 10)')].iloc[0]
            ax.errorbar(i + off, r.beta_SD, yerr=1.96 * r.SE, fmt=mk, color=col, ms=4.3, capsize=1.7, lw=1.3, mec='white', mew=.4,
                        label=TL[t] if i == 0 else None)
    ax.axhline(0, color=INK, lw=.7)
    ax.set_xticks(range(4)); ax.set_xticklabels([g[1] for g in grps], size=6.4)
    ax.set_xlabel('Duplication carriers removed'); ax.set_ylabel('Effect per 10 genes duplicated (SD)')
    ax.set_title('Duplication burden effect', size=8.4, pad=5, weight='bold')
    ax.set_ylim(top=.045); ax.legend(frameon=False, fontsize=6.6, loc='upper left', ncol=2)


# ================================================================== Figure 2
F2 = Fig(15.0, 10.4, [1.0, 1.15, 1.05], [1.05, 1.15], top=.885, bot=.115)
LX, RX = (.072, .475), (.565, .988)
a = F2.section(0, *LX, 'AB', 'Genome-wide CNV burden is non-linear, while recurrent CNVs have mirror effects (Fig. 1F)', wspace=.32)
panel_burden(a[0], 'height'); panel_burden(a[1], 'bmi')
a = F2.section(0, *RX, 'CD', 'Within mirror CNVs, effect size scales with gene count', wspace=.30)
panel_genecount(a)
a = F2.section(1, *LX, 'EF', 'Common-variant gene effects predict CNV dosage effects', wspace=.42, wr=[1.55, 1.0])
panel_twas(*a)
a = F2.section(1, *RX, 'GH', 'Genome-wide, the mirror is carried by large recurrent CNVs', wspace=.36)
panel_grid_mirror(a[0]); panel_callsize(a[1])
a = F2.section(2, .072, .655, 'IJK', 'Constraint: CNVs of constrained genes shorten stature and raise BMI', wspace=.36)
panel_hc_dose(a[0], 'height'); panel_hc_dose(a[1], 'bmi'); panel_constraint_grid(a[2])
a = F2.section(2, .775, .988, 'L', 'Mirror and constraint are additive')
panel_additivity(a[0])
NOTE2 = (f'A, B, G–L: {SRC}' + (' (N = 311,874; All of Us to be added).' if SRC == 'UK Biobank' else ' (meta-analysis).') + ' A, B: carriers of deletions only or duplications only, binned by genes '
         'deleted (left) or duplicated (right); bars = 95% CI. C–F: 33 recurrent CNV loci (Sacks et al.); blue = mirror, red = asymmetric.\n'
         'E: dosage effect = (DUP − DEL)/2; TWAS z summed over the genes in the CNV (Mancuso et al. 2017); open square = Williams–Beuren region '
         '(7q11.23), excluded from the line and the correlations; $p$ from circular genomic permutation. F: bars = medians; non-overlapping CNVs only.\n'
         'G, H, K: z against 1,000 burden-matched circular shifts of the gene annotation. I, J: carriers of one CNV. '
         'L: fitted per-person mirror and constraint scores; slope in carriers of 2+ CNVs divided by the slope in carriers of one; bars = 95% CI.')
F2.finish('Mirror and constraint: two additive components of CNV dose response', NOTE2, OUTDIR / 'Fig2_reorg')

# ================================================================== Figure 3
F3 = Fig(15.0, 10.6, [1.0, 1.25, 1.05], [1.05, 1.25], top=.885, bot=.11)
a = F3.section(0, *LX, 'AB', 'The perpendicular offset from the mirror line is bimodal', wspace=.30)
panel_bimodal(*a)
a = F3.section(0, *RX, 'CDEF', 'Mirror and asymmetric CNVs differ by size and DUP effect, not DEL', wspace=.62)
panel_class_diff(a)
a = F3.section(1, *LX, 'GHI', 'Asymmetric duplications act through a single constrained gene', wspace=.55, wr=[.8, 1.15, 1.1])
panel_hit_offset(a[0]); panel_presence(a[1]); panel_subadd(a[2])
a = F3.section(1, *RX, 'JK', 'What the two directional models produce', wspace=.30)
panel_clouds(a)
a = F3.section(2, .072, .40, 'LM', 'Only the single-hit model produces two groups', wspace=.36)
b = F3.section(2, .475, .988, 'NOP', 'Mirror + single hit explains the U-shaped genome-wide burden', wspace=.40, wr=[.8, 1.35, .9])
panel_models(b[0], a[0], a[1]); panel_u_resolved(b[1]); panel_dup_slopes(b[2])
NOTE3 = ('A–G, J–N: 33 recurrent CNV loci with published height effects (Sacks et al.); blue = mirror (21), red = asymmetric (12). '
         f'H, I, O, P: {SRC}' + (' (All of Us to be added).' if SRC == 'UK Biobank' else ' (meta-analysis).') + ' G: DUP offset below the mirror line by whether the interval contains a top-5% s$_{het}$ '
         'or ClinGen haploinsufficient gene.\nH: carriers of one duplication. I: interaction of carrying 2+ CNVs with the constrained-gene count. '
         'L, M: each model’s 2,000 simulated datasets (real gene counts and SEs); median $p$ = the typical dataset’s $p$ against the unimodal null in B.\n'
         'N: each model is calibrated so its deletion burden equals the published −0.41 cm/Mb (Macé et al. 2017); its duplication slope is an '
         'out-of-sample prediction; grey band = published 95% CI. O, P: bars = 95% CI.')
F3.finish('Asymmetric CNVs are monogenic: a rare single hit on the duplication axis', NOTE3, OUTDIR / 'Fig3_reorg')
