#!/usr/bin/env python
"""Dosage paper, opening figures (Jon 2026-11-22, author notes 1-3 in the working Doc): the genome-wide story before the CNV GWAS section,
as three separate figures for now.
  Fig. 1  The paradox. a,b single recurrent loci have mirror effects (adult UK Biobank + All of Us meta); c,d genome-wide burden of all
          genes is U-shaped.
  Fig. 2  Dissection by TWAS, constraint and CNV size. a,b genes with weak TWAS effects: U-shaped, carried by constrained genes;
          c-h genes with strong TWAS effects (z > +3, z < -3, |z| > 3), less-constrained vs constrained; i-l the same components by the
          number of genes the CNV call hits (1, 2-5, 6-20, > 20): mirror z of the signed TWAS count in less-constrained and constrained
          genes, and the constrained-gene effect of deletions and duplications.
  Fig. 3  Joint model of gene dosage architecture. a univariate heatmap (each term + total gene count); b joint model (six classes that
          partition all genes); c the joint classes by CNV size; d schematic of the components per trait (after Shanta et al. Fig. 3I).
Inputs (aggregate tables only, so the same code draws UK Biobank alone or a meta-analysis across cohorts):
  GW2_ROOT/class_bins/<trait>/class_dose_bins.tsv, class_dose_shape.tsv   (class_dose_bins.py; analytic OLS, null model: none)
  GW2_GRID/<trait>/<variant>/grid_summary.tsv                              (conga_grid_qt.py, design class_models_2026-11-22.json;
                                                                            null model: genome shift of the gene annotation, 1,000 shared
                                                                            replicates; LOO = each major recurrent locus removed)
  GW2_LOCI/meta_locus_effects_<trait>.tsv                                  (meta_locus_effects.py, UK Biobank + All of Us)
  GW2_LOCI34 (optional): 34-locus federated DEL / DUP effects for the locus mirror arrow of Fig. 3d
  GW2_ROOT=... GW2_GRID=... GW2_LOCI=... FIG_OUT=... python figs_dosage_story.py [1] [2] [3]
"""
import os, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle

# Fig. 3 from a federated meta-analysis (2026-11-23): GW2_FIG3_META=<meta_fig3.py output> points the four Fig. 3 inputs at the meta tables
# (same file names and columns as the UK Biobank run folders) and labels the figure with the cohorts (meta_info.json). Unset = UK Biobank runs.
META3 = os.environ.get("GW2_FIG3_META")
if META3:
    os.environ.update(GW2_SIZE2=str(Path(META3) / "size2t_mirror"), GW2_SIZE2VAR="size2t", GW2_ALLCNV=str(Path(META3) / "all2t_mirror"),
                      GW2_TOTDOSE=str(Path(META3) / "totdose"), GW2_GMIRROR=str(Path(META3) / "gene_mirror"))
ROOT = Path(os.environ["GW2_ROOT"]); OUT = Path(os.environ.get("FIG_OUT", ROOT)); OUT.mkdir(parents=True, exist_ok=True)
GRID = Path(os.environ.get("GW2_GRID", ROOT.parent / "class_models_2026-11-22"))
LOCI = Path(os.environ.get("GW2_LOCI", ROOT / "locus_meta"))
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial"], "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 7.6,
                     "axes.linewidth": .7, "xtick.labelsize": 6.6, "ytick.labelsize": 6.6, "legend.fontsize": 6,
                     "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42})
INK = "#1a1a1a"; GREY = "#8a8a8a"; LIGHT = "#dcdcdc"; BLUE = "#1f5c8b"; RED = "#b03a2e"
LCC, HCC = "#4f8fc0", RED                                   # less-constrained, constrained
TL = {"height": "Height", "bmi": "BMI"}; TRAITS = ["height", "bmi"]
fp = lambda p: "n/a" if not np.isfinite(p) else (f"{p:.2g}" if p >= 1e-4 else f"{p:.0e}")


def letter(fig, ax, L, dx=.055, dy=.008):
    p = ax.get_position(); fig.text(p.x0 - dx, p.y1 + dy, L, size=11, weight="bold")


# ------------------------------------------------------------------ panels
SHOW_LOCI = [("16p11.2_BP4_BP5", "16p11.2\nBP4–BP5"), ("16p12.2_BP2_BP3", "16p12.2\nBP2–BP3"), ("1q21.1", "1q21.1\ndistal")]


def loci_panel(ax, trait):
    """Mirror effects of example recurrent loci: fixed-effect meta of adult UK Biobank and All of Us (bars, 95% CI); symbols = each cohort."""
    M = pd.read_csv(LOCI / f"meta_locus_effects_{trait}.tsv", sep="\t"); labs = []
    for i, (nm, lab) in enumerate(SHOW_LOCI):
        for j, (dos, col) in enumerate([("DEL", RED), ("DUP", BLUE)]):
            r = M[(M.NAME == nm) & (M.CNV == dos)].iloc[0]; x = i + (j - .5) * .38
            ax.bar(x, r.beta, .36, color=col, alpha=.85, label=("deletion" if dos == "DEL" else "duplication") if i == 0 else None)
            ax.errorbar(x, r.beta, yerr=1.96 * r.SE, color=INK, lw=.8, capsize=1.5)
            for c, mk in (("UKBB", "o"), ("AoU", "^")):
                if f"{c}_beta" in r and np.isfinite(r[f"{c}_beta"]):
                    ax.scatter(x + (.08 if c == "AoU" else -.08), r[f"{c}_beta"], marker=mk, s=9, facecolor="white", edgecolor=INK, lw=.5, zorder=5)
            s = "***" if r.P < 1e-3 else ("**" if r.P < .01 else ("*" if r.P < .05 else ""))
            ax.text(x, r.beta + np.sign(r.beta) * (1.96 * r.SE + .07), s, ha="center", va="bottom" if r.beta > 0 else "top", size=6)
        g = M[M.NAME == nm].set_index("CNV"); only = " (UKB only)" if g.k_cohorts.min() < 2 else ""
        labs.append(f"{lab}{only}\nn = {int(g.loc['DEL'].n_carriers)} / {int(g.loc['DUP'].n_carriers)}")
    ax.axhline(0, color=INK, lw=.6); ax.set_xticks(range(len(SHOW_LOCI))); ax.set_xticklabels(labs, size=6)
    ax.set_ylim(-1.55, 1.55); ax.set_ylabel(f"{TL[trait]} (SD), CNV effect")


def bins(trait):
    d = ROOT / "class_bins" / trait
    return pd.read_csv(d / "class_dose_bins.tsv", sep="\t"), pd.read_csv(d / "class_dose_shape.tsv", sep="\t")


def curve_panel(ax, trait, part, classes, noloci=True):
    """Dose response of gene classes from one joint model: adjusted difference from people with no gene of the class (95% CI) against the
    mean number of genes deleted (left) or duplicated (right). classes = [(class name, label, colour, marker)]."""
    B, _ = bins(trait); B = B[B.partition == part]
    for k, (c, lab, col, mk) in enumerate(classes):
        for rn in (("full", "no_loci") if noloci else ("full",)):
            d = B[(B.cls == c) & (B.run == rn)]; d = d[np.isfinite(d.beta)].sort_values("order")
            x = np.where(d.side == "DEL", -1, 1) * d.mean_genes.fillna(0) * (1 + (k - (len(classes) - 1) / 2) * .06)
            ax.errorbar(x, d.beta, yerr=1.96 * d.SE, color=col, marker=mk, ms=3.4, lw=1.3 if rn == "full" else .9, ls="-" if rn == "full" else ":",
                        capsize=1.4, mec="white", mew=.4, alpha=1 if rn == "full" else .55, label=lab if rn == "full" else None)
    ax.set_xscale("symlog", linthresh=1); ax.axhline(0, color=LIGHT, lw=.7, zorder=0); ax.axvline(0, color=LIGHT, lw=.7, zorder=0)
    t = [-30, -10, -3, -1, 0, 1, 3, 10, 30]; ax.set_xlim(-45, 45); ax.set_xticks(t); ax.set_xticklabels([str(abs(v)) for v in t])
    ax.set_xlabel("genes deleted          genes duplicated"); ax.set_ylabel(f"{TL[trait]} (SD), adjusted")


def grid(trait, variant):
    f = GRID / trait / variant / "grid_summary.tsv"
    if not f.exists():
        return None
    g = pd.read_csv(f, sep="\t"); g[["term", "arm"]] = g.readout.str.rsplit("|", n=1, expand=True); return g


SIZES = [("N1_1", "1"), ("N2_5", "2–5"), ("N6_20", "6–20"), ("N21_inf", "> 20")]


def size_panel(ax, trait, kind):
    """z against the genome-shift null by the number of genes the CNV call hits. kind 'mirror': mirror z (DUP arm - DEL arm) of the signed
    TWAS count in less-constrained and constrained genes (variant size4s); 'constraint': DEL and DUP z of the other constrained genes."""
    g = grid(trait, "size4s")
    if g is None:
        ax.text(.5, .5, "running", transform=ax.transAxes, ha="center", color=GREY); return
    ser = ([("cls:sig:LC", "mirror", LCC, "less constrained"), ("cls:sig:HC", "mirror", HCC, "constrained")] if kind == "mirror"
           else [("cls:other:HC", "DEL", RED, "deletions"), ("cls:other:HC", "DUP", BLUE, "duplications")])
    w = .38
    for j, (m, arm, col, lab) in enumerate(ser):
        for i, (c, _) in enumerate(SIZES):
            r = g[(g.term == f"{c}~{m}") & (g.arm == arm)]
            if not len(r):
                continue
            r = r.iloc[0]; x = i + (j - .5) * w
            ax.bar(x, r.z_emp, w * .95, color=col, alpha=.9, label=lab if i == 0 else None)
            if np.isfinite(r.noloci_z_emp):
                ax.plot([x - w * .4, x + w * .4], [r.noloci_z_emp] * 2, color=INK, lw=1.1, zorder=5)
            s = "*" if r.P_emp < .05 else ""
            if s and r.loo_held == r.loo_n:
                s = "**"
            if s:
                ax.text(x, r.z_emp + np.sign(r.z_emp) * .25, s, ha="center", va="center", size=6.5)
    ax.axhline(0, color=INK, lw=.6)
    for v in (-1.96, 1.96):
        ax.axhline(v, color=GREY, lw=.4, ls=":")
    ax.set_xticks(range(len(SIZES))); ax.set_xticklabels([s for _, s in SIZES]); ax.set_xlabel("genes hit by the CNV call")
    ax.set_ylabel("mirror z (DUP − DEL)" if kind == "mirror" else "constrained-gene z"); ax.text(1, -.32, "* P < 0.05; ** and held without each major locus;\nbar: z; dash: 34 recurrent intervals removed" if kind == "constraint" and trait == "bmi" else "", transform=ax.transAxes, ha="right", va="top", size=5)


# ------------------------------------------------------------------ Figure 1
def fig1():
    fig = plt.figure(figsize=(7.2, 4.9)); gs = fig.add_gridspec(2, 2, left=.085, right=.985, top=.93, bottom=.1, hspace=.72, wspace=.28)
    for j, t in enumerate(TRAITS):
        ax = fig.add_subplot(gs[0, j]); loci_panel(ax, t); ax.set_title(f"{TL[t]}: single recurrent CNVs are mirrors", loc="left", weight="bold")
        letter(fig, ax, "ab"[j])
        if j == 0:
            ax.legend(frameon=False, loc="upper left")
        else:
            ax.scatter([], [], marker="o", s=9, facecolor="white", edgecolor=INK, lw=.5, label="UK Biobank")
            ax.scatter([], [], marker="^", s=9, facecolor="white", edgecolor=INK, lw=.5, label="All of Us"); ax.legend(frameon=False, loc="upper right")
        ax = fig.add_subplot(gs[1, j]); curve_panel(ax, t, "all", [("all", "all genes", INK, "o")])
        ax.set_title(f"{TL[t]}: genome-wide burden is U-shaped", loc="left", weight="bold"); letter(fig, ax, "cd"[j])
        h, l = ax.get_legend_handles_labels(); h.append(Line2D([], [], color=GREY, ls=":", lw=.9)); l.append("34 recurrent intervals removed")
        ax.legend(h, l, frameon=False, loc="lower center" if t == "height" else "upper center")
    save(fig, "Fig1_paradox")


# ------------------------------------------------------------------ Figure 2
def fig2():
    con = lambda tw: [(f"{tw}_LC", "less constrained (s$_{het}$ ≤ 0.029)", LCC, "o"), (f"{tw}_HC", "constrained (s$_{het}$ > 0.029)", HCC, "s")]
    rows = [("other", "twas_x_con", "weak TWAS effect (|z| ≤ 3 or none)"), ("up", "twas_x_con", "TWAS z > +3 (higher expression, higher trait)"),
            ("down", "twas_x_con", "TWAS z < −3 (higher expression, lower trait)"), ("abs3", "abs3_x_con", "|TWAS z| > 3, either direction")]
    fig = plt.figure(figsize=(7.2, 11.2))
    gs = fig.add_gridspec(5, 4, left=.085, right=.985, top=.965, bottom=.045, hspace=.62, wspace=.75, height_ratios=[1, 1, 1, 1, .9])
    k = 0
    for i, (tw, part, ttl) in enumerate(rows):
        for j, t in enumerate(TRAITS):
            ax = fig.add_subplot(gs[i, 2 * j:2 * j + 2]); curve_panel(ax, t, part, con(tw))
            ax.set_title(f"{TL[t]}: {ttl}", loc="left", weight="bold", size=7.4); letter(fig, ax, "abcdefgh"[k]); k += 1
            h, l = ax.get_legend_handles_labels(); h.append(Line2D([], [], color=GREY, ls=":", lw=.9)); l.append("34 recurrent intervals removed")
            ax.legend(h, l, frameon=False, loc="lower center" if (tw == "other" and t == "height") else "upper left", fontsize=5.4)
    for j, (t, kind) in enumerate([("height", "mirror"), ("height", "constraint"), ("bmi", "mirror"), ("bmi", "constraint")]):
        ax = fig.add_subplot(gs[4, j]); size_panel(ax, t, kind)
        ax.set_title(f"{TL[t]}: {'TWAS mirror' if kind == 'mirror' else 'constraint'}\nby CNV size", loc="left", weight="bold", size=7.2)
        letter(fig, ax, "ijkl"[j], dx=.06); ax.legend(frameon=False, fontsize=5.2, loc="best")
    save(fig, "Fig2_dissection")


# ------------------------------------------------------------------ Figure 3
CL = {"cls:up:LC": "TWAS z > +3, less constrained", "cls:up:HC": "TWAS z > +3, constrained", "cls:down:LC": "TWAS z < −3, less constrained",
      "cls:down:HC": "TWAS z < −3, constrained", "cls:sig:LC": "signed TWAS count, less constrained", "cls:sig:HC": "signed TWAS count, constrained",
      "cls:sig:all": "signed TWAS count, all genes", "cls:abs:LC": "|z| > 3 count, less constrained", "cls:abs:HC": "|z| > 3 count, constrained",
      "cls:other:LC": "other genes, less constrained", "cls:other:HC": "other genes, constrained", "none_H": "all constrained genes",
      "n_all": "total gene count"}
UNI = [("u_up_LC", "cls:up:LC"), ("u_up_HC", "cls:up:HC"), ("u_down_LC", "cls:down:LC"), ("u_down_HC", "cls:down:HC"), ("u_sig_LC", "cls:sig:LC"),
       ("u_sig_HC", "cls:sig:HC"), ("u_sig_all", "cls:sig:all"), ("u_abs_LC", "cls:abs:LC"), ("u_abs_HC", "cls:abs:HC"),
       ("u_other_LC", "cls:other:LC"), ("u_other_HC", "cls:other:HC"), ("u_none_H", "none_H")]   # total gene count alone is not testable: the shift null leaves it unchanged
J6 = ["cls:up:LC", "cls:up:HC", "cls:down:LC", "cls:down:HC", "cls:other:LC", "cls:other:HC"]
S4 = ["cls:up:LC", "cls:up:HC", "cls:down:LC", "cls:down:HC", "cls:other:HC"]
ARMS = ["DEL", "DUP", "mirror"]


def cell(g, term, arm):
    if g is None:
        return None
    r = g[(g.term == term) & (g.arm == arm)]
    return r.iloc[0] if len(r) else None


def heat(ax, rows, title, ylabels=True):
    """rows = [(label, {trait: grid df}, term)]; columns = trait x (DEL, DUP, mirror) for all CNVs, then the same with the 34 recurrent
    intervals removed. Fill = z against the genome-shift null (capped at +-5); * P < 0.05, ** BH q < 0.05 within the block;
    box = P < 0.05 and held (P < 0.05, same sign) when each major recurrent locus is removed."""
    cols = [(t, a, blk) for blk in ("all", "noloci") for t in TRAITS for a in ARMS]
    Zm = np.full((len(rows), len(cols)), np.nan); P = np.full(Zm.shape, np.nan); H = np.zeros(Zm.shape, bool)
    for i, (_, G, term) in enumerate(rows):
        for j, (t, a, blk) in enumerate(cols):
            r = cell(G.get(t), term, a)
            if r is None:
                continue
            Zm[i, j] = r.z_emp if blk == "all" else r.noloci_z_emp; P[i, j] = r.P_emp if blk == "all" else r.noloci_P_emp
            H[i, j] = blk == "all" and r.P_emp < .05 and r.loo_n > 0 and r.loo_held == r.loo_n
    ax.imshow(np.clip(Zm, -5, 5), cmap="RdBu_r", vmin=-5, vmax=5, aspect="auto")
    nb = len(TRAITS) * len(ARMS)
    for blk in (0, 1):
        pv = P[:, blk * nb:(blk + 1) * nb]; ok = np.isfinite(pv); q = np.full(pv.shape, np.nan)
        if ok.any():
            v = pv[ok]; o = np.argsort(v); qq = v[o] * len(v) / np.arange(1, len(v) + 1); qq = np.minimum.accumulate(qq[::-1])[::-1]
            tmp = np.empty(len(v)); tmp[o] = qq; q[ok] = tmp
        for i in range(pv.shape[0]):
            for j in range(nb):
                if ok[i, j] and pv[i, j] < .05:
                    ax.text(blk * nb + j, i, "**" if q[i, j] < .05 else "*", ha="center", va="center", size=6.5, color=INK)
    for i, j in zip(*np.where(H)):
        ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, fill=False, ec=INK, lw=1.1))
    for i, j in zip(*np.where(~np.isfinite(Zm))):
        ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, fc="#eeeeee", ec="white", lw=.5))
    for x in (2.5, 8.5):
        ax.axvline(x, color="white", lw=1.5)
    ax.axvline(5.5, color=INK, lw=1)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([a for _, a, _ in cols], rotation=60, ha="right", size=5.8)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows] if ylabels else [], size=5.9)
    for k, lab in enumerate(["Height", "BMI", "Height", "BMI"]):
        ax.text(1 + 3 * k, -.75, lab, ha="center", va="bottom", size=6.4, weight="bold")
    ax.text(2.5, -1.75, "all CNVs", ha="center", va="bottom", size=6.4); ax.text(8.5, -1.75, "34 recurrent intervals removed", ha="center", va="bottom", size=6.4)
    ax.set_title(title, loc="left", weight="bold", size=7.4, pad=24)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)


def locus_mirror(trait):
    """Correlation of deletion and duplication effects across the recurrent reciprocal loci (federated meta-analysis; GW2_LOCI34)."""
    base = os.environ.get("GW2_LOCI34")
    if not base:
        return None
    f = Path(base) / ("height/inputs/sacks_height_all_loci_exact.csv" if trait == "height" else "bmi/inputs/sacks_bmi_complete_reciprocal_loci_exact.csv")
    if not f.exists():
        return None
    d = pd.read_csv(f); dc = [c for c in d.columns if c.lower().startswith("deletion_beta")][0]; uc = [c for c in d.columns if c.lower().startswith("duplication_beta")][0]
    x, y = d[dc].to_numpy(float), d[uc].to_numpy(float); ok = np.isfinite(x) & np.isfinite(y); x, y = x[ok], y[ok]
    r = np.corrcoef(x, y)[0, 1]; n = len(x); z = np.arctanh(r) * np.sqrt(n - 3)
    return dict(r=r, n=n, z=z, sign=np.sign(np.mean(y) - np.mean(x)))


GM = Path(os.environ.get("GW2_GMIRROR", ROOT.parent / "gene_mirror_2026-11-22"))   # gene_mirror_qt.py output, one folder per trait
GMD = Path(os.environ.get("GW2_GMDIAG", ROOT.parent / "gene_mirror_twas_diag_2026-11-22"))   # gene_mirror_twas_diag.py output


def shanta_arrow(ax, v, col, sig, thick, L, off=(0, 0), z=3):
    """Arrow drawn as in Shanta et al. Fig. 3I (fig2_shanta.R epanel): butt-ended shaft that stops at the head; significant = one filled sharp
    triangle head, not significant = thin shaft with an open two-line chevron. Jon 2026-11-22: "Use the arrow formats from Shanta. If cnv-twas
    convergence is significant, the orange arrow should be just as thick as the yellow arrow. Thin arrows are for trends that are not
    significant." thick: 'full' (yellow, and significant orange), 'arm' (green constraint). Lengths below a quarter of the half-width L are
    drawn at that minimum (as Shanta, EMIN = 0.25); longer than L are capped."""
    from matplotlib.patches import Polygon
    v = np.asarray(v, float); n = np.linalg.norm(v)
    if not np.isfinite(n) or n == 0:
        return
    u = v / n; n = min(max(n, .25 * L), L - .3); o = np.asarray(off, float); tip = o + u * n; w = np.array([-u[1], u[0]])
    if sig:
        hl, hw, lw = .10 * L, .05 * L, 2.6     # Jon 2026-11-22: "Make all significant arrows the same thickness as the green arrow"
    else:
        hl, hw, lw = .07 * L, .04 * L, .8
    base = tip - u * hl
    if sig:
        ax.plot(*np.c_[o, base + u * .01], color=col, lw=lw, solid_capstyle="butt", zorder=z)
        ax.add_patch(Polygon([tip, base + w * hw, base - w * hw], closed=True, fc=col, ec="none", zorder=z))
    else:
        ax.plot(*np.c_[o, tip], color=col, lw=lw, solid_capstyle="butt", zorder=z)
        ax.plot(*np.c_[base + w * hw, tip, base - w * hw], color=col, lw=.8, solid_joinstyle="miter", zorder=z)


def schematic(ax, trait, caption=True):
    """Components of the dose response in the plane of deletion and duplication effects (after Shanta et al. Fig. 3I). Arrow length = |z|.
    yellow: mirror of the recurrent loci (DEL-DUP correlation across loci, Fisher z); orange: TWAS mirror, mirror z of the signed TWAS count
    in less-constrained (light) and constrained (dark) genes, drawn along the mirror line; green: constraint, DEL and DUP z of the other
    constrained genes (joint model). Filled = P < 0.05; open, thin = not significant."""
    L = 7.5; ax.plot([-L, L], [L, -L], color=GREY, lw=.7, ls="--", zorder=1); ax.axhline(0, color=INK, lw=.9); ax.axvline(0, color=INK, lw=.9)
    j6 = grid(trait, "joint6"); lm = locus_mirror(trait)
    sgn = lm["sign"] if lm else (1 if trait == "height" else -1)          # +1: duplication raises / deletion lowers the trait
    d = np.array([-1, 1]) / np.sqrt(2) * sgn

    def arrow(v, col, sig, lw, off=(0, 0), z=3):
        v = np.asarray(v, float); n = np.linalg.norm(v)
        if not np.isfinite(n) or n == 0:
            return
        v = v / n * min(max(n, .8), L - .4); o = np.asarray(off, float)
        ax.add_patch(FancyArrowPatch(o, o + v, arrowstyle="-|>" if sig else "->", mutation_scale=7 + lw * 1.6, lw=lw if sig else .9,
                                     color=col, zorder=z, shrinkA=0, shrinkB=0))
    gmf = GM / trait / "gene_mirror_summary.tsv"
    Q = []   # Jon 2026-11-23: "always draw the arrows aligned perfectly with the axis (Del, dup, mirror) with shorter arrow on top": no offsets;
    #          arrows are queued and drawn longest first, so every shorter arrow lies on top of the longer ones on the same axis
    if gmf.exists():
        # Jon 2026-11-22 ("as per the Shanta paper"): yellow = TOTAL mirror, gene-level r(zDEL, zDUP) (all genes, all CNV sizes), length |z| vs
        # the gene-circle shift null. A correlation has no direction along its line, so it is drawn toward 'duplication raises' (DUP = -DEL
        # when r < 0), the same frame as the TWAS arrows.
        # Orange = CNV-TWAS convergence from the four tiles, on FIXED axes (Jon 2026-11-22: "Do you really think it's a good idea to set the
        # direction of the arrows, the direction probably has a pretty wide confidence interval. It might be better just to stick with one mirror
        # direction and then the vertical and horizontal axis"): the four-tile mirror M4 along DUP = -DEL (positive = up-left), the deletion arm
        # (t_DEL up - t_DEL down)/sqrt2 on the horizontal axis and the duplication arm (t_DUP up - t_DUP down)/sqrt2 on the vertical axis, each z vs
        # the genome-shift null (raw sign: positive = deleting / duplicating genes predicted to raise the trait raises it). Light = less-constrained,
        # dark = constrained genes; filled = all CNVs, outlined = single-gene CNVs; thick = P < 0.05. Offsets keep parallel arrows apart.
        g = pd.read_csv(gmf, sep="\t"); g = g[(g["size"] == "all sizes") & (g.stratum == "all")].set_index("stat")
        tot = g.loc["total"]
        dd = np.array([-1, 1]) / np.sqrt(2) if tot.obs < 0 else np.array([1, 1]) / np.sqrt(2)
        if S2VAR == "size2t":
            total_arrows(Q, trait, dose=True)   # Jon 2026-11-23 'apply this change to both versions': yellow = total on both diagonals
            ax.plot([-L, L], [-L, L], color=GREY, lw=.7, ls="--", zorder=1)
        else:
            Q.append((dd * abs(tot.z), "#f2c230", tot.P < .05, "full"))
        # orange components: significant = as thick as the yellow arrow; trends (1.5 <= |z| < 1.96) thin with an open head; |z| < 1.5 not drawn.
        # Mirror on top of the yellow arrow; arms beside the green constraint arm of the same axis.
        if _s3sum() is not None:
            # Jon 2026-11-23: "Height has a significant twas mirror for single genes. I think we can add a significant orange arrow to the height
            # mirror axis and I don't think we need to distinguish single gene from >=2 gene." Orange components now come from the same model as
            # Fig. 3a (size3, matched TWAS, per constraint stratum): for each statistic the CNV-size row (single-gene or >= 2-gene) with the larger
            # |z| is drawn, with no size marking; significance uses P x 2 (Bonferroni over the two size rows) so that taking the larger of two
            # is not anticonservative. Thin open-headed = trend (|z| >= 1.5, P x 2 >= 0.05); |z| < 1.5 not drawn.
            mdir = np.array([-1, 1]) / np.sqrt(2)
            for st, col in (("LC", "#f28e2b"), ("HC", "#7b3294")):
                # arrow rule (pick_arrow) + one arrow per colour (component): mirror M4 on the diagonal, else the DEL / DUP arms
                g_ = lambda stat: (lambda t, s, st_: s2get(t, s, st_, stat))
                component(Q, trait, (st, "HC" if st == "LC" else "LC"), col, "full", (g_("M4"), mdir, ("M4", st)),
                          [(g_("armDEL"), (1, 0)), (g_("armDUP"), (0, 1))])
            if S2VAR == "size2t":   # grey: total-genes term of the joint model (ordinary gene), Jon 2026-11-23
                o_ = lambda arm: (lambda t, s, st_: s2term(t, f"{s}~n_all", arm))
                pass   # Jon 2026-11-23: "Keep the term in the model and in the figure. Just remove the gray arrows"
        j6 = xgrid(trait, trait, "joint6") if xgrid(trait, trait, "joint6") is not None else j6
    elif (XG / "mirror4_cross_summary.tsv").exists():
        if lm:
            arrow(d * abs(lm["z"]), "#f0c030", abs(lm["z"]) > 1.96, 7, z=2)
        # Jon 2026-11-22: TWAS mirror = four-tile mirror (DEL and DUP of TWAS z > 3 and z < -3 genes) of the trait on its own TWAS;
        # solid = all CNVs, outlined = single-gene CNVs; light = less-constrained, dark = constrained genes
        X = xsum(); X = X[(X.trait == f"{trait}_x_{trait}") & (X.stat == "M4")]
        for st, col, side in (("LC", "#f28e2b", 1), ("HC", "#7b3294", -1)):
            for size, solid in (("all sizes", True), ("N1_1", False)):
                r = X[(X.stratum == st) & (X["size"] == size)]
                if not len(r):
                    continue
                r = r.iloc[0]; perp = np.array([1, 1]) / np.sqrt(2) * side * (.25 if solid else .6)
                if solid:
                    arrow(d * r.z, col, r.P < .05, 3.2, off=perp, z=4)
                else:
                    v = d * r.z; n = np.linalg.norm(v); v = v / n * min(max(n, .8), L - .4)
                    ax.add_patch(FancyArrowPatch(perp, perp + v, arrowstyle="-|>", mutation_scale=9, lw=1.0, ec=col, fc="white", ls="-", zorder=4,
                                                 shrinkA=0, shrinkB=0))
        j6 = xgrid(trait, trait, "joint6") if xgrid(trait, trait, "joint6") is not None else j6
    else:
        for var, term, col, off in (("u_sig_LC", "cls:sig:LC", "#f28e2b", .0), ("u_sig_HC", "cls:sig:HC", "#7b3294", .0)):
            r = cell(grid(trait, var), term, "mirror")
            if r is not None and np.isfinite(r.z_emp):
                perp = np.array([1, 1]) / np.sqrt(2) * (.22 if "LC" in term else -.22)
                arrow(d * r.z_emp, col, r.P_emp < .05, 3.2, off=perp, z=4)       # negative z points against the locus mirror
    if _s3sum() is not None:
        # green, one arrow (component): the mirror of the 'weak or no TWAS' HC term ((t_DUP - t_DEL)/sqrt2, readout '|mirror') on the
        # diagonal, else its DEL / DUP terms on the axes
        gcls = "other" if S2VAR == "size2t" else "none"
        g_ = lambda arm: (lambda t, s, st_: s2term(t, f"{s}~cls:{gcls}:{st_}", arm))
        component(Q, trait, ("HC",), "#1b7837", "arm", (g_("mirror"), np.array([-1, 1]) / np.sqrt(2), ("mirOth", "HC")),
                  [(g_("DEL"), (1, 0)), (g_("DUP"), (0, 1))])
    for arm, vec in (("DEL", (1, 0)), ("DUP", (0, 1))) if _s3sum() is None else ():   # older runs only (no size2t / size3 results)
        r = cell(j6, "cls:other:HC", arm)
        if r is not None and np.isfinite(r.z_emp):
            shanta_arrow(ax, np.array(vec) * r.z_emp, "#1b7837", r.P_emp < .05, "arm", L, z=5)
    flush_arrows(ax, Q, L)
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.text(L, -.5, f"deletion raises {TL[trait] if trait == 'bmi' else 'height'}", ha="right", va="top", size=5.6); ax.text(-L, -.5, "lowers", ha="left", va="top", size=5.6)
    ax.text(.35, L, f"duplication raises {TL[trait] if trait == 'bmi' else 'height'}", ha="left", va="top", size=5.6); ax.text(.35, -L, "lowers", ha="left", va="bottom", size=5.6)
    ax.set_title(TL[trait], loc="left", weight="bold", size=9)
    if gmf.exists() and caption:
        sh_txt = "n.a."
        fd = GMD / trait / "twas_diag.tsv"
        if fd.exists():
            dg = pd.read_csv(fd, sep="\t").set_index("analysis")
            if "explained B4" in dg.index:
                e4 = dg.loc["explained B4"]; sh_txt = f"{100 * max(e4.share, 0):.0f}%" + (f" (z = {e4.z:.1f})" if e4.P < .05 else " (n.s.)")
        ax.text(-L, -L - .5, f"total mirror r = {tot.obs:.2f}, z = {tot.z:.1f} ({int(tot.loo_held)}/{int(tot.n_loo)})\n"
                f"without the 34 loci z = {tot.noloci_z:.1f}\nshare explained by the four TWAS tiles: {sh_txt}", ha="left", va="top", size=5.2, color="#7a5a00")


def figS_heatmaps():
    G = {v: {t: grid(t, v) for t in TRAITS} for v in [u for u, _ in UNI] + ["joint6", "size4"]}
    fig = plt.figure(figsize=(7.2, 10.6))
    gs = fig.add_gridspec(3, 2, left=.3, right=.985, top=.94, bottom=.04, hspace=.42, wspace=.9, height_ratios=[13, 20, 13], width_ratios=[1, 1])
    axA = fig.add_axes([.3, .745, .3, .185]); heat(axA, [(CL[m], G[v], m) for v, m in UNI], "Each term alone (+ total gene count)")
    axB = fig.add_axes([.69, .835, .295, .095]); heat(axB, [(CL[m].replace("less constrained", "LC").replace("constrained", "HC"), G["joint6"], m) for m in J6], "Joint model: six gene classes")
    axC = fig.add_axes([.3, .345, .3, .33])
    heat(axC, [(f"{CL[m]}, {lab} gene{'s' if lab != '1' else ''}", G["size4"], f"{c}~{m}") for m in S4 for c, lab in SIZES], "Joint model by CNV size (genes hit by the call)")
    for y in (3.5, 7.5, 11.5, 15.5):
        axC.axhline(y, color="white", lw=1.5)
    cax = fig.add_axes([.69, .79, .2, .008]); sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(-5, 5))
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal"); cb.set_label("z against genome-shift null", size=6); cb.ax.tick_params(labelsize=5.6)
    fig.text(.69, .748, "* P < 0.05;  ** BH q < 0.05;\nbox: P < 0.05 and held without\neach major recurrent locus", size=5.6, color=INK, va="top")
    axD = [fig.add_axes([.14 + .43 * j, .045, .36, .245]) for j in range(2)]
    for ax, t in zip(axD, TRAITS):
        schematic(ax, t)
    hs = [Line2D([], [], color="#f0c030", lw=6), Line2D([], [], color="#f28e2b", lw=3), Line2D([], [], color="#7b3294", lw=3), Line2D([], [], color="#1b7837", lw=3)]
    fig.legend(hs, ["mirror of recurrent loci (DEL–DUP correlation)", "TWAS mirror, less-constrained genes", "TWAS mirror, constrained genes",
                    "constraint (other constrained genes)"], loc="lower center", ncol=2, frameon=False, fontsize=6, bbox_to_anchor=(.55, 0))
    for ax, Lt, dx in ((axA, "a", .28), (axB, "b", .05), (axC, "c", .28), (axD[0], "d", .04)):
        p = ax.get_position(); fig.text(p.x0 - dx, p.y1 + (.028 if ax is not axD[0] else .005), Lt, size=11, weight="bold")
    save(fig, "FigS_class_heatmaps")   # Jon 2026-11-22: superseded in Fig. 3 by the Shanta-style mirror and dosage grids; kept as a supplementary sheet


# ------------------------------------------------------------------ Figure 3 (Shanta et al. Fig. 3E-F style, Jon 2026-11-22)
XG = Path(os.environ.get("GW2_XGRID", ROOT.parent / "cross_mirror_2026-11-22"))
SZ2 = [("all sizes", "joint6", "all CNVs"), ("N1_1", "size4", "single-gene CNVs")]


def xsum():
    return pd.read_csv(XG / "mirror4_cross_summary.tsv", sep="\t")


def xgrid(cnv, twas, variant):
    f = XG / "grid" / f"{cnv}_x_{twas}" / variant / "grid_summary.tsv"
    if not f.exists():
        return None
    g = pd.read_csv(f, sep="\t"); g[["term", "arm"]] = g.readout.str.rsplit("|", n=1, expand=True); return g


def tile_grid(ax, M, rows, cols, title, sub=None, col_groups=None, row_colors=None):
    """M: rows x cols array of (z, P, held_all) or None. Fill = z (capped +-4, as Shanta Fig. 3E-F); * P < 0.05; ** BH q < 0.05 over the
    panel; box = P < 0.05 and held when each major recurrent locus is removed."""
    Z = np.array([[np.nan if c is None else c[0] for c in r] for r in M]); P = np.array([[np.nan if c is None else c[1] for c in r] for r in M])
    ax.imshow(np.clip(Z, -4, 4), cmap="RdBu_r", vmin=-4, vmax=4, aspect="equal")
    ok = np.isfinite(P); q = np.full(P.shape, np.nan)
    if ok.any():
        v = P[ok]; o = np.argsort(v); qq = np.minimum.accumulate((v[o] * len(v) / np.arange(1, len(v) + 1))[::-1])[::-1]; t = np.empty(len(v)); t[o] = qq; q[ok] = t
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            c = M[i][j]
            if c is None:
                ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, fc="#eeeeee", ec="white")); continue
            if P[i, j] < .05:
                ax.text(j, i + .08, "**" if q[i, j] < .05 else "*", ha="center", va="center", size=8)
            if c[2]:
                ax.add_patch(Rectangle((j - .48, i - .48), .96, .96, fill=False, ec=INK, lw=1.4))
    for k in range(Z.shape[1] + 1):
        ax.axvline(k - .5, color="white", lw=1)
    for k in range(Z.shape[0] + 1):
        ax.axhline(k - .5, color="white", lw=1)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=0, ha="center", size=6.2)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, size=7)
    if row_colors:
        for tl, c in zip(ax.get_yticklabels(), row_colors):
            tl.set_color(c)
    if col_groups:
        for x0, x1, lab in col_groups:
            ax.text((x0 + x1) / 2, -.75, lab, ha="center", va="bottom", size=6.6, weight="bold")
            if x1 < len(cols) - 1:
                ax.axvline(x1 + .5, color=INK, lw=1.2)
    ax.set_title(title + (f"\n{sub}" if sub else ""), weight="bold", size=7.8, pad=22 if col_groups else 6)
    for s_ in ax.spines.values():
        s_.set_visible(False)
    ax.tick_params(length=0)


def figS_cross_mirror():
    """Jon 2026-11-22: 'see the four tiles separately for Mirror, and then the final estimate ... Row 1 +3: Dup up, Del up; Row 2 -3: dup
    down, Del down; next to this 2x2 grid one merged tile (same size as the 2x2 grid) that gives the mirror estimate; summary Del estimate
    (Del up - Del down) etc. find a place'. a: one unit per CNV trait x CNV size (rows) and TWAS trait x constraint stratum (columns):
    2x2 tiles (z of each signed-count term in the joint model), a strip under each column with the dosage arm (up - down)/sqrt2, and the
    merged four-tile mirror M4 = (-DEL up + DUP up + DEL down - DUP down)/2. b: constraint (other HC genes, |TWAS z| < 3). c: schematic.
    null model: genome shift of the gene annotation, shared shifts; every statistic recomputed in each replicate (mirror4_summary.py)."""
    X = xsum(); cmap = plt.cm.RdBu_r; norm = plt.Normalize(-4, 4)
    def mget(cnv, twas, size, st, stat):
        r = X[(X.trait == f"{cnv}_x_{twas}") & (X["size"] == size) & (X.stratum == st) & (X.stat == stat)]
        if not len(r):
            return None
        r = r.iloc[0]; return (r.z, r.P, r.P < .05 and r.loo_n > 0 and r.loo_held == r.loo_n)
    def tget(g, term, arm):
        if g is None:
            return None
        rr = g[(g.term == term) & (g.arm == arm)]
        if not len(rr):
            return None
        rr = rr.iloc[0]; return (rr.z_emp, rr.P_emp, rr.P_emp < .05 and rr.loo_n > 0 and rr.loo_held == rr.loo_n)
    def tile(ax, x, y, w, h, c, fs=5.6):
        if c is None or not np.isfinite(c[0]):
            ax.add_patch(Rectangle((x, y), w, h, fc="#eeeeee", ec="white", lw=1)); return
        ax.add_patch(Rectangle((x, y), w, h, fc=cmap(norm(np.clip(c[0], -4, 4))), ec="white", lw=1))
        dark = abs(c[0]) > 2.6
        ax.text(x + w / 2, y + h / 2, f"{c[0]:.1f}" + ("*" if c[1] < .05 else ""), ha="center", va="center", size=fs,
                color="white" if dark else INK, weight="bold" if c[1] < .05 else "normal")
        if c[2]:
            ax.add_patch(Rectangle((x + .04, y + .04), w - .08, h - .08, fill=False, ec=INK, lw=1.3))
    fig = plt.figure(figsize=(7.2, 5.4))
    ax = fig.add_axes([.02, .11, .96, .80]); ax.set_aspect("equal"); ax.axis("off")
    UW, GAP, SH, CG, RG = 4.35, .35, .55, 1.0, 1.15        # unit width, tile-to-mirror gap, arm strip height, column gap, row gap
    UH = 2 + .12 + SH
    cols = [(t, s) for t in TRAITS for s in ("LC", "HC")]
    rows = [(c, sz, var, lab) for c in TRAITS for sz, var, lab in SZ2]
    ytop = 0; TS = Tiles()   # FDR family = one statistic across the 16 CNV trait x size x TWAS trait x LC/HC units
    for i, (cnv, size, var, slab) in enumerate(rows):
        y0 = -(i * (UH + RG)) - 2           # bottom of the 2x2
        ax.text(-1.25, y0 + 1, f"{TL[cnv]} CNVs\n{slab}", ha="right", va="center", size=7, weight="bold")
        for j, (twas, st) in enumerate(cols):
            x0 = j * (UW + CG); g = xgrid(cnv, twas, var); pre = "N1_1~" if size == "N1_1" else ""
            for r_, (band, cls) in enumerate((("z > +3", "up"), ("z < −3", "down"))):
                for k, dos in enumerate(("DUP", "DEL")):
                    TS.add(ax, x0 + k, y0 + 1 - r_, 1, 1, tget(g, f"{pre}cls:{cls}:{st}", dos), (cls, dos))
            for k, arm in enumerate(("armDUP", "armDEL")):
                TS.add(ax, x0 + k, y0 - .12 - SH, 1, SH, mget(cnv, twas, size, st, arm), (arm,), fs=5)
            TS.add(ax, x0 + 2 + GAP, y0, 2, 2, mget(cnv, twas, size, st, "M4"), ("M4",), fs=8)
            if i == 0:
                ax.text(x0 + (2 + GAP + 2) / 2, y0 + 2.75, f"{TL[twas]} TWAS genes\n{'less constrained' if st == 'LC' else 'constrained'}",
                        ha="center", va="bottom", size=6.8, weight="bold")
                ax.text(x0 + .5, y0 + 2.08, "DUP", ha="center", va="bottom", size=6, color=BLUE)
                ax.text(x0 + 1.5, y0 + 2.08, "DEL", ha="center", va="bottom", size=6, color=RED)
                ax.text(x0 + 3 + GAP, y0 + 2.08, "mirror", ha="center", va="bottom", size=6)
            if j == 0:
                ax.text(x0 - .08, y0 + 1.5, "z>+3", ha="right", va="center", size=5.4)
                ax.text(x0 - .08, y0 + .5, "z<−3", ha="right", va="center", size=5.4)
                ax.text(x0 - .08, y0 - .12 - SH / 2, "up−down", ha="right", va="center", size=5.4)
        if i == 1:
            ax.axhline(y0 - .12 - SH - RG / 2, color=GREY, lw=.6, xmin=0, xmax=1)
    TS.draw()
    ax.set_xlim(-4.6, 4 * UW + 3 * CG + .1)
    ax.set_ylim(-(3 * (UH + RG)) - 2 - .12 - SH - .2, 3.9)
    cax = fig.add_axes([.12, .075, .2, .012]); sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal"); cb.set_label("z against genome-shift null", size=6); cb.ax.tick_params(labelsize=5.6)
    fig.text(.36, .08, "tile = z of the signed-count term (TWAS z > +3 or < −3 genes, joint model);  up−down = (up − down)/√2 per dosage;\n"
             "mirror = (−DEL up + DUP up + DEL down − DUP down)/2;  * P < 0.05;  ** FDR < 5% (BH within each statistic);\n"
             "box: P < 0.05 and held when each major recurrent locus is removed",
             size=5.4, va="center")
    fig.text(.55, .95, "CNV–TWAS mirror, every CNV trait by every TWAS trait (all CNVs and single-gene CNVs)", ha="center", size=8.5, weight="bold")
    save(fig, "FigS_cross_trait_mirror")   # Jon 2026-11-22: "You can put the all CNVs by all TWAS in the supplementary"


# ------------------------------------------------------------------ Figure 3, matched TWAS only, CNVs split 1 gene vs >= 2 genes (Jon 2026-11-22 mockup)
# Jon 2026-11-22: Fig. 3 and the supplementary full grid from the size3 run (conga_grid_qt.py, per CNV size 1 gene / >= 2 genes, no total-gene
# term): size3c = TWAS up (z > +3) / down (z < -3) / weak (has a z, -3 <= z <= 3) / none (no TWAS z) x LC / HC; size3a = the same, all genes.
# 'weak' replaces 'other' so that 'genes with weak TWAS signal, -3 < z < +3' excludes genes without a TWAS z. Earlier runs size2 / size2j
# (where 'other' included genes without a z) are sensitivity runs.
S2VAR = os.environ.get("GW2_SIZE2VAR", "size3c")
XS2 = Path(os.environ.get("GW2_SIZE2", ROOT.parent / "size3_mirror_2026-11-22"))
WEAK = "weak"
S2ROWS = [("N2_inf", "≥2-gene CNVs"), ("N1_1", "Single-gene CNVs")]
CMAP3 = plt.cm.RdBu_r; NORM3 = plt.Normalize(-4, 4)


from functools import lru_cache


@lru_cache(None)
def _s3sum():
    f = XS2 / "mirror4_cross_summary.tsv"
    return pd.read_csv(f, sep="\t") if f.exists() else None


@lru_cache(None)
def _s3grid(cnv, twas, variant):
    f = XS2 / "grid" / f"{cnv}_x_{twas}" / variant / "grid_summary.tsv"
    if not f.exists():
        return None
    g = pd.read_csv(f, sep="\t"); g[["term", "arm"]] = g.readout.str.rsplit("|", n=1, expand=True); return g


def g3get(cnv, twas, variant, size, st, stat):
    X = _s3sum()
    if X is None:
        return None
    r = X[(X.trait == f"{cnv}_x_{twas}") & (X.variant == variant) & (X["size"] == size) & (X.stratum == st) & (X.stat == stat)]
    if not len(r):
        return None
    r = r.iloc[0]; return (r.z, r.P, r.P < .05 and r.loo_n > 0 and r.loo_held == r.loo_n)


@lru_cache(None)
def _obs(d):
    f = Path(d) / "observed_terms.tsv"
    return pd.read_csv(f, sep="\t") if f.exists() else None


def analytic_held(d, readout, t_full):
    """Locus robustness of an analytic (total-genes) term (Jon 2026-11-23 'confirm that all tiles have been tested by LOO'): |t| >= 1.96 with the
    sign of the full fit in every leave-one-major-locus-out run (the 8 major loci of run_info.json) and in minor_all (all minor loci removed),
    the same 9 runs as the mirror / dose-aligned summaries. null model: analytic t of the joint model in each run."""
    import json
    from scipy.stats import norm
    o = _obs(str(d)); ri = Path(d) / "run_info.json"
    if o is None or not ri.exists() or 2 * norm.sf(abs(t_full)) >= .05:
        return False
    runs = [f"loo_{x}" for x in json.load(open(ri)).get("major_loci", [])] + ["minor_all"]
    x = o[(o.readout == readout) & o.run.isin(runs)].set_index("run").t
    return len(x) == len(runs) and bool(((x.abs() >= 1.96) & (np.sign(x) == np.sign(t_full))).all())


def held9(rr):
    """Locus box for a class-term tile (Jon 2026-11-23 'Let's be consistent. use 9'): P < 0.05 and held (P < 0.05, same sign, own matched null)
    in each of the 8 major-locus removals AND with all minor loci removed (minor_all) = the same 9 runs as the mirror / dose-aligned summaries."""
    ok = rr.P_emp < .05 and rr.loo_n > 0 and rr.loo_held == rr.loo_n
    return bool(ok and rr.minor_all_P_emp < .05 and np.sign(rr.minor_all_z_emp) == np.sign(rr.z_emp))


def g3term(cnv, twas, variant, term, arm):
    g = _s3grid(cnv, twas, variant)
    if g is None:
        return None
    rr = g[(g.term == term) & (g.arm == arm)]
    if not len(rr):
        return None
    rr = rr.iloc[0]
    if term.split("~")[-1] == "n_all":
        # total genes hit does not depend on any gene annotation, so the genome-shift null cannot test it: analytic t of the joint model, two-sided
        # normal P, no locus box (Jon 2026-11-23, "Ok add total")
        from scipy.stats import norm
        t = float(rr.t_analytic); return (t, 2 * norm.sf(abs(t)), analytic_held(XS2 / "grid" / f"{cnv}_x_{twas}" / variant, f"{term}|{arm}", t))
    return (rr.z_emp, rr.P_emp, held9(rr))


XA = Path(os.environ.get("GW2_ALLCNV", ROOT.parent / "all2t_mirror_2026-11-23"))   # pooled model (all CNVs), Jon 2026-11-23 'add the all CNVs row'


@lru_cache(None)
def _asum():
    f = XA / "mirror4_cross_summary.tsv"
    return pd.read_csv(f, sep="\t") if f.exists() else None


@lru_cache(None)
def _agrid(cnv):
    f = XA / "grid" / f"{cnv}_x_{cnv}" / "all2t" / "grid_summary.tsv"
    if not f.exists():
        return None
    g = pd.read_csv(f, sep="\t"); g[["term", "arm"]] = g.readout.str.rsplit("|", n=1, expand=True); return g


def s2get(trait, size, st, stat):
    if size == "N1_inf":                  # all-CNVs row: pooled all2t run
        X = _asum()
        if X is None:
            return None
        r = X[(X.trait == f"{trait}_x_{trait}") & (X.variant == "all2t") & (X["size"] == size) & (X.stratum == st) & (X.stat == stat)]
        if not len(r):
            return None
        r = r.iloc[0]; return (r.z, r.P, r.P < .05 and r.loo_n > 0 and r.loo_held == r.loo_n)
    return g3get(trait, trait, S2VAR, size, st, stat)


def s2term(trait, term, arm):
    if term.startswith("N1_inf~"):
        g = _agrid(trait)
        if g is None:
            return None
        rr = g[(g.term == term) & (g.arm == arm)]
        if not len(rr):
            return None
        rr = rr.iloc[0]
        if term.split("~")[-1] == "n_all":
            from scipy.stats import norm
            t = float(rr.t_analytic); return (t, 2 * norm.sf(abs(t)), analytic_held(XA / "grid" / f"{trait}_x_{trait}" / "all2t", f"{term}|{arm}", t))
        return (rr.z_emp, rr.P_emp, held9(rr))
    return g3term(trait, trait, S2VAR, term, arm)


@lru_cache(None)
def _dose():
    # dose_summary.py: contrasts of the size2t joint model, plus the pooled all2t model (all-CNVs row, size N1_inf)
    fs = [f for f in (XS2 / "dose_summary.tsv", XA / "dose_summary.tsv") if f.exists()]
    return pd.concat([pd.read_csv(f, sep="\t") for f in fs], ignore_index=True) if fs else None


def dget(trait, size, st, stat):
    """Dose-aligned contrast (D4 / doseOth / doseTot) of the matched pair -> (z, P, held)."""
    X = _dose()
    if X is None:
        return None
    r = X[(X.trait == f"{trait}_x_{trait}") & (X["size"] == size) & (X.stratum == st) & (X.stat == stat)]
    if not len(r):
        return None
    r = r.iloc[0]; return (r.z, r.P, r.P < .05 and r.loo_n > 0 and r.loo_held == r.loo_n)


GREYC = "#8c8c8c"   # total-genes term of the joint model (ordinary gene), Jon 2026-11-23
TD = Path(os.environ.get("GW2_TOTDOSE", ROOT.parent / "totdose_2026-11-23"))   # totdose fit + dose_summary.py (total genes, all CNV sizes)


@lru_cache(None)
def _totdose():
    f = TD / "dose_summary.tsv"
    if not f.exists():
        return None
    X = pd.read_csv(f, sep="\t"); X = X[(X["size"] == "N1_inf") & (X.stat == "doseTot")].copy()
    X["q"] = bh(X.P.to_numpy()) if len(X) else []      # FDR family: the two traits
    return X


def total_tiles(trait):
    """Yellow total tiles of the all-CNVs row (Jon 2026-11-23 'where is the yellow text'): mirror = gene-level total mirror, shown as -z of
    r(zDEL, zDUP) so that positive = deletion and duplication effects opposite (gene-circle shift null); dose-aligned = total genes, all CNV sizes,
    totdose fit, analytic. -> ((z, P, held) mirror, (z, P, held) dose)."""
    m = d = None
    gmf = GM / trait / "gene_mirror_summary.tsv"
    if gmf.exists():
        g = pd.read_csv(gmf, sep="\t"); g = g[(g["size"] == "all sizes") & (g.stratum == "all")].set_index("stat"); t = g.loc["total"]
        m = (-t.z, t.P, t.P < .05 and t.n_loo > 0 and t.loo_held == t.n_loo)
    X = _totdose()
    if X is not None and len(X[X.trait == f"{trait}_x_{trait}"]):
        r = X[X.trait == f"{trait}_x_{trait}"].iloc[0]; d = (r.z, r.P, r.P < .05 and r.loo_n > 0 and r.loo_held == r.loo_n)
    return m, d


def total_arrows(Q, trait, dose):
    """Yellow = total (Jon 2026-11-23). Mirror diagonal: gene-level total mirror r(zDEL, zDUP), all genes and CNV sizes, z vs the gene-circle
    shift null (as the main figure). Dose-aligned diagonal (alternate figure only): total genes hit, all CNV sizes, (t_DEL + t_DUP)/sqrt2 of the
    totdose fit, analytic; thick = FDR < 5% across the two traits, thin = P < 0.05 and held in every leave-one-locus-out run."""
    gmf = GM / trait / "gene_mirror_summary.tsv"
    if gmf.exists():
        g = pd.read_csv(gmf, sep="\t"); g = g[(g["size"] == "all sizes") & (g.stratum == "all")].set_index("stat"); tot = g.loc["total"]
        if tot.P < .05:
            dd = np.array([-1, 1]) / np.sqrt(2) if tot.obs < 0 else np.array([1, 1]) / np.sqrt(2)
            Q.append((dd * abs(tot.z), "#f2c230", True, "full")); DRAWN.add((("tot", "mirror"), trait, "N1_inf"))
    X = _totdose() if dose else None
    if X is not None:
        r = X[X.trait == f"{trait}_x_{trait}"]
        if len(r):
            r = r.iloc[0]; held = r.P < .05 and r.loo_n > 0 and r.loo_held == r.loo_n
            if r.q < .05 or held:
                Q.append((np.array([1, 1]) / np.sqrt(2) * r.z, "#f2c230", bool(r.q < .05), "full")); DRAWN.add((("tot", "dose"), trait, "N1_inf"))


def schematic_dose(ax, trait):
    """Alternate Fig. 3c (Jon 2026-11-23: 'make an alternate version of the figure with that joint model ... yellow (CNV) and orange/brown (twas)
    and green (weak/no twas) arrows ... that would follow a dose-aligned axis'). All arrows on the dose-aligned axis DUP = DEL (up-right = both
    dosages raise the trait), from the size2t joint model: yellow = total genes (DEL + DUP)/sqrt2, analytic; orange / brown = D4 of the 2x2 TWAS
    tiles, LC / HC; green = 'weak or no TWAS' HC (DEL + DUP)/sqrt2. Arrow rule as Fig. 3c (pick_arrow); shorter arrows on top.
    null model: D4, green: genome shift of the gene annotation (dose_summary.py); yellow: analytic (no annotation, no shift null)."""
    L = 7.5; ax.plot([-L, L], [-L, L], color=GREY, lw=.7, ls="--", zorder=1); ax.axhline(0, color=INK, lw=.9); ax.axvline(0, color=INK, lw=.9)
    u = np.array([1, 1]) / np.sqrt(2); Q = []
    d_ = lambda stat: (lambda t, s, st_: dget(t, s, st_, stat))
    a_ = lambda stat: (lambda t, s, st_: s2get(t, s, st_, stat))
    o_ = lambda term, arm: (lambda t, s, st_: s2term(t, f"{s}~{term}", arm))
    # one arrow per colour (component): dose-aligned diagonal, else the DEL / DUP axis arrows of the same component
    m = np.array([-1, 1]) / np.sqrt(2)          # mirror diagonal (up-left = duplication raises, deletion lowers)
    ax.plot([-L, L], [L, -L], color=GREY, lw=.7, ls="--", zorder=1)
    # each colour: dose-aligned and mirror diagonals (every one that passes), else the DEL / DUP axes. Total genes: mirror readout of n_all,
    # analytic like its dose contrast (g3term); TWAS: M4; weak or no TWAS HC: the term's '|mirror' readout
    # Jon 2026-11-23 ('yes'): yellow = the TOTAL on both diagonals (mirror: gene-level total mirror, as the main figure; dose-aligned: total
    # genes, all CNV sizes, from the totdose fit, analytic); the total-genes term of the joint model (ordinary gene) is a component in grey
    total_arrows(Q, trait, dose=True)
    # grey total-genes arrows removed (Jon 2026-11-23: "Just remove the gray arrows"); the term stays in the model and in panel a
    for st, col in (("LC", "#f28e2b"), ("HC", "#7b3294")):
        component(Q, trait, (st, "HC" if st == "LC" else "LC"), col, "full", [(d_("D4"), u, ("D4", st)), (a_("M4"), m, ("M4", st))],
                  [])   # Jon 2026-11-23: 'Purple and Orange are meant to indicate an axis where there is CNV-TWAS convergence' -> diagonals only, no DEL/DUP arm fallback
    component(Q, trait, ("HC",), "#1b7837", "arm", [(d_("doseOth"), u, ("doseOth", "HC")), (o_("cls:other:HC", "mirror"), m, ("mirOth", "HC"))],
              [(o_("cls:other:HC", "DEL"), (1, 0), ("othDEL", "HC")), (o_("cls:other:HC", "DUP"), (0, 1), ("othDUP", "HC"))])
    flush_arrows(ax, Q, L)
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    tl = TL[trait] if trait == "bmi" else "height"
    # Jon 2026-11-23: 'not clear that vertical axis is DUP and horizontal axis is DEL ... you need to read 4 text labels': each axis gets one
    # title in the panel a colours (DEL red, horizontal; DUP blue, vertical) at its positive end, + / - at the ends (+ = raises the trait)
    ax.text(L, .35, f"DEL effect on {tl}", ha="right", va="bottom", size=6.4, weight="bold", color=RED)
    ax.text(.35, L, f"DUP effect on {tl}", ha="left", va="top", size=6.4, weight="bold", color=BLUE)
    for x_, y_, t_, ha_, va_ in ((L, -.35, "+", "right", "top"), (-L, -.35, "−", "left", "top"), (-.35, L, "+", "right", "top"),
                                 (-.35, -L, "−", "right", "bottom")):
        ax.text(x_, y_, t_, ha=ha_, va=va_, size=11, weight="bold", color=INK)   # Jon: larger + / -
    ax.set_title(TL[trait], loc="left", weight="bold", size=9)


def flush_arrows(ax, Q, L):
    """Draw queued arrows longest first (shorter on top). Jon 2026-11-23: 'When 2 overlapping arrows are the same length, make one longer just so we
    can see it': an arrow on the same axis and direction whose drawn length is within one head length of a shorter one on top of it is lengthened
    to that length + 1.3 head lengths, so both heads show (a display offset only; lengths are |z| otherwise)."""
    dl = lambda v: min(max(np.linalg.norm(v), .25 * L), L - .3)      # drawn length (shanta_arrow caps)
    A = [(np.asarray(v, float), c, s, k) for v, c, s, k in Q if np.isfinite(np.linalg.norm(v)) and np.linalg.norm(v) > 0]
    A = [(v / np.linalg.norm(v) * dl(v), *r) for v, *r in A]          # vectors at their drawn length
    A.sort(key=lambda a: -np.linalg.norm(a[0]))
    head, cap = .10 * L, L - .3
    for i in range(len(A) - 2, -1, -1):          # each longer arrow must clear every shorter one on the same axis and direction
        u = A[i][0] / np.linalg.norm(A[i][0])
        for j in range(i + 1, len(A)):
            w = A[j][0]; nw = np.linalg.norm(w)
            if np.allclose(u, w / nw, atol=1e-6) and np.linalg.norm(A[i][0]) - nw < head:
                if nw + 1.3 * head <= cap:
                    A[i] = (u * (nw + 1.3 * head), *A[i][1:])          # lengthen the lower arrow
                else:
                    A[j] = (w / nw * (np.linalg.norm(A[i][0]) - 1.3 * head), *A[j][1:])   # at the cap: shorten the upper one instead
    for k, (v, col, sig, kind) in enumerate(A):
        shanta_arrow(ax, v, col, sig, kind, L, z=2 + k)


def component(Q, trait, strata, col, kind, diag, arms):
    """Jon 2026-11-23 ('If weak/no twas is significant for Del or dup and mirror/dose aligned, should we draw multiple arrows or just the diagonal
    arrow?' -> 'Yes' to: one arrow per colour): draw the diagonal (mirror or dose-aligned) if it passes pick_arrow; only if it does not, draw any
    DEL / DUP axis arrow that passes. diag = (fetch, unit vector); arms = [(fetch, unit vector), ...]."""
    # diag may be a list of diagonals (alternate figure: dose-aligned and mirror, Jon 2026-11-23 'mirror diagonal'): every one that passes is
    # drawn; the DEL / DUP axis arrows only when none passes
    hit = False
    for d_item in (diag if isinstance(diag, list) else [diag]) if diag else []:
        fetch, u = d_item[:2]; key = d_item[2] if len(d_item) > 2 else None
        pk = pick_arrow(fetch, trait, strata, key)
        if pk is not None:
            Q.append((np.asarray(u, float) * pk[0], col, pk[1], kind)); hit = True
    if hit:
        return
    for a_item in arms:
        fetch, u = a_item[:2]; key = a_item[2] if len(a_item) > 2 else None   # key: the panel a tiles behind an axis arrow
        pk = pick_arrow(fetch, trait, strata, key)
        if pk is not None:
            Q.append((np.asarray(u, float) * pk[0], col, pk[1], kind))


DRAWN = set()   # (key, trait, size) of every drawn arrow; panel a colours only these tiles (Jon 2026-11-23 '(only for the arrows that are drawn)')


def pick_arrow(fetch, trait, strata, key=None):
    """Arrow rule for Fig. 3c (Jon 2026-11-23): '1. When single gene and >2 gene trend in the same direction but only one is significant, draw
    the thick arrow in either case. 2. Only draw thin arrows when there is nominal P<0.05 that is locus-robust.' fetch(trait, size, stratum) ->
    (z, P, held). Thick = FDR < 5% in either CNV size (BH within the statistic across trait x CNV size x strata, the same family as the Fig. 3a
    tiles); thin = P < 0.05 and held in every leave-one-locus-out run, not FDR; nothing else is drawn. Returns (z, thick) or None; z = the
    qualifying size with the larger |z|. strata = the strata of the family (e.g. ('LC', 'HC')); the arrow's own stratum is strata[0]."""
    fam = [(t, s, st) for t in TRAITS for s, _ in S2ROWS for st in strata]
    vals = {k: fetch(*k) for k in fam}
    ok = [k for k in fam if vals[k] is not None and np.isfinite(vals[k][1])]
    q = dict(zip(ok, bh([vals[k][1] for k in ok]))) if ok else {}
    mine = [(vals[k], q[k], k) for k in ok if k[0] == trait and k[2] == strata[0]]
    for sel, thick_ in (([m for m in mine if m[1] < .05], True), ([m for m in mine if m[0][1] < .05 and m[0][2]], False)):
        if sel:
            v, _, k = max(sel, key=lambda m: abs(m[0][0]))
            if key is not None:
                DRAWN.add((key, trait, k[1]))
            return v[0], thick_
    return None



def bh(P):
    """Benjamini-Hochberg q values."""
    P = np.asarray(P, float); n = len(P); o = np.argsort(P)
    q = np.minimum.accumulate((P[o] * n / np.arange(1, n + 1))[::-1])[::-1]; out = np.empty(n); out[o] = np.minimum(q, 1); return out


class Tiles:
    """Collects tiles, then draws them with '**' where BH q < 0.05 within the tile's family (Shanta et al., fig2_shanta.R: star(P, q) =
    '**' if q < 0.05, '*' if P < 0.05; BH within each statistic across its rows). Jon 2026-11-22: '* p<0.05 under null and ** fdr < 5%'."""
    def __init__(self):
        self.items = []

    def add(self, ax, x, y, w, h, c, fam, fs=5.6, sigcol=None):
        self.items.append((ax, x, y, w, h, c, fam, fs, sigcol))

    def draw(self):
        fams = {}
        for i, it in enumerate(self.items):
            if it[5] is not None and np.isfinite(it[5][1]):
                fams.setdefault(it[6], []).append(i)
        q = {}
        for ix in fams.values():
            q.update(zip(ix, bh([self.items[i][5][1] for i in ix])))
        for i, (ax, x, y, w, h, c, fam, fs, sigcol) in enumerate(self.items):
            tile3(ax, x, y, w, h, c, fs, q=q.get(i, np.nan), sigcol=sigcol)


def tile3(ax, x, y, w, h, c, fs=5.6, q=np.nan, sigcol=None):
    """Tile coloured by z (capped at +-4); text = z, '*' P < 0.05, '**' BH q < 0.05 in its family; box = P < 0.05 and held when each major
    recurrent locus is removed."""
    if c is None or not np.isfinite(c[0]):
        ax.add_patch(Rectangle((x, y), w, h, fc="#eeeeee", ec="white", lw=1)); return
    ax.add_patch(Rectangle((x, y), w, h, fc=CMAP3(NORM3(np.clip(c[0], -4, 4))), ec="white", lw=1))
    st = "**" if q < .05 else ("*" if c[1] < .05 else "")
    t = ax.text(x + w / 2, y + h / 2, f"{c[0]:.1f}" + st, ha="center", va="center", size=fs,
                color="white" if abs(c[0]) > 2.6 else INK, weight="bold" if c[1] < .05 else "normal")
    if sigcol and (sigcol[1] in DRAWN if isinstance(sigcol, tuple) else q < .05):
        sigcol = sigcol[0] if isinstance(sigcol, tuple) else sigcol
        # Jon 2026-11-23: 'make the z value in the tile match the color [of its arrow] when it's significant' (FDR < 5%); white outline
        # keeps the coloured text legible on the tile colour
        import matplotlib.patheffects as pe
        # Jon 2026-11-23: 'Instead of making the text yellow orange green ... place a small colored box around the text to signify the arrow it
        # corresponds to': text keeps its normal colour; a box in the arrow's colour marks the tile whose arrow is drawn in panel b
        t.set_bbox(dict(boxstyle="round,pad=0.18", fc="none", ec=sigcol, lw=1.6))
    if c[2]:
        ax.add_patch(Rectangle((x + .04, y + .04), w - .08, h - .08, fill=False, ec=INK, lw=1.3))


def fig3(dose=False):
    """Jon 2026-11-22 mockup: 'For now let's only show matched TWAS. And instead of "all cnvs" split them by 1 gene and >1 gene.'
    a: per trait (height CNV effects x height TWAS genes; BMI x BMI), LC and HC genes, rows = >= 2-gene and single-gene CNVs. Each unit: 2x2 tiles
    (rows TWAS z > +3 'up' and z < -3 'down', columns DUP and DEL; z of the signed-count term, joint size2 model), up-down strip ((up - down)/sqrt2
    per dosage) and the merged four-tile mirror M4 = (-DEL up + DUP up + DEL down - DUP down)/2. b: constraint, other constrained genes
    (|TWAS z| < 3) by trait and CNV size. c: schematic (total mirror, CNV-TWAS convergence, constraint), no caption (Jon: 'remove the caption').
    null model: genome shift of the gene annotation, 2,000 shared shifts; every statistic recomputed in each replicate (mirror4_summary.py)."""
    # Jon 2026-11-23: "I think we need everything from the joint model in the Main figure. drop the up-down row, we don't use it, and replace it
    # with 'No-TWAS' row. Insert the -3<z<+3 row between z>3 and z<-3." Rows: z > +3, -3 < z < +3, z < -3, no TWAS z = all 8 terms per CNV size
    # of the size3c joint model; the mirror tile (a contrast of the up / down terms of the same fit) sits beside the three TWAS-z rows.
    FH = 8.6 if S2VAR == "size2t" else 8.4; fy = 7.4 / FH    # figure taller (all-CNVs row); lower panels keep their size in inches
    fig = plt.figure(figsize=(7.2, FH))
    ax = fig.add_axes([.01, 4.45 / FH, .98, (FH - 4.62) / FH]); ax.set_aspect("equal"); ax.axis("off")   # Jon: lower the title, fill the whitespace
    # Jon 2026-11-22 second mockup: one frame per trait x CNV-size row (frames share borders), size label and TWAS band labels inside the
    # height frame, 'TWAS' outside on the left, 'CNV' over the LC / HC gene headers.
    UW, GAP, SH, CG, LA = 4.35, .35, 0, .9, 2.7            # unit width, tile-mirror gap, (no strip), LC-HC gap, label column (height frame)
    ROWS4 = [("up", "z>+3"), (WEAK, "−3<z<+3"), ("down", "z<−3"), ("none", "no TWAS z")]
    if S2VAR == "size2t":
        # Jon 2026-11-23: "The no TWAS and the weak TWAS always trend in the same direction ... collapse them ... put them in the model row where
        # weak TWAS currently sits. Also ... total gene count term ... in bottom row where no-TWAS currently sits." Model size2t, per CNV size:
        # total genes hit + up / down x LC / HC + other (weak or no TWAS z) HC; other LC genes are the reference absorbed by the total term.
        ROWS4 = [("up", "z>+3"), ("down", "z<−3"), ("other", "weak or no TWAS"), ("n_all", "total genes")]   # Jon 2026-11-23: weak or no TWAS to row 3
        # Jon 2026-11-23: "the figure needs to have a weak or no TWAS tile for mirror and for dose-aligned. So there would be mirror CNV-TWAS,
        # mirror CNV, dose aligned CNV-TWAS and dose-aligned CNV, to the right of each 2x8 grid." Right block: columns mirror | dose-aligned,
        # rows CNV-TWAS (M4 / D4 of the up and down terms) | CNV (weak or no TWAS HC term; in LC, whose weak-or-no-TWAS genes are the reference,
        # the total-genes term that carries them, analytic).
        SW = 1.45; RL = .45; UW = 2 + GAP + 2 * SW + RL   # RL: column for the block's row labels (Jon: 'add row labels', coloured as the arrows)
    UH = 4.0; BW = 2 * UW + CG; PT, PB, PS = .55, .3, .35   # frame padding top / bottom / sides
    STEP = UH + PT + PB                                    # rows abut: next frame top = this frame bottom
    TA = Tiles()   # FDR family = one statistic (tile position, arm or mirror) across the 8 trait x CNV-size x LC/HC units
    # Jon 2026-11-23: 'add a full grid row for all CNVs above ... order it all genes, single gene CNVs, >2 gene CNVs'. The all-CNVs row is the
    # pooled model (all2t); its tiles form their own FDR families (tag 'all') so the size-row stars do not change; arrows still come from the size rows
    AROWS = ([("N1_inf", "All CNVs")] if S2VAR == "size2t" else []) + [("N1_1", "Single-gene CNVs"), ("N2_inf", "≥2-gene CNVs")]
    tag = lambda size, fam: fam + (("all",) if size == "N1_inf" else ())
    fx = {"height": (-LA - PS, BW + PS)}; fx["bmi"] = (fx["height"][1], fx["height"][1] + BW + 2 * PS)
    for b, trait in enumerate(TRAITS):
        x_l, x_r = fx[trait]; bx = x_r - PS - BW
        # Jon 2026-11-23: 'Remove CNV from the subtitle'
        ax.text((bx + bx + BW) / 2, UH + PT + .75, TL[trait], ha="center", va="bottom", size=10, weight="bold")
        for i, (size, slab) in enumerate(AROWS):
            y0 = -i * STEP
            ax.add_patch(Rectangle((x_l, y0 - PB), x_r - x_l, STEP, fill=False, ec="#c8c8c8", lw=1.1))
        for j, st in enumerate(("LC", "HC")):
            x0 = bx + j * (UW + CG)
            ax.text(x0 + UW / 2, UH + PT + .1, f"{st} genes", ha="center", va="bottom", size=7.4, weight="bold")
            ax.text(x0 + .5, UH + .06, "DUP", ha="center", va="bottom", size=6, color=BLUE)
            ax.text(x0 + 1.5, UH + .06, "DEL", ha="center", va="bottom", size=6, color=RED)
            if S2VAR == "size2t":
                ax.text(x0 + 2 + GAP + SW / 2, UH + .06, "mirror", ha="center", va="bottom", size=5.4)
                ax.text(x0 + 2 + GAP + 1.5 * SW, UH + .06, "dose-aligned", ha="center", va="bottom", size=5)
            else:
                ax.text(x0 + 3 + GAP, UH + .06, "mirror", ha="center", va="bottom", size=6)
            for i, (size, _) in enumerate(AROWS):
                y0 = -i * STEP
                for r_, (cls, _) in enumerate(ROWS4):
                    for k, dos in enumerate(("DUP", "DEL")):
                        yy = y0 + UH - 1 - r_
                        if cls == "n_all":            # one total-count term per CNV size and dosage, shared by LC and HC: drawn once, under LC
                            if st == "LC":
                                TA.add(ax, x0 + k, yy, 1, 1, s2term(trait, f"{size}~n_all", dos), tag(size, ("a", "n_all", dos)))
                            continue
                        if cls == "other" and st == "LC":   # reference class (absorbed into the total term)
                            ax.add_patch(Rectangle((x0 + k, yy), 1, 1, fc="#f2f2f2", ec="white", lw=1))
                            ax.text(x0 + k + .5, yy + .5, "ref.", ha="center", va="center", size=5.2, color=GREY)
                            continue
                        # Jon 2026-11-23 ('If there is a purple arrow, there should be a purple box on the grid', via ConGA): an arrow that
                        # falls back to the DEL / DUP axis outlines the class tiles it comes from (arm = up - down of that dosage; or the
                        # weak-or-no-TWAS HC term itself)
                        sc = None
                        if cls == "other" and st == "HC":   # green only: orange / purple never fall back to an axis (Jon 2026-11-23)
                            sc = ("#1b7837", (("oth" + dos, "HC"), trait, size))
                        TA.add(ax, x0 + k, yy, 1, 1, s2term(trait, f"{size}~cls:{cls}:{st}", dos), tag(size, ("a", cls, dos)), sigcol=sc)
                if S2VAR == "size2t":
                    xm = x0 + 2 + GAP
                    if st == "LC" and size == "N1_inf":   # all-CNVs row: the yellow totals (Jon 2026-11-23)
                        (mC, dC), fC, lC = total_tiles(trait), "tot", "CNV total"
                    elif st == "LC":
                        mC, dC, fC, lC = s2term(trait, f"{size}~n_all", "mirror"), dget(trait, size, "total", "doseTot"), "total", "CNV (total)"
                    else:
                        mC, dC, fC, lC = s2term(trait, f"{size}~cls:other:HC", "mirror"), dget(trait, size, "HC", "doseOth"), "otherHC", "CNV"
                    for cx, top, bot, fam in ((xm, s2get(trait, size, st, "M4"), mC, "mirror"), (xm + SW, dget(trait, size, st, "D4"), dC, "dose")):
                        kt = ("M4" if fam == "mirror" else "D4", st); kb = ("mirOth" if fam == "mirror" else "doseOth", "HC")
                        TA.add(ax, cx, y0 + 2.04, SW - .04, 1.96, top, tag(size, ("a", fam, "TWAS")), fs=7,
                               sigcol=("#f28e2b" if st == "LC" else "#7b3294", (kt, trait, size)))
                        sb = (("#f2c230", (("tot", fam), trait, size)) if fC == "tot" else None if st == "LC" else ("#1b7837", (kb, trait, size)))
                        TA.add(ax, cx, y0, SW - .04, 1.96, bot, tag(size, ("a", fam, fC)), fs=7, sigcol=sb)
                    # row labels of the block, coloured as the panel b arrows (Jon 2026-11-23): CNV-TWAS = orange (LC) / brown (HC);
                    # CNV only = grey (LC: total-genes term, ordinary gene) / green (HC: weak or no TWAS)
                    xr = xm + 2 * SW + .08
                    ax.text(xr, y0 + 3.02, "CNV–TWAS", rotation=90, ha="left", va="center", size=5.2, weight="bold",
                            color=INK)   # Jon 2026-11-23: 'remove all green, yellow orange text fonts'
                    lab_, col_ = ("CNV total" if (st == "LC" and size == "N1_inf") else "CNV only"), INK
                    ax.text(xr, y0 + .98, lab_, rotation=90, ha="left", va="center", size=5.2, weight="bold", color=col_)
                elif dose:   # alternate figure: mirror (top) and dose-aligned D4 (bottom) of the same 2x2, stacked beside the four rows
                    TA.add(ax, x0 + 2 + GAP, y0 + 2.05, 2, 1.95, s2get(trait, size, st, "M4"), ("a", "M4"), fs=8)
                    TA.add(ax, x0 + 2 + GAP, y0, 2, 1.95, dget(trait, size, st, "D4"), ("a", "D4"), fs=8)
                    ax.text(x0 + 3 + GAP, y0 + .12, "dose-aligned", ha="center", va="bottom", size=4.6, color=GREY)
                else:
                    TA.add(ax, x0 + 2 + GAP, y0 + 1.5, 2, 2, s2get(trait, size, st, "M4"), ("a", "M4"), fs=8)   # beside the three TWAS-z rows
    for i, (_, slab) in enumerate(AROWS):
        y0 = -i * STEP
        ax.text(-LA - PS + .2, y0 + UH + PT - .12, slab, ha="left", va="top", size=6.8, weight="bold")
        # Jon 2026-11-23: 'the 2x2 at the top are the only rows that are based solely on TWAS ... move the label here. Use the same font size as
        # CNV-TWAS': rotated label beside the z > +3 / z < -3 rows
        ax.text(-1.35, y0 + UH - 1, "TWAS", ha="center", va="center", size=5.2, weight="bold", rotation=90)
        for r_, (cls_, lab) in enumerate(ROWS4):
            ax.text(-.15, y0 + UH - .5 - r_, lab, ha="right", va="center", size=5.6,
                    color=INK)   # Jon 2026-11-23: no coloured label text
    ybot = -(len(AROWS) - 1) * STEP - PB
    ax.set_xlim(-LA - PS - .9, fx["bmi"][1] + .1); ax.set_ylim(ybot - .2, UH + PT + 1.75)
    cax = fig.add_axes([.08, .535 * fy, .2, .009 * fy]); sm = plt.cm.ScalarMappable(cmap=CMAP3, norm=NORM3)
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal"); cb.set_label("z against genome-shift null", size=6); cb.ax.tick_params(labelsize=5.6)
    foot = ("tile = z of each gene-count term of one joint model per trait (8 per CNV size: TWAS z > +3, −3 < z < +3, z < −3 or no TWAS z,\n"
            "x LC / HC = less / highly constrained genes);  mirror = (−DEL up + DUP up + DEL down − DUP down)/2 of the z > +3 (up) and z < −3 (down) terms;\n")
    if S2VAR == "size2t":
        foot = ("tile = z of each term of one joint model per trait and CNV size (all CNVs: the same model pooled over size): total genes hit,\n"
                "and genes with TWAS z > +3 (up), "
                "z < −3 (down) or weak or no TWAS (−3 ≤ z ≤ +3 or no z), by LC / HC = less / highly constrained genes;\n"
                "ref. = weak-or-no-TWAS LC genes, carried by the total;  mirror = (−DEL up + DUP up + DEL down − DUP down)/2;\n"
                "dose-aligned = (DEL up + DUP up + DEL down + DUP down)/2;  CNV tiles: mirror (DUP − DEL)/√2 and dose-aligned (DUP + DEL)/√2\n"
                "of the weak-or-no-TWAS HC term (LC: of the total-genes term); CNV total: gene-level DEL–DUP anticorrelation (mirror)\n"
                "and a model of total genes hit only (dose-aligned); coloured box: tile shown as an arrow in b;\n" +
                "total genes: analytic t (no gene annotation, so no shift null); other tiles: z against the null;\n")
    fig.text(.33, .538 * fy, foot + "* P < 0.05;  ** FDR < 5% (BH within each statistic);  box: P < 0.05 and held when each of 8 major loci, or all minor loci, is removed",
             size=5.3, va="center")
    SHOW_B = False   # Jon 2026-11-23: 'Drop panel b' (its weak-TWAS genes are now the 'weak or no TWAS' row of panel a); kept for reference
    if SHOW_B:
        # b: constraint, other HC genes (|TWAS z| < 3), own TWAS, by CNV size
        # b: other genes (|TWAS z| < 3), LC and HC (Jon 2026-11-22), rows = trait x CNV size, columns = LC (DUP, DEL) | HC (DUP, DEL)
        axB = fig.add_axes([.02, .07 * fy, .40, .36 * fy]); axB.set_aspect("equal"); axB.axis("off")
        brows = [(t, size, slab) for t in TRAITS for size, slab in S2ROWS]; TB = Tiles()   # FDR family = each stratum x dosage across the 4 rows
        for i, (t, size, slab) in enumerate(brows):
            y = len(brows) - 1 - i
            for j, st in enumerate(("LC", "HC")):
                for k, dos in enumerate(("DUP", "DEL")):
                    TB.add(axB, j * 2.3 + k, y, 1, 1, s2term(t, f"{size}~cls:{WEAK}:{st}", dos), ("b", st, dos), fs=6)
            axB.text(-.1, y + .5, slab.replace(" CNVs", ""), ha="right", va="center", size=6)
            if size == S2ROWS[0][0]:
                axB.text(-2.25, y, TL[t], ha="center", va="center", size=6.8, weight="bold", rotation=90)
        TB.draw()
        axB.plot([-2.45, 4.3], [2, 2], color=GREY, lw=.6)
        for j, st in enumerate(("LC", "HC")):
            axB.text(j * 2.3 + 1, 4.55, f"{st} genes", ha="center", va="bottom", size=7, weight="bold")
            axB.text(j * 2.3 + .5, 4.08, "DUP", ha="center", va="bottom", size=6.2, color=BLUE)
            axB.text(j * 2.3 + 1.5, 4.08, "DEL", ha="center", va="bottom", size=6.2, color=RED)
        axB.set_xlim(-2.6, 4.4); axB.set_ylim(-.1, 6.3)
        axB.text(.9, 5.15, "Genes with weak TWAS signal\n−3 < z < +3", ha="center", va="bottom", size=7.4, weight="bold")   # Jon 2026-11-22
    axC = [fig.add_axes([.07 + .47 * j, .085 * fy, .4, .37 * fy]) for j in range(2)]   # panel b dropped: arrows panel takes the full width
    DRAWN.clear()
    for a_, t in zip(axC, TRAITS):
        schematic_dose(a_, t) if dose else schematic(a_, t, caption=False)
    TA.draw()   # after panel b, so tile values take their arrow colour only when that arrow is drawn
    hs = [Line2D([], [], color="#f2c230", lw=3), Line2D([], [], color="#f28e2b", lw=3), Line2D([], [], color="#7b3294", lw=3), Line2D([], [], color="#1b7837", lw=3)]
    lb = ["total mirror (DEL vs DUP effects across genes)", "CNV–TWAS convergence, LC genes", "CNV–TWAS convergence, HC genes",
          "constraint   (thick, filled head: FDR < 5%; thin, open head: P < 0.05, locus-robust)"]
    if dose or S2VAR == "size2t":
        lb[0] = "total (all genes: mirror = gene-level DEL vs DUP; dose-aligned = total genes, analytic t)"
    if dose:
        lb = ["total (all genes: mirror = gene-level DEL vs DUP; dose-aligned = total genes, analytic t)", "CNV–TWAS, LC genes (TWAS z > +3 and z < −3)",
              "CNV–TWAS, HC genes (TWAS z > +3 and z < −3)", "weak or no TWAS, HC genes   (thick: FDR < 5%; thin: P < 0.05, locus-robust)"]
    if S2VAR == "size2t":   # grey component: total-genes term (ordinary gene) of the joint model
        pass   # no grey arrows, so no legend entry (Jon 2026-11-23)
    fig.legend(hs, lb, loc="lower center", ncol=2, frameon=False, fontsize=5.8, bbox_to_anchor=(.5, -.003))
    fig.text(.01, .975, "a", size=11, weight="bold"); fig.text(.01, .465 * fy, "b", size=11, weight="bold")
    lab3 = os.environ.get("GW2_FIG3_LABEL")     # cohort label: meta runs only, from meta_info.json; GW2_FIG3_LABEL overrides ("" = none)
    if META3 and lab3 is None:
        import json
        mi = Path(META3) / "meta_info.json"; lab3 = json.load(open(mi)).get("label", "") if mi.exists() else ""
    if lab3:
        fig.text(.99, .992, lab3, ha="right", va="top", size=6.4, color=GREY)
    save(fig, ("Fig3_joint_model_dose" if dose else "Fig3_joint_model") + ("_meta" if META3 else ""))


def figS_full_grid():
    """Jon 2026-11-22: 'This supplemental should do a full grid of LC genes, HC genes, all genes x z > +3, z < -3, -3 < z < 3'. One unit per
    CNV trait x CNV size (rows) and TWAS trait (columns); within a unit, three blocks (LC, HC, all genes), each with tiles for TWAS z > +3,
    z < -3 and -3 < z < +3 (rows) x DUP, DEL (columns), the up-down strip ((up - down)/sqrt2) and the four-tile mirror. LC / HC from model
    size3c, all genes from model size3a (genes without a TWAS z are a fourth class). * P < 0.05; ** BH q < 0.05 within
    each statistic across the 16 units; box: P < 0.05 and held when each major recurrent locus is removed.
    null model: genome shift of the gene annotation, 2,000 shared shifts; every statistic recomputed in each replicate (mirror4_summary.py)."""
    fig = plt.figure(figsize=(7.2, 6.7))
    ax = fig.add_axes([.01, .13, .98, .82]); ax.set_aspect("equal"); ax.axis("off")
    BL3 = [("up", "z>+3"), ("down", "z<−3"), (WEAK, "−3<z<+3"), ("none", "no TWAS z")]   # Jon: 'maybe add a category of other genes (no TWAS z)'
    SH, G1, G2, MW = .5, .25, .7, 1.6         # strip height, tiles-mirror gap, block gap, mirror tile width (= height of the two TWAS rows / 1.25)
    BW = 2 + G1 + MW; UW = 3 * BW + 2 * G2; UG = 1.6; RH = 4 + .12 + SH; RG = 1.5
    blocks = [("LC", "size3c", "LC genes"), ("HC", "size3c", "HC genes"), ("all", "size3a", "all genes")]
    rows = [(c, size, slab) for c in TRAITS for size, slab in S2ROWS]
    T = Tiles()
    for i, (cnv, size, slab) in enumerate(rows):
        y0 = -i * (RH + RG)                    # bottom of the three TWAS rows
        ax.text(-2.3, y0 + 2, f"{TL[cnv]} CNVs\n{slab.replace(' CNVs', '')}", ha="right", va="center", size=6.4, weight="bold")
        for j, twas in enumerate(TRAITS):
            ux = j * (UW + UG)
            if i == 0:
                ax.text(ux + UW / 2, 5.55, f"{TL[twas]} TWAS genes", ha="center", va="bottom", size=7.4, weight="bold")
            for b_, (st, var, blab) in enumerate(blocks):
                x0 = ux + b_ * (BW + G2)
                if i == 0:
                    ax.text(x0 + BW / 2, 4.75, blab, ha="center", va="bottom", size=6.2, weight="bold")
                    ax.text(x0 + .5, 4.08, "DUP", ha="center", va="bottom", size=5, color=BLUE)
                    ax.text(x0 + 1.5, 4.08, "DEL", ha="center", va="bottom", size=5, color=RED)
                    ax.text(x0 + 2 + G1 + MW / 2, 4.08, "mirror", ha="center", va="bottom", size=5)
                for r_, (cls, _) in enumerate(BL3):
                    for k, dos in enumerate(("DUP", "DEL")):
                        T.add(ax, x0 + k, y0 + 3 - r_, 1, 1, g3term(cnv, twas, var, f"{size}~cls:{cls}:{st}", dos), (st, cls, dos), fs=4.6)
                for k, arm in enumerate(("armDUP", "armDEL")):
                    T.add(ax, x0 + k, y0 - .12 - SH, 1, SH, g3get(cnv, twas, var, size, st, arm), (st, arm), fs=4.4)
                T.add(ax, x0 + 2 + G1, y0 + 2, MW, 2, g3get(cnv, twas, var, size, st, "M4"), (st, "M4"), fs=6.4)
            if j == 0:
                for r_, (_, lab) in enumerate(BL3):
                    ax.text(ux - .08, y0 + 3.5 - r_, lab, ha="right", va="center", size=4.8)
                ax.text(ux - .08, y0 - .12 - SH / 2, "up−down", ha="right", va="center", size=4.6)
        if i == 1:
            ax.plot([-6.6, 2 * UW + UG], [y0 - .12 - SH - RG / 2] * 2, color=GREY, lw=.6)
    T.draw()
    ax.axvline(UW + UG / 2, color="#c8c8c8", lw=.8)
    ax.set_xlim(-6.7, 2 * UW + UG + .1); ax.set_ylim(-(len(rows) - 1) * (RH + RG) - .12 - SH - .3, 6.3)
    cax = fig.add_axes([.08, .075, .2, .01]); sm = plt.cm.ScalarMappable(cmap=CMAP3, norm=NORM3)
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal"); cb.set_label("z against genome-shift null", size=6); cb.ax.tick_params(labelsize=5.6)
    fig.text(.33, .078, "tile = z of the signed-count term for genes with TWAS z > +3, z < −3, −3 < z < +3, or without a TWAS z (no TWAS z),\n"
             "joint model by CNV size (LC / HC: one model; all genes: a second model);  up−down = (up − down)/√2 per dosage;\n"
             "mirror = (−DEL up + DUP up + DEL down − DUP down)/2;  * P < 0.05;  ** FDR < 5% (BH within each statistic);\n"
             "box: P < 0.05 and held when each major recurrent locus is removed", size=5.2, va="center")
    fig.text(.5, .965, "CNV–TWAS tiles by gene constraint, every CNV trait by every TWAS trait", ha="center", size=8.5, weight="bold")
    save(fig, "FigS_full_grid")


# ------------------------------------------------------------------ TWAS z-band convergence (Shanta et al. Fig. 3D style, Jon 2026-11-22)
ZB = Path(os.environ.get("GW2_ZBAND", ROOT.parent / "zband_by_constraint_2026-11-22"))
BL = ["≤ −3", "−3 to −2", "−2 to −1", "−1 to 1", "1 to 2", "2 to 3", "≥ 3"]
FAMS = [("TWAS z band, all genes", "all genes"), ("TWAS z band, less constrained", "less-constrained genes"), ("TWAS z band, constrained", "constrained genes")]


def zband_panel(axD, axU, trait, fam, lab):
    """CNV burden z of the genes in each TWAS z band (each band alone, adjusted for total gene count; z against the genome-shift null),
    deletions (red) and duplications (blue), shading by band order; title = pre-specified trend: mirror = (trend DUP - trend DEL)/sqrt2,
    P and runs held when each major recurrent locus is removed (as fig2_scz_zband.R)."""
    f = ZB / f"band_trend_{trait}.tsv"
    if not f.exists():
        for ax in (axD, axU):
            ax.text(.5, .5, "running", transform=ax.transAxes, ha="center", color=GREY)
        return
    T = pd.read_csv(f, sep="\t"); T = T[T.family == fam]
    if not len(T):
        for ax in (axD, axU):
            ax.text(.5, .5, "running", transform=ax.transAxes, ha="center", color=GREY); ax.set_axis_off()
        return
    ymax = max(2.2, np.nanmax(np.abs(T[T.stat.str.match(r"^b\d$")].z)) + .8)
    for ax, dos, col in ((axD, "DEL", RED), (axU, "DUP", BLUE)):
        d = T[(T.dosage == dos) & T.stat.str.match(r"^b\d$")].copy(); d["k"] = d.stat.str[1].astype(int); d = d.sort_values("k")
        for _, r in d.iterrows():
            ax.bar(r.k, r.z, .75, color=col, alpha=.3 + .7 * (r.k - 1) / 6)
            s = "**" if r.q < .05 else ("*" if r.P < .05 else "")
            if s:
                ax.text(r.k, r.z + np.sign(r.z) * .35, s, ha="center", va="center", size=6.5)
        ax.axhline(0, color="#666666", lw=.5); ax.set_ylim(-ymax, ymax)
        ax.set_xticks(range(1, 8)); ax.set_xticklabels(BL, rotation=90, size=5.4)
        ax.set_title("Deletions" if dos == "DEL" else "Duplications", size=6.6, weight="bold", pad=2)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
    axU.set_yticklabels([]); axD.set_ylabel("CNV burden, z vs null", size=6.6)
    m = T[(T.stat == "trend") & (T.dosage == "mirror")].iloc[0]
    held = f", {int(m.loo_held)}/{int(m.loo_n)}" if np.isfinite(m.get("loo_held", np.nan)) else ""
    nl = f"\nno recurrent intervals: z = {m.noloci_z:.1f}" if np.isfinite(m.get("noloci_z", np.nan)) else ""
    axD.text(0, 1.12, f"{TL[trait]}: {lab}\nmirror trend z = {m.z:.1f}, P = {fp(m.P)}{held}{nl}", transform=axD.transAxes,
             size=6.8, weight="bold", va="bottom")
    axD.text(1.05, -.36, f"{TL[trait]} TWAS z band", transform=axD.transAxes, ha="center", size=6.4)


def fig_zband():
    fig = plt.figure(figsize=(7.2, 6.2)); W = .135; gap = .012
    for i, t in enumerate(TRAITS):
        for j, (fam, lab) in enumerate(FAMS):
            x0 = .07 + j * .325; y0 = .6 - i * .46
            axD = fig.add_axes([x0, y0, W, .21]); axU = fig.add_axes([x0 + W + gap, y0, W, .21])
            zband_panel(axD, axU, t, fam, lab)
            fig.text(x0 - .05, y0 + .3, "abcdef"[3 * i + j], size=11, weight="bold")
    fig.text(.5, .005, "Each band alone, adjusted for the total number of genes hit; z against 2,000 shifts of the gene annotation along the genome "
             "(UK Biobank).\n* P < 0.05; ** q < 0.05 within the panel. Trend: pre-specified linear trend over band order; mirror = duplication "
             "trend − deletion trend; runs held when each major recurrent locus is removed.", ha="center", va="bottom", size=5.6, color=GREY)
    save(fig, "FigS_twas_zband_convergence")


BINS3 = Path(META3) / "dose_bins" if META3 else Path(os.environ.get("GW2_FIG3_BINS", ROOT.parent / "fig3_dose_bins_2026-11-23"))


def figS_dose_bins():
    """Jon 2026-11-23: 'dose response plots like this [Fig. 2a-h] for ALL of the major strata in new joint model (Fig. 3)'. One figure per trait.
    Rows = CNV size (all CNVs, single-gene CNVs, >= 2-gene CNVs); columns = gene class (TWAS z > +3, z < -3, weak or no TWAS; LC blue, HC red)
    and total genes (grey). x = mean genes deleted (left) / duplicated (right) in each fixed bin; y = adjusted difference in the trait (SD) from
    people with no gene of that class in that CNV-size class, 95% CI; dotted = 34 recurrent intervals removed. Input: fig3_dose_bins.py tables
    (GW2_FIG3_BINS/<trait>/fig3_dose_bins.tsv), or their meta-analysis (GW2_FIG3_META/dose_bins). UNIVARIATE (Jon 2026-11-23: 'This comes
    before the joint model, so these are univariate tests of stratified CNVs'): every curve is its own model, trait ~ covariates + the bins of
    that class in that CNV-size stratum (partition 'uni'). Each panel also prints the univariate DEL z and DUP z of the class (linear gene
    count, conga_grid_qt.py fits uni_<size>; <bins folder>/../uni_<size>/grid/<trait>_x_<trait>/univariate/grid_summary.tsv).
    null model: bin points and 95% CI: none (analytic OLS, descriptive; bins with 1-20 people are not shown). DEL z / DUP z: genome shift of the
    gene annotation, 2,000 shared shifts; * P < 0.05; † also held in all 9 locus-removal runs. Total genes: analytic t (no annotation)."""
    import json
    from scipy.stats import norm
    ROWS = [("N1_inf", "All CNVs", "uni", "uni"), ("N1_1", "Single-gene CNVs", "uni", "uni"), ("N2_inf", "≥2-gene CNVs", "uni", "uni")]

    def uni_z(trait, size, cls):
        """-> {'DEL': (z, P, held), 'DUP': ...} of the univariate fit of this class and stratum, or {}."""
        d = BINS3.parent / f"uni_{size}" / "grid" / f"{trait}_x_{trait}" / "univariate"; f = d / "grid_summary.tsv"
        if not f.exists():
            return {}
        g = pd.read_csv(f, sep="\t"); term = f"{size}~n_all" if cls == "all" else f"{size}~cls:{cls.replace('_', ':')}"; out = {}
        for arm in ("DEL", "DUP"):
            r = g[g.readout == f"{term}|{arm}"]
            if not len(r):
                continue
            r = r.iloc[0]
            if cls == "all":
                t = float(r.t_analytic); out[arm] = (t, 2 * norm.sf(abs(t)), analytic_held(d, f"{term}|{arm}", t))
            else:
                out[arm] = (r.z_emp, r.P_emp, held9(r))
        return out
    COLS = [("up", "TWAS z > +3"), ("down", "TWAS z < −3"), ("other", "weak or no TWAS"), ("all", "total genes")]
    for trait in TRAITS:
        f = BINS3 / trait / "fig3_dose_bins.tsv"
        if not f.exists():
            print("no dose bins for", trait, f); continue
        B = pd.read_csv(f, sep="\t"); B["beta"] = pd.to_numeric(B.beta, errors="coerce"); B["SE"] = pd.to_numeric(B.SE, errors="coerce")
        lab = os.environ.get("GW2_FIG3_LABEL")
        if lab is None:
            mi = Path(META3) / "meta_info.json" if META3 else None
            lab = json.load(open(mi)).get("label", "") if mi and mi.exists() else {"UKBB": "UK Biobank"}.get(str(B.cohort.iloc[0]), str(B.cohort.iloc[0]))
        fig, axs = plt.subplots(3, 4, figsize=(7.2, 6.0), sharex=True)
        for i, (size, slab, pcls, pall) in enumerate(ROWS):
            for j, (cls, clab) in enumerate(COLS):
                ax = axs[i, j]
                sets = [("all", "total genes", GREY, "o")] if cls == "all" else [(f"{cls}_LC", "LC genes", LCC, "o"), (f"{cls}_HC", "HC genes", HCC, "s")]
                P = B[(B.partition == (pall if cls == "all" else pcls)) & (B["size"] == size)]
                for k, (c, l_, col, mk) in enumerate(sets):
                    for rn in ("full", "no_loci"):
                        d = P[(P.cls == c) & (P.run == rn)]; d = d[np.isfinite(d.beta)].sort_values("order")
                        x = np.where(d.side == "DEL", -1, 1) * d.mean_genes.fillna(0) * (1 + (k - (len(sets) - 1) / 2) * .06)
                        ax.errorbar(x, d.beta, yerr=1.96 * d.SE, color=col, marker=mk, ms=3, lw=1.2 if rn == "full" else .8, ls="-" if rn == "full" else ":",
                                    capsize=1.2, mec="white", mew=.4, alpha=1 if rn == "full" else .55, label=l_ if rn == "full" else None)
                ax.set_xscale("symlog", linthresh=1); ax.axhline(0, color=LIGHT, lw=.7, zorder=0); ax.axvline(0, color=LIGHT, lw=.7, zorder=0)
                t = [-30, -10, -3, -1, 0, 1, 3, 10, 30]; ax.set_xlim(-45, 45); ax.set_xticks(t); ax.set_xticklabels([str(abs(v)) for v in t], fontsize=5.4)
                ax.tick_params(axis="y", labelsize=5.6)
                for k, (c, l_, col, mk) in enumerate(sets):       # univariate z of the class: DEL bottom-left, DUP bottom-right
                    zz = uni_z(trait, size, c)
                    for arm, xa, ha in (("DEL", .03, "left"), ("DUP", .97, "right")):
                        if arm in zz and np.isfinite(zz[arm][0]):
                            z_, p_, h_ = zz[arm]
                            ax.text(xa, .03 + .085 * (len(sets) - 1 - k), f"z {z_:+.1f}" + ("*" if p_ < .05 else "") + ("†" if h_ else ""), transform=ax.transAxes,
                                    ha=ha, va="bottom", size=5.2, color=col, weight="bold" if p_ < .05 else "normal")
                if i == 0:
                    ax.set_title(clab, fontsize=7.4, weight="bold")
                    if j in (0, 3):
                        ax.legend(frameon=False, fontsize=5.4, loc="best", handlelength=1.4)
                if j == 0:
                    ax.set_ylabel(f"{slab}\n{TL[trait]} (SD), adjusted", fontsize=6.4)
                if i == 2:
                    ax.set_xlabel("genes deleted     genes duplicated", fontsize=5.8)
        fig.suptitle(f"{TL[trait]}: univariate dose response of each gene class, by CNV size", x=.02, ha="left", fontsize=8.6, weight="bold")
        if lab:
            fig.text(.99, .992, lab, ha="right", va="top", size=6.4, color=GREY)
        fig.text(.02, .006, "points: difference from people with no gene of that class in that CNV-size stratum, each class in its own model (95% CI, analytic);  "
                 "dotted: 34 recurrent intervals removed;  bins with 1–20 people not shown\n"
                 "z: univariate z of the class's gene count per dosage against the genome-shift null (total genes: analytic t);  * P < 0.05;  "
                 "† P < 0.05 and held when each of 8 major loci, or all minor loci, is removed", size=5.0, color=INK, va="bottom")
        fig.tight_layout(rect=(0, .035, 1, .965))
        save(fig, f"FigS_fig3_dose_response_{trait}" + ("_meta" if META3 else ""))


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=300)
    plt.close(fig); print("written", OUT / f"{name}.png")


if __name__ == "__main__":
    want = sys.argv[1:] or ["1", "2", "3"]
    if "1" in want: fig1()
    if "2" in want: fig2()
    if "3" in want: fig3()
    if "D" in want: fig3(dose=True)   # alternate Fig. 3, dose-aligned axis (Jon 2026-11-23)
    if "B" in want: figS_dose_bins()  # binned dose response of every Fig. 3 stratum (Jon 2026-11-23)
    if "S" in want: figS_heatmaps()
    if "Z" in want: fig_zband()
    if "X" in want: figS_cross_mirror()
    if "F" in want: figS_full_grid()
