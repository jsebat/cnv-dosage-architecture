# Checks for the person-level ConGA port (conga_ctx_qt.py, run_pipeline.py step 'ctx', export, meta_ctx.py).
# null model: n/a (test assertions; the runs under test use the genome-shift null with the shared shifts)
"""Subcommands
  pick    --ukbb_full <UKBB full_<t>_<d> dirs, comma> --sets <membership> --n 3 --out <subset membership + set list>
          (the n smallest-P signal groups and the n smallest-P contexts of the UKBB tables, for the reproduction check)
  repro   --ukbb <UKBB run dirs, comma> --port <ported run dir> --mode full|loo   (z, z_emp, P_emp equal to the UKBB run)
  export  --results <pipeline results ctx dir> --export <export ctx dir>          (no carrier counts, float32 nulls, same z)
  meta    --a <export ctx dir A> --aa <meta A+A dir> --ab <meta A+B dir>          (A+A: z = sqrt2 z_A, z_emp = z_emp_A; A+B finite)
Each prints PASS/FAIL lines and exits non-zero on any FAIL."""
import argparse
import glob
import os
import sys
import numpy as np
import pandas as pd

FAIL = []


def check(ok, msg):
    print(("PASS " if ok else "FAIL ") + msg, flush=True)
    if not ok:
        FAIL.append(msg)


def tab(d):
    f = glob.glob(os.path.join(d, "conga_ctx_*.tsv.gz"))
    return pd.read_csv(f[0], sep="\t") if f else None


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("cmd")
    for k in ["ukbb_full", "sets", "out", "ukbb", "port", "mode", "results", "export", "a", "aa", "ab"]:
        ap.add_argument(f"--{k}", default="")
    ap.add_argument("--n", type=int, default=3)
    a = ap.parse_args()
    if a.cmd == "pick":
        T = pd.concat([tab(d) for d in a.ukbb_full.split(",")])
        T = T[T.run == "full"].sort_values("P_emp")
        sg = T[T.set_id.str.startswith("SG:")].set_id.drop_duplicates().head(a.n)
        cx = T[~T.set_id.str.startswith("SG:")].set_id.drop_duplicates().head(a.n)
        pick = list(sg) + list(cx)
        sm = pd.read_csv(a.sets, sep="\t", dtype=str); sm = sm[sm.set_id.isin(pick)]
        os.makedirs(a.out, exist_ok=True); sm.to_csv(os.path.join(a.out, "membership.tsv.gz"), sep="\t", index=False)
        pd.Series(pick).to_csv(os.path.join(a.out, "sets.txt"), index=False, header=False)
        print("picked", pick); return
    if a.cmd == "repro":
        P = tab(a.port); U = pd.concat([t for t in (tab(d) for d in a.ukbb.split(",")) if t is not None])
        m = P.merge(U, on=["run", "set_id"], suffixes=("_port", "_ukbb"))
        check(len(m) == len(P) and len(m) > 0, f"{a.mode}: every ported set x run found in the UKBB run ({len(m)} of {len(P)})")
        for c, tol in [("z", 2e-3), ("z_emp", 2e-3), ("P_emp", 1e-9), ("n_genes", 0)]:
            dmax = float(np.nanmax(np.abs(m[f"{c}_port"].astype(float) - m[f"{c}_ukbb"].astype(float)))) if len(m) else np.inf
            check(dmax <= tol, f"{a.mode}: max |{c} port - UKBB| = {dmax:.3g} (tolerance {tol}; tables are written to 5 significant digits)")
        print(m[["run", "set_id", "z_port", "z_ukbb", "z_emp_port", "z_emp_ukbb", "P_emp_port", "P_emp_ukbb"]].head(12).to_string(index=False))
    if a.cmd == "export":
        for sub in sorted(os.listdir(a.results)):
            r, e = os.path.join(a.results, sub), os.path.join(a.export, sub)
            if not os.path.exists(os.path.join(r, "info.json")) or tab(r) is None:
                continue
            R, E = tab(r), tab(e)
            check(E is not None and "n_carriers" not in E.columns, f"export {sub}: table present, no carrier count column")
            check(np.allclose(R.z.to_numpy(float), E.z.to_numpy(float), equal_nan=True), f"export {sub}: z unchanged")
            for f in glob.glob(os.path.join(r, "null_z_*.npy")):
                x = np.load(os.path.join(e, os.path.basename(f)))
                check(x.dtype == np.float32 and x.shape == np.load(f).shape, f"export {sub}/{os.path.basename(f)}: float32, same shape")
    if a.cmd == "meta":
        for sub in sorted(os.listdir(a.aa)):
            if not os.path.isdir(os.path.join(a.aa, sub)):
                continue
            A, AA, AB = tab(os.path.join(a.a, sub)), tab(os.path.join(a.aa, sub)), tab(os.path.join(a.ab, sub))
            m = AA.merge(A, on=["run", "set_id"], suffixes=("_aa", "_a"))
            check(np.allclose(m.z_aa, np.sqrt(2) * m.z_a, atol=1e-3, equal_nan=True), f"meta A+A {sub}: z = sqrt(2) z_A")
            check(np.allclose(m.z_emp_aa, m.z_emp_a, atol=2e-3, equal_nan=True), f"meta A+A {sub}: z_emp = z_emp_A (matched null)")
            check(bool(AB.shared_shifts.all()) and np.isfinite(AB.z).mean() > 0.9, f"meta A+B {sub}: shared shifts, z finite")
    if FAIL:
        sys.exit(f"{len(FAIL)} check(s) failed")


if __name__ == "__main__":
    main()
