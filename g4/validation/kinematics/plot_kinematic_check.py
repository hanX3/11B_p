#!/usr/bin/env python3
"""Kinematics validation figure: simulated (theta_lab, E) of the primary
alpha from the 162-keV alpha0 branch, overlaid with the ANALYTIC
relativistic two-body curve p + 11B -> alpha + 8Be(g.s.).

The 8Be ground state is narrow (~eV), so the primary alpha of the
alpha0 branch must lie exactly on the two-body kinematic curve --
a one-figure proof that the four-vector bookkeeping (invariants,
c.m. split, boost) is exact.

Usage:
    python3 plot_kinematic_check.py FILE.root
    python3 plot_kinematic_check.py FILE.root --ep 0.162 --stem kincheck
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import uproot

plt.rcParams.update({
    "font.size": 14, "axes.titlesize": 15, "axes.labelsize": 14,
    "xtick.labelsize": 12, "ytick.labelsize": 12,
    "axes.linewidth": 1.1, "legend.fontsize": 12, "figure.dpi": 100,
})

FIGSIZE = (16 / 3, 8 / 2)

M_P = 938.272
M_B11 = 10252.548
M_ALPHA = 3727.379
M_BE8 = 7454.850

FIELDS = ["e_alpha1", "theta_lab_alpha1", "branch_id"]


def analytic_curve(ep_lab_mev):
    Ep = M_P + ep_lab_mev
    pp = np.sqrt(Ep ** 2 - M_P ** 2)
    s = (Ep + M_B11) ** 2 - pp ** 2
    rs = np.sqrt(s)
    beta = pp / (Ep + M_B11)
    gamma = 1.0 / np.sqrt(1.0 - beta ** 2)
    e_star = (s + M_ALPHA ** 2 - M_BE8 ** 2) / (2.0 * rs)
    p_star = np.sqrt(np.maximum(e_star ** 2 - M_ALPHA ** 2, 0.0))
    ct = np.linspace(-1.0, 1.0, 800)
    st = np.sqrt(1.0 - ct ** 2)
    e_lab = gamma * (e_star + beta * p_star * ct)
    pz = gamma * (p_star * ct + beta * e_star)
    pt = p_star * st
    return np.degrees(np.arctan2(pt, pz)), e_lab - M_ALPHA


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rootfile")
    ap.add_argument("--ep", type=float, default=0.162,
                    help="proton lab kinetic energy in MeV (default 0.162)")
    ap.add_argument("--stem", default="kincheck")
    args = ap.parse_args()

    with uproot.open(args.rootfile) as f:
        a = f["tr"].arrays(FIELDS, library="np")

    m = np.asarray(a["branch_id"]) == 0
    if m.sum() == 0:
        sys.exit("no alpha0 events (branch_id == 0) found.")

    e = np.asarray(a["e_alpha1"], dtype=float)[m]
    th = np.degrees(np.asarray(a["theta_lab_alpha1"], dtype=float)[m])
    ok = np.isfinite(e) & np.isfinite(th) & (e > 0)
    e, th = e[ok], th[ok]

    fig, ax = plt.subplots(figsize=FIGSIZE)
    hb = ax.hist2d(th, e, bins=[90, 90],
                   range=[[0, 180], [5.2, 6.5]],
                   cmin=1, cmap="viridis", rasterized=True)
    fig.colorbar(hb[3], ax=ax, fraction=0.046, pad=0.03)

    tc, ec = analytic_curve(args.ep)
    ax.plot(tc, ec, "-", lw=3.2, color="white", zorder=5)
    ax.plot(tc, ec, "--", lw=1.4, color="#C00000", zorder=6)

    ax.set_xlabel(r"$\theta_{\rm lab}$ (deg)")
    ax.set_ylabel(r"$E_{\alpha}$ LAB (MeV)")
    ax.set_xlim(0, 180)
    ax.set_ylim(5.2, 6.5)
    ax.set_title(r"$\alpha_0$ primary: exact kinematics check")
    plt.tight_layout(pad=0.8)
    plt.savefig(f"{args.stem}.png", dpi=200, bbox_inches="tight")
    print(f"wrote {args.stem}.png   ({int(m.sum())} alpha0 events)")


if __name__ == "__main__":
    main()
