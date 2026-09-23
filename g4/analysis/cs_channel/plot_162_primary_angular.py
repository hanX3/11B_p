#!/usr/bin/env python3
"""Lab-frame primary-alpha angular distribution: anisotropic alpha0 vs
isotropic alpha1, from a single run.

One figure, lab frame only. Two curves, split by branch_id:

    branch_id == 0  ->  alpha0  (8Be g.s.)  : anisotropic (Becker a2=0.6)
    branch_id == 1  ->  alpha1  (8Be 2+)    : isotropic   (Spraker, a2=0)

Feed a run with the 162 primary angular distribution enabled, so that the
alpha0 channel carries its Becker anisotropy while alpha1 stays isotropic.
The angle is cos(theta_lab) of the primary alpha; in the lab frame the
distribution is boosted, so it is NOT the analytic W(theta) of the c.m. frame.

Dual use:
    terminal :  python3 plot_162_primary_angular.py FILE.root
                (saves the PNG, no window pop-up)
    jupyter  :  from plot_162_primary_angular import plot_162_primary_angular
                fig = plot_angular_alpha0("FILE.root")
                (displays inline automatically under %matplotlib inline)
"""

import numpy as np
import matplotlib.pyplot as plt
import uproot

FIXED_MARGINS = dict(left=0.18, right=0.86, bottom=0.15, top=0.90)

# tree name -- change if your file uses a different tree name
TREE = "reaction"

C_A0 = "#0E7C6B"   # teal   = alpha0 (anisotropic)
C_A1 = "#C77D2E"   # orange = alpha1 (isotropic)


def plot_162_primary_angular(filename, save=False):
    """Plot the lab-frame primary-alpha angular distribution for both channels.

    Parameters
    ----------
    filename : str or pathlib.Path
        ROOT file from a run with the 162 primary angular distribution enabled.
    save : bool or str, optional
        False (default) -> do not write a file.
        True            -> write the default PNG.
        str             -> use it as the output path.

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    with uproot.open(filename) as f:
        a = f[TREE].arrays(["theta_lab_alpha1", "branch_id"], library="np")

    bid = a["branch_id"]
    cos_lab = np.cos(a["theta_lab_alpha1"])   # theta_lab is in radians

    def clean(mask):
        c = cos_lab[mask]
        return c[np.isfinite(c) & (c >= -1) & (c <= 1)]

    a0 = clean(bid == 0)   # anisotropic
    a1 = clean(bid == 1)   # isotropic

    bins = np.linspace(-1, 1, 41)
    ctr = 0.5 * (bins[:-1] + bins[1:])

    def hist(c):
        h, _ = np.histogram(c, bins=bins, density=True)
        return h

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    ax.axhline(0.5, color="#999999", lw=0.8, ls=":")   # isotropic reference
    ax.step(ctr, hist(a0), where="mid", lw=1.8, color=C_A0,
            label=rf"$\alpha_0$ anisotropic (8Be g.s.), N={a0.size}")
    ax.step(ctr, hist(a1), where="mid", lw=1.8, color=C_A1,
            label=rf"$\alpha_1$ isotropic (8Be 2$^+$), N={a1.size}")
    ax.set_xlabel(r"$\cos\theta_{\mathrm{lab}}$  (primary $\alpha$)")
    ax.set_ylabel("normalised yield")
    ax.set_xlim(-1, 1)
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=9)
    fig.subplots_adjust(**FIXED_MARGINS)

    if save:
        path = "angular_alpha0_alpha1_lab.png" if save is True else str(save)
        fig.savefig(path, dpi=200)
        print(f"wrote {path}")

    return fig


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Lab-frame primary-alpha angular distribution: "
                    "anisotropic alpha0 vs isotropic alpha1.")
    parser.add_argument("filename", help="ROOT file (162 primary angular on)")
    args = parser.parse_args()

    plot_162_primary_angular(args.filename, save=True)
