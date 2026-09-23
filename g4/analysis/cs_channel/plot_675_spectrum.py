#!/usr/bin/env python3
"""Single-alpha LAB energy spectrum, decomposed into the primary alpha and
the secondary alphas (sequential channel, branch_id == 1).

The three alpha slots are tagged by role:
    alpha_role == 1  (FirstStepAlpha)  -> primary   (1 per event)
    alpha_role == 2  (BeDecayAlpha)     -> secondary (2 per event)
The figure overlays: total (all three), primary only, secondary only.

IMPORTANT -- permutation symmetrization must be OFF for this split to mean
anything. With /h11b/675StrictPermutationSymmetrized true (the default of
validation_675_seq.mac) the three slots are interchangeable
bookkeeping labels, so all three marginal distributions are identical and
the "primary" and "secondary" curves come out with the same shape. Generate
the input with validation_675_seq_no_symmetrized.mac to get a physical
primary-vs-secondary decomposition.

Dual use:
    terminal :  python3 plot_675_spectrum.py FILE.root
                python3 plot_675_spectrum.py FILE.root --out my_name.png
                (saves the PNG)
    jupyter  :  from plot_675_spectrum import plot_675_spectrum
                fig = plot_675_spectrum("FILE.root")
                fig = plot_675_spectrum("FILE.root", save="my_name.png")
                (figure displays inline automatically under %matplotlib inline)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import uproot

plt.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 15,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "axes.linewidth": 1.1,
    "legend.fontsize": 11,
    "figure.dpi": 100,
})

FIGSIZE = (16 / 3, 8 / 2)   # 4:3
FIXED_MARGINS = dict(left=0.18, right=0.86, bottom=0.15, top=0.90)

TREE = "reaction"

# particle-role codes (include/H11BParticleLabel.hh)
FIRST_STEP_ALPHA = 1   # primary
BE_DECAY_ALPHA   = 2   # secondary

FIELDS = [
    "e_alpha1", "e_alpha2", "e_alpha3",
    "alpha1_role", "alpha2_role", "alpha3_role",
    "branch_id",
]

C_TOTAL = "#333333"   # dark  = total
C_PRIM  = "#0E7C6B"   # teal  = primary
C_SEC   = "#C77D2E"   # orange = secondary


def plot_675_spectrum(rootfile, save=False, nbins=120):
    """Single-alpha LAB spectrum split into primary and secondary alphas.

    Parameters
    ----------
    rootfile : str or pathlib.Path
        ROOT file containing the ``reaction`` tree.
    save : bool or str, optional
        False (default) -> do not write a file.
        True            -> write spectrum_primary_secondary.png.
        str             -> write to that path (".png" appended if missing).
    nbins : int, optional
        Number of energy bins (default 120).

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    with uproot.open(rootfile) as f:
        a = f[TREE].arrays(FIELDS, library="np")

    m = np.asarray(a["branch_id"]) == 1
    n = int(m.sum())
    if n == 0:
        raise ValueError(
            "no sequential events (branch_id == 1) found. "
            "Did you run a 675 sequential macro?")

    e_slots    = [np.asarray(a["e_alpha1"], dtype=float)[m],
                  np.asarray(a["e_alpha2"], dtype=float)[m],
                  np.asarray(a["e_alpha3"], dtype=float)[m]]
    role_slots = [np.asarray(a["alpha1_role"])[m],
                  np.asarray(a["alpha2_role"])[m],
                  np.asarray(a["alpha3_role"])[m]]

    prim_parts, sec_parts = [], []
    for e, role in zip(e_slots, role_slots):
        good = np.isfinite(e) & (e > 0)
        prim_parts.append(e[good & (role == FIRST_STEP_ALPHA)])
        sec_parts.append(e[good & (role == BE_DECAY_ALPHA)])

    e_prim = np.concatenate(prim_parts)
    e_sec  = np.concatenate(sec_parts)
    e_all  = np.concatenate([e_prim, e_sec])

    emax = e_all.max() * 1.05
    ebins = np.linspace(0, emax, nbins)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.hist(e_all,  bins=ebins, histtype="step", lw=2.4, color=C_TOTAL)
    ax.hist(e_prim, bins=ebins, histtype="step", lw=2.0, color=C_PRIM)
    ax.hist(e_sec,  bins=ebins, histtype="step", lw=2.0, color=C_SEC)

    handles = [
        Line2D([], [], color=C_TOTAL, lw=2.4, label=f"total"),
        Line2D([], [], color=C_PRIM,  lw=2.0,
               label=rf"primary $\alpha$"),
        Line2D([], [], color=C_SEC,   lw=2.0,
               label=rf"secondary $\alpha$"),
    ]
    ax.legend(handles=handles, loc="upper left", frameon=False)

    ax.set_xlabel(r"$\alpha$ energy, LAB (MeV)")
    ax.set_ylabel("counts")
    ax.set_xlim(0, emax)
    fig.subplots_adjust(**FIXED_MARGINS)

    if save:
        path = ("spectrum_primary_secondary.png" if save is True
                else str(save))
        if not path.lower().endswith(".png"):
            path += ".png"
        fig.savefig(path, dpi=200)
        print(f"wrote {path}   (sequential events: {n})")

    return fig


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("rootfile")
    ap.add_argument("--out", default=None,
                    help="output png name (default: "
                         "spectrum_primary_secondary.png)")
    ap.add_argument("--nbins", type=int, default=120)
    args = ap.parse_args()

    out = True if args.out is None else args.out
    plot_675_spectrum(args.rootfile, save=out, nbins=args.nbins)
