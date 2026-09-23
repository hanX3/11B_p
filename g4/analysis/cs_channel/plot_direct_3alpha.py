#!/usr/bin/env python3
"""Direct 3-alpha (democratic phase space) figure set, resonance-agnostic.

Direct decay is the same democratic phase-space process regardless of which
resonance seeded it, so this script does NOT distinguish 162 from 675: it
selects every direct-3-alpha event (reaction_channel == 2) and pools them.

Produces three single 4:3 panels:

  <stem>_spectrum.png : pooled single-alpha LAB energy spectrum
      (all three alphas; smooth phase-space distribution, no
      alpha0/alpha1 peak structure).
  <stem>_angular.png  : cos(theta) of alpha1 vs the beam axis,
      3-alpha c.m. (flat) overlaid with the LAB frame (small forward
      tilt from the c.m. boost).
  <stem>_opening.png  : pairwise c.m. opening-angle correlation
      (12/13/23 pooled) -- pure phase-space + momentum-conservation
      baseline, peak near 130-140 deg.

Dual use:
    terminal :  python3 plot_direct_3alpha.py FILE.root
                python3 plot_direct_3alpha.py FILE.root --stem d3a
                (saves the three PNGs, no window pop-up)
    jupyter  :  from plot_direct_3alpha import plot_direct_3alpha
                fig_spec, fig_ang, fig_open = plot_direct_3alpha("FILE.root")
                (figures display inline automatically under %matplotlib inline)
"""

import numpy as np
import matplotlib.pyplot as plt
import uproot

plt.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 15,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "axes.linewidth": 1.1,
    "legend.fontsize": 12,
    "figure.dpi": 100,
})

FIGSIZE = (4, 4)
FIXED_MARGINS = dict(left=0.18, right=0.86, bottom=0.15, top=0.90)

# reaction_channel == 2 is DirectDecay3Alpha (both 162- and 675-seeded direct
# decay); this excludes the gamma-capture channel, which also uses branch_id=-1.
DIRECT_CHANNEL = 2

TREE = "reaction"

FIELDS = [
    "e_alpha1", "e_alpha2", "e_alpha3",
    "theta_lab_alpha1", "cos_theta_primary_cm",
    "opening_angle_alpha12_cm", "opening_angle_alpha13_cm",
    "opening_angle_alpha23_cm",
    "reaction_channel",
]

C_LAB = "#0E7C6B"   # teal
C_CM  = "#C77D2E"   # orange


def clean_cos(c):
    c = np.asarray(c, dtype=float)
    return c[np.isfinite(c) & (c >= -1.0) & (c <= 1.0)]


def plot_direct_3alpha(filename, save=False, stem="direct3alpha"):
    """Make the three direct-3-alpha figures.

    Parameters
    ----------
    filename : str or pathlib.Path
        ROOT file containing the ``reaction`` tree.
    save : bool, optional
        False (default) -> do not write files.
        True            -> write {stem}_spectrum.png, {stem}_angular.png,
                           {stem}_opening.png.
    stem : str, optional
        Output filename stem (default "direct3alpha").

    Returns
    -------
    (fig_spectrum, fig_angular, fig_opening) : tuple of matplotlib Figures
    """
    with uproot.open(filename) as f:
        a = f[TREE].arrays(FIELDS, library="np")

    m = np.asarray(a["reaction_channel"]) == DIRECT_CHANNEL
    n = int(m.sum())
    if n == 0:
        raise ValueError(
            "no direct-decay events (reaction_channel == 2) found. "
            "Did you run a direct-3-alpha macro?")

    def col(k):
        v = np.asarray(a[k], dtype=float)[m]
        return v[np.isfinite(v)]

    # ============ figure 1: pooled LAB energy spectrum =================
    e = np.concatenate([col("e_alpha1"), col("e_alpha2"), col("e_alpha3")])
    e = e[e > 0]
    emax = e.max() * 1.05
    ebins = np.linspace(0, emax, 120)

    fig_spec, ax = plt.subplots(figsize=FIGSIZE)
    ax.hist(e, bins=ebins, histtype="step", lw=2.2, color=C_LAB)
    ax.set_xlabel(r"$\alpha$ energy, LAB (MeV)")
    ax.set_ylabel("counts")
    ax.set_xlim(0, emax)
    ax.set_title(r"direct $3\alpha$ - pooled LAB spectrum")
    fig_spec.subplots_adjust(**FIXED_MARGINS)

    # ============ figure 2: angular, c.m. vs lab =======================
    abins = np.linspace(-1, 1, 41)
    actr = 0.5 * (abins[:-1] + abins[1:])

    cos_cm = clean_cos(np.asarray(a["cos_theta_primary_cm"])[m])
    cos_lab = clean_cos(np.cos(np.asarray(a["theta_lab_alpha1"])[m]))
    h_cm, _ = np.histogram(cos_cm, bins=abins, density=True)
    h_lab, _ = np.histogram(cos_lab, bins=abins, density=True)

    fig_ang, ax = plt.subplots(figsize=FIGSIZE)
    ax.step(actr, h_cm, where="mid", lw=2.4, color=C_CM,
            label=r"c.m. ($3\alpha$ frame)")
    ax.step(actr, h_lab, where="mid", lw=2.4, color=C_LAB, label="lab")
    ax.set_xlabel(r"$\cos\theta$  ($\alpha_1$ vs beam)")
    ax.set_ylabel("normalised yield")
    ax.set_xlim(-1, 1)
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, loc="upper left")
    ax.set_title("angular vs beam: c.m. vs lab")
    fig_ang.subplots_adjust(**FIXED_MARGINS)

    # ============ figure 3: pairwise opening-angle correlation =========
    opening = np.degrees(np.concatenate([
        col("opening_angle_alpha12_cm"),
        col("opening_angle_alpha13_cm"),
        col("opening_angle_alpha23_cm")]))
    obins = np.linspace(0, 180, 61)
    octr = 0.5 * (obins[:-1] + obins[1:])
    h_op, _ = np.histogram(opening, bins=obins, density=True)

    fig_open, ax = plt.subplots(figsize=FIGSIZE)
    ax.step(octr, h_op, where="mid", lw=2.4, color=C_CM)
    ax.set_xlabel(r"pairwise opening angle $\theta_{ij}$ (deg, c.m.)")
    ax.set_ylabel("normalised yield")
    ax.set_xlim(0, 180)
    ax.set_ylim(bottom=0)
    ax.set_title(r"direct $3\alpha$ opening-angle correlation")
    fig_open.subplots_adjust(**FIXED_MARGINS)

    # ------------------------------------------------------------------- save
    if save:
        fig_spec.savefig(f"{stem}_spectrum.png", dpi=200)
        fig_ang.savefig(f"{stem}_angular.png", dpi=200)
        fig_open.savefig(f"{stem}_opening.png", dpi=200)
        print(f"wrote {stem}_spectrum.png")
        print(f"wrote {stem}_angular.png")
        print(f"wrote {stem}_opening.png")
        print(f"  direct-decay events: {n}")

    return fig_spec, fig_ang, fig_open


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Direct 3-alpha (phase space) figure set, "
                    "resonance-agnostic.")
    parser.add_argument("rootfile", help="ROOT file with the 'reaction' tree")
    parser.add_argument("--stem", default="direct3alpha",
                        help="output filename stem (default: direct3alpha)")
    args = parser.parse_args()

    plot_direct_3alpha(args.rootfile, save=True, stem=args.stem)
