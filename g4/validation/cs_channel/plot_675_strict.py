#!/usr/bin/env python3
"""675-keV strict-model (symmetrized coherent l1 = 1/3) figure set.

Reads ONE strict-run reaction file (validation_675_primary_spectrum.mac)
and produces the four slide figures, each a single 4:3 panel sized for
a 2-row x 3-column grid on a 16:8 PowerPoint slide:

  <stem>_spectrum.png : symmetrized single-alpha LAB spectrum with the
      three region annotations ('primary-like' + forward secondaries /
      backward secondaries / sideways secondaries). The regions map
      one-to-one onto the three peaks of the cos(chi) correlation.
  <stem>_primary_angular.png : cos(theta) of alpha slot 1 vs the beam
      axis (c.m.) -- expect FLAT (equal-m sum, rotational invariance).
  <stem>_coschi.png : secondary-alpha correlation cos(chi) vs the 8Be
      recoil axis -- the three-peak structure carved by the l1 = 1/3
      coherent mixing (isotropic reference at 0.5). Note: (2,3) is not
      always the resonant pair after symmetrization; read the shape
      qualitatively -- the clean observable is the Dalitz density.
  <stem>_dalitz.png : Dalitz plot in Kuhlwein 2022 coordinates,
      x = sqrt(3)(E2-E3)/Sum E, y = (2E1-E2-E3)/Sum E, energies in the
      3-alpha c.m. frame, ALL SIX permutations filled per event
      (identical bosons). Empty bins are white; the kinematic boundary
      (unit circle) is dashed.

Selection: branch_id == 1 (sequential channel).

Usage:
    python3 plot_675_strict.py FILE.root
    python3 plot_675_strict.py FILE.root --stem s675 --nbins 160
"""

import argparse
import sys
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

FIGSIZE = (16 / 3, 8 / 2)   # 4:3

FIELDS = [
    "e_alpha1", "e_alpha2", "e_alpha3",
    "cos_theta_primary_cm",
    "cos_chi_secondary_8be", "cos_theta_secondary_correlation",
    "e_3alpha_cm_alpha1", "e_3alpha_cm_alpha2", "e_3alpha_cm_alpha3",
    "branch_id",
]

C_MAIN = "#990000"   # 675 red
LIM = 1.15           # Dalitz axis range


def clean_cos(c):
    c = np.asarray(c, dtype=float)
    return c[np.isfinite(c) & (c >= -1.0) & (c <= 1.0)]


def dalitz_xy_symmetrized(E):
    """E: (N, 3) c.m. kinetic energies -> x, y with all 6 permutations."""
    perms = [(0, 1, 2), (0, 2, 1), (1, 0, 2),
             (1, 2, 0), (2, 0, 1), (2, 1, 0)]
    s = E.sum(axis=1)
    ok = np.isfinite(s) & (s > 0)
    E, s = E[ok], s[ok]
    xs, ys = [], []
    for i, j, k in perms:
        xs.append(np.sqrt(3.0) * (E[:, j] - E[:, k]) / s)
        ys.append((2.0 * E[:, i] - E[:, j] - E[:, k]) / s)
    return np.concatenate(xs), np.concatenate(ys)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rootfile")
    ap.add_argument("--stem", default="alpha675")
    ap.add_argument("--nbins", type=int, default=200,
                    help="Dalitz 2D bins per axis (default 200)")
    args = ap.parse_args()

    with uproot.open(args.rootfile) as f:
        a = f["tr"].arrays(FIELDS, library="np")

    m = np.asarray(a["branch_id"]) == 1
    n = int(m.sum())
    if n == 0:
        sys.exit("no sequential events (branch_id == 1) found.\n"
                 "  Did you run validation_675_primary_spectrum.mac?")

    def col(k):
        v = np.asarray(a[k], dtype=float)[m]
        return v[np.isfinite(v)]

    # ============ figure 1: symmetrized single-alpha spectrum ==========
    e1 = col("e_alpha1")
    e1 = e1[e1 > 0]
    emax = e1.max() * 1.05
    ebins = np.linspace(0, emax, 120)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.hist(e1, bins=ebins, histtype="step", lw=2.2, color=C_MAIN)
    ax.set_xlabel(r"$\alpha$ energy, LAB (MeV)")
    ax.set_ylabel("counts")
    ax.set_xlim(0, emax)
    ax.set_title(r"single-$\alpha$ spectrum (symmetrized)")
    ax.text(0.32, 0.88, "'primary-like'\n+ forward secondaries",
            transform=ax.transAxes, fontsize=11, color="#555555",
            ha="center")
    ax.text(0.13, 0.48, "backward\nsecondaries",
            transform=ax.transAxes, fontsize=11, color="#555555",
            ha="center")
    ax.text(0.38, 0.05, "sideways\nsecondaries",
            transform=ax.transAxes, fontsize=10, color="#999999",
            ha="center")
    plt.tight_layout(pad=0.8)
    plt.savefig(f"{args.stem}_spectrum.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {args.stem}_spectrum.png")

    # ============ figure 2: primary angular (expect flat) ===============
    abins = np.linspace(-1, 1, 41)
    actr = 0.5 * (abins[:-1] + abins[1:])

    c = clean_cos(np.asarray(a["cos_theta_primary_cm"])[m])
    h, _ = np.histogram(c, bins=abins, density=True)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.step(actr, h, where="mid", lw=2.4, color=C_MAIN)
    ax.set_xlabel(r"$\cos\theta$  ($\alpha$ vs beam, c.m.)")
    ax.set_ylabel("normalised yield")
    ax.set_xlim(-1, 1)
    ax.set_ylim(0, 1)
    ax.set_title("angular vs beam (expect isotropic)")
    plt.tight_layout(pad=0.8)
    plt.savefig(f"{args.stem}_primary_angular.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {args.stem}_primary_angular.png")

    # ============ figure 3: cos(chi) secondary correlation ==============
    chi_a = clean_cos(np.asarray(a["cos_chi_secondary_8be"])[m])
    chi_b = clean_cos(np.asarray(a["cos_theta_secondary_correlation"])[m])
    chi = chi_b if chi_b.size >= chi_a.size else chi_a
    hh, _ = np.histogram(chi, bins=abins, density=True)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.step(actr, hh, where="mid", lw=2.4, color=C_MAIN)
    ax.set_xlabel(r"$\cos\chi$  (secondary $\alpha$ vs $^8$Be axis)")
    ax.set_ylabel("normalised yield")
    ax.set_xlim(-1, 1)
    ax.set_ylim(bottom=0)
    ax.set_title(r"internal correlation ($l_1$ = 1, 3 coherent)")
    plt.tight_layout(pad=0.8)
    plt.savefig(f"{args.stem}_coschi.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {args.stem}_coschi.png")

    # ============ figure 4: symmetrized Dalitz plot ======================
    E = np.column_stack([col("e_3alpha_cm_alpha1"),
                         col("e_3alpha_cm_alpha2"),
                         col("e_3alpha_cm_alpha3")])
    x, y = dalitz_xy_symmetrized(E)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    edges = np.linspace(-LIM, LIM, args.nbins + 1)
    hist2d, _, _ = np.histogram2d(x, y, bins=[edges, edges])
    hist2d = np.ma.masked_equal(hist2d, 0)      # empty bins -> white
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("white")
    im = ax.pcolormesh(edges, edges, hist2d.T, cmap=cmap, rasterized=True)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(t), np.sin(t), color="#CCCCCC", lw=1.0, ls="--")
    ax.set_xlabel(r"$x=\sqrt{3}\,(E_2-E_3)/\sum E_i$")
    ax.set_ylabel(r"$y=(2E_1-E_2-E_3)/\sum E_i$")
    ax.set_xlim(-LIM, LIM)
    ax.set_ylim(-LIM, LIM)
    ax.set_aspect("equal")
    ax.set_title(r"Dalitz plot ($l_1$ = 1, 3 coherent)")
    plt.tight_layout(pad=0.8)
    plt.savefig(f"{args.stem}_dalitz.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {args.stem}_dalitz.png")

    print(f"  sequential events: {n}")


if __name__ == "__main__":
    main()
