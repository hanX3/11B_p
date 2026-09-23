#!/usr/bin/env python3
"""675-keV direct 3-alpha (democratic phase space) figure set.

Reads the output of validation_675_direct3alpha.mac and produces four
single 4:3 panels (2-row x 3-column grid on a 16:8 slide):

  <stem>_spectrum.png : pooled single-alpha LAB energy spectrum
      (all three alphas; smooth phase-space distribution, no
      alpha0/alpha1 peak structure).
  <stem>_angular.png : cos(theta) of alpha1 vs the beam axis,
      3-alpha c.m. (flat) overlaid with the LAB frame (small forward
      tilt from the c.m. boost).
  <stem>_opening.png : pairwise c.m. opening-angle correlation
      (12/13/23 pooled) -- pure phase-space + momentum-conservation
      baseline, peak near 130-140 deg.
  <stem>_dalitz.png : Dalitz plot, Kuhlwein coordinates, all six
      permutations filled -- uniform disk (no structure): the
      no-mechanism baseline against the sequential-strict pattern.

Selection: branch_id == -1 (direct decay).

Usage:
    python3 plot_675_direct.py FILE.root
    python3 plot_675_direct.py FILE.root --stem d675 --nbins 160
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
    "theta_lab_alpha1", "cos_theta_primary_cm",
    "opening_angle_alpha12_cm", "opening_angle_alpha13_cm",
    "opening_angle_alpha23_cm",
    "e_3alpha_cm_alpha1", "e_3alpha_cm_alpha2", "e_3alpha_cm_alpha3",
    "branch_id",
]

C_LAB = "#0E7C6B"   # teal
C_CM  = "#C77D2E"   # orange
LIM = 1.15          # Dalitz axis range


def clean_cos(c):
    c = np.asarray(c, dtype=float)
    return c[np.isfinite(c) & (c >= -1.0) & (c <= 1.0)]


def dalitz_xy_symmetrized(E):
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
    ap.add_argument("--stem", default="direct675")
    ap.add_argument("--nbins", type=int, default=200,
                    help="Dalitz 2D bins per axis (default 200)")
    args = ap.parse_args()

    with uproot.open(args.rootfile) as f:
        a = f["tr"].arrays(FIELDS, library="np")

    m = np.asarray(a["branch_id"]) == -1
    n = int(m.sum())
    if n == 0:
        sys.exit("no direct-decay events (branch_id == -1) found.\n"
                 "  Did you run validation_675_direct3alpha.mac?")

    def col(k):
        v = np.asarray(a[k], dtype=float)[m]
        return v[np.isfinite(v)]

    # ============ figure 1: pooled LAB energy spectrum =================
    e = np.concatenate([col("e_alpha1"), col("e_alpha2"),
                        col("e_alpha3")])
    e = e[e > 0]
    emax = e.max() * 1.05
    ebins = np.linspace(0, emax, 120)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.hist(e, bins=ebins, histtype="step", lw=2.2, color=C_LAB)
    ax.set_xlabel(r"$\alpha$ energy, LAB (MeV)")
    ax.set_ylabel("counts")
    ax.set_xlim(0, emax)
    ax.set_title(r"direct $3\alpha$ - pooled LAB spectrum")
    plt.tight_layout(pad=0.8)
    plt.savefig(f"{args.stem}_spectrum.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {args.stem}_spectrum.png")

    # ============ figure 2: angular, c.m. vs lab =======================
    abins = np.linspace(-1, 1, 41)
    actr = 0.5 * (abins[:-1] + abins[1:])

    cos_cm = clean_cos(np.asarray(a["cos_theta_primary_cm"])[m])
    cos_lab = clean_cos(np.cos(np.asarray(a["theta_lab_alpha1"])[m]))
    h_cm, _ = np.histogram(cos_cm, bins=abins, density=True)
    h_lab, _ = np.histogram(cos_lab, bins=abins, density=True)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.step(actr, h_cm, where="mid", lw=2.4, color=C_CM,
            label=r"c.m. ($3\alpha$ frame)")
    ax.step(actr, h_lab, where="mid", lw=2.4, color=C_LAB, label="lab")
    ax.set_xlabel(r"$\cos\theta$  ($\alpha_1$ vs beam)")
    ax.set_ylabel("normalised yield")
    ax.set_xlim(-1, 1)
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, loc="upper left")
    ax.set_title("angular vs beam: c.m. vs lab")
    plt.tight_layout(pad=0.8)
    plt.savefig(f"{args.stem}_angular.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {args.stem}_angular.png")

    # ============ figure 3: pairwise opening-angle correlation =========
    opening = np.degrees(np.concatenate([
        col("opening_angle_alpha12_cm"),
        col("opening_angle_alpha13_cm"),
        col("opening_angle_alpha23_cm")]))
    obins = np.linspace(0, 180, 61)
    octr = 0.5 * (obins[:-1] + obins[1:])
    h_op, _ = np.histogram(opening, bins=obins, density=True)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.step(octr, h_op, where="mid", lw=2.4, color=C_CM)
    ax.set_xlabel(r"pairwise opening angle $\theta_{ij}$ (deg, c.m.)")
    ax.set_ylabel("normalised yield")
    ax.set_xlim(0, 180)
    ax.set_ylim(bottom=0)
    ax.set_title(r"direct $3\alpha$ opening-angle correlation")
    plt.tight_layout(pad=0.8)
    plt.savefig(f"{args.stem}_opening.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {args.stem}_opening.png")

    # ============ figure 4: symmetrized Dalitz plot =====================
    E = np.column_stack([col("e_3alpha_cm_alpha1"),
                         col("e_3alpha_cm_alpha2"),
                         col("e_3alpha_cm_alpha3")])
    x, y = dalitz_xy_symmetrized(E)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    edges = np.linspace(-LIM, LIM, args.nbins + 1)
    hist2d, _, _ = np.histogram2d(x, y, bins=[edges, edges])
    hist2d = np.ma.masked_equal(hist2d, 0)
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
    ax.set_title(r"Dalitz plot (phase space)")
    plt.tight_layout(pad=0.8)
    plt.savefig(f"{args.stem}_dalitz.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {args.stem}_dalitz.png")

    print(f"  direct-decay events: {n}")


if __name__ == "__main__":
    main()
