#!/usr/bin/env python3
"""675-keV gamma-channel validation figures, 11B(p,gamma)12C (16.62, 2-).

Reads the output of validation_675_gamma.mac and produces two single
4:3 panels (2-row x 3-column grid on a 16:8 slide):

  <stem>_spectrum.png : LAB gamma energy spectrum split by branch.
      gamma_branch ids (Constants.hh):
        3 = to g.s. (0+)      ~16.6 MeV   (M2, 15.7 rel.)
        4 = to 4.44 (2+)      ~12.2 MeV   (dominant E1; + 4.44 cascade)
        5 = to 7.65 (0+)      upper limit
        6 = to 12.71 (1+)     ~3.9 MeV
        7 = to 15.11 (1+)     ~1.5 MeV
      Relative intensities in the generator: Kelley 2017 (Nucl. Phys. A
      968, 71), Table 12.14 footnote x; original measurement Zijderhand
      1990 (Nucl. Instrum. Methods A 286, 490) -- single-angle relative
      intensities at 55 deg (near the P2 magic angle), used here as
      quasi angle-integrated branching ratios.
      Both prompt gammas of an event are filled (gamma1 + gamma2).
  <stem>_angular.png : primary-gamma cos(theta) vs beam (c.m.).
      Expect FLAT: s-wave-formed 2- has equal-m population -- the same
      rotational-invariance theorem as the alpha channel. No angular
      data exists for the 16.62 gamma decay (single-angle measurement
      only); isotropy is the theory-motivated default (A1 = A2 = 0).

Usage:
    python3 plot_675_gamma.py FILE.root
    python3 plot_675_gamma.py FILE.root --stem g675
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
    "legend.fontsize": 10,
    "figure.dpi": 100,
})

FIGSIZE = (16 / 3, 8 / 2)

FIELDS = [
    "gamma_branch", "n_prompt_gammas",
    "gamma1_energy", "gamma2_energy",
    "cos_theta_gamma_cm",
]

# branch id -> (label, color)
BRANCHES = {
    4: (r"to 4.44 ($2^+$), E1", "#990000"),
    3: (r"to g.s. ($0^+$), M2", "#0E7C6B"),
    6: (r"to 12.71 ($1^+$)", "#C77D2E"),
    7: (r"to 15.11 ($1^+$)", "#2E5FA3"),
    5: (r"to 7.65 ($0^+$)", "#888888"),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rootfile")
    ap.add_argument("--stem", default="gamma675")
    args = ap.parse_args()

    with uproot.open(args.rootfile) as f:
        a = f["tr"].arrays(FIELDS, library="np")

    gb = np.asarray(a["gamma_branch"])
    m675 = (gb >= 3) & (gb <= 7)
    if m675.sum() == 0:
        sys.exit("no 675 gamma events (gamma_branch 3-7) found.\n"
                 "  Did you run validation_675_gamma.mac?")

    e1 = np.asarray(a["gamma1_energy"], dtype=float)
    e2 = np.asarray(a["gamma2_energy"], dtype=float)

    # ---------------- figure 1: spectrum by branch ----------------
    emax = 17.5
    ebins = np.linspace(0, emax, 875)   # 20 keV bins

    fig, ax = plt.subplots(figsize=FIGSIZE)
    from matplotlib.lines import Line2D
    handles = []
    for bid, (label, color) in BRANCHES.items():
        mb = gb == bid
        if mb.sum() == 0:
            continue
        e = np.concatenate([e1[mb], e2[mb]])
        e = e[np.isfinite(e) & (e > 0)]
        ax.hist(e, bins=ebins, histtype="step", lw=1.6, color=color)
        handles.append(Line2D([], [], color=color, lw=1.6,
                              label=f"{label}  ({int(mb.sum())})"))
    ax.set_xlabel(r"$\gamma$ energy, LAB (MeV)")
    ax.set_ylabel("counts")
    ax.set_xlim(0, emax)
    ax.set_yscale("log")
    ymax = ax.get_ylim()[1]
    ax.set_ylim(top=ymax * 20)   # headroom for the legend
    ax.legend(handles=handles, loc="upper right", frameon=False)
    ax.set_title(r"$^{11}$B(p,$\gamma$)$^{12}$C - 675 keV")
    plt.tight_layout(pad=0.8)
    plt.savefig(f"{args.stem}_spectrum.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {args.stem}_spectrum.png")

    # ---------------- figure 2: primary-gamma angular ----------------
    c = np.asarray(a["cos_theta_gamma_cm"], dtype=float)[m675]
    c = c[np.isfinite(c) & (c >= -1.0) & (c <= 1.0)]
    abins = np.linspace(-1, 1, 41)
    actr = 0.5 * (abins[:-1] + abins[1:])
    h, _ = np.histogram(c, bins=abins, density=True)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.step(actr, h, where="mid", lw=2.4, color="#990000")
    ax.set_xlabel(r"$\cos\theta_{\gamma}$  (c.m., vs beam)")
    ax.set_ylabel("normalised yield")
    ax.set_xlim(-1, 1)
    ax.set_ylim(0, 1)
    ax.set_title(r"primary $\gamma$ angular (expect isotropic)")
    plt.tight_layout(pad=0.8)
    plt.savefig(f"{args.stem}_angular.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {args.stem}_angular.png")

    # ---------------- console: branching check ----------------
    print("  branch counts (relative to 'to 4.44' = 100):")
    n444 = max(int((gb == 4).sum()), 1)
    for bid, (label, _) in BRANCHES.items():
        nb = int((gb == bid).sum())
        print(f"    id {bid:d} {label:28s}: {nb:>8d}  "
              f"rel = {100.0 * nb / n444:6.2f}")


if __name__ == "__main__":
    main()
