#!/usr/bin/env python3
"""Plot the theoretical primary-alpha0 angular distribution

    W(theta) = 1 + a1 * P1(cos theta) + a2 * P2(cos theta)

for visual comparison against Becker et al. 1987, Fig. 10 (the E = 148 keV
panel of the 11B(p, alpha0) 8Be reaction).

Two panels are produced:
  (left)  Cartesian, theta on x-axis, log y  -> matches Becker Fig. 10 layout
  (right) polar, beam axis vertical          -> shows the emission pattern

Coefficients are in the *code convention* (theta relative to the proton beam
axis in the 12C c.m. frame). Edit A1, A2 below to explore.

Becker Fig. 10, E = 148 keV, digitized -> a1 ~ 0.0, a2 ~ +0.6 (a2 = 0.6 +/- 0.15;
Becker flags the 150 keV data as uncertain). The panel shows a clear 90-deg
dip with near forward-backward symmetry.
"""

from __future__ import annotations

import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# coefficients to plot (code convention; edit these)
A1 = 0.0
A2 = 0.6
LABEL = f"$W(\\theta)=1{A1:+.2f}\\,P_1{A2:+.2f}\\,P_2$"
# optional: overlay Becker's digitized data points for the 148 keV panel
# (theta_deg, relative yield). Set to None to hide.
BECKER_148 = [(0, 9.5), (30, 7.5), (60, 5.0), (90, 4.2),
              (120, 5.5), (150, 8.0), (180, 9.5)]
# scale factor to place the (arbitrary-unit) Becker points on top of W(theta):
# choose so that the mean matches; adjusted automatically below.
# ----------------------------------------------------------------------

TEAL = "#0E7C6B"
INK = "#1A2E35"
AMBER = "#C77D2E"


def W(theta_rad: np.ndarray, a1: float, a2: float) -> np.ndarray:
    c = np.cos(theta_rad)
    p1 = c
    p2 = 0.5 * (3.0 * c**2 - 1.0)
    return 1.0 + a1 * p1 + a2 * p2


def main() -> None:
    th_deg = np.linspace(0, 180, 361)
    th = np.radians(th_deg)
    w = W(th, A1, A2)

    fig = plt.figure(figsize=(11, 4.6))

    # -------- left: Cartesian (Becker Fig.10 style) --------
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.plot(th_deg, w, "-", color=TEAL, lw=2.4, label=LABEL)

    if BECKER_148 is not None:
        bt = np.array([p[0] for p in BECKER_148])
        by = np.array([p[1] for p in BECKER_148])
        # scale Becker points so their mean matches the curve's mean
        scale = w.mean() / by.mean()
        ax1.plot(bt, by * scale, "o", color=INK, ms=6,
                 label="Becker Fig.10, 148 keV (digitized)")

    ax1.set_yscale("log")
    ax1.set_xlim(0, 180)
    ax1.set_xticks([0, 30, 60, 90, 120, 150, 180])
    ax1.set_xlabel(r"$\theta_{\mathrm{c.m.}}$ (deg)")
    ax1.set_ylabel(r"$W(\theta)$  (rel. units)")
    ax1.set_title(r"$^{11}$B(p,$\alpha_0$) primary angular distribution")
    ax1.grid(alpha=0.25, which="both")
    ax1.legend(fontsize=8, loc="lower center")

    # -------- right: polar (beam axis vertical) --------
    ax2 = fig.add_subplot(1, 2, 2, projection="polar")
    # full 0..360 for a closed lobe; W is symmetric in phi
    tfull = np.linspace(0, 2 * np.pi, 720)
    # map polar angle measured from +y (beam) : theta here = angle from beam axis
    wfull = W(tfull, A1, A2)
    wfull = np.clip(wfull, 0, None)
    # set 0 deg at top (beam direction), clockwise
    ax2.set_theta_zero_location("N")
    ax2.set_theta_direction(-1)
    ax2.plot(tfull, wfull, "-", color=TEAL, lw=2.4)
    ax2.fill(tfull, wfull, color=TEAL, alpha=0.12)
    ax2.set_title("polar view (beam axis = up)", pad=18)
    ax2.set_rlabel_position(135)

    fig.suptitle(f"a1 = {A1:+.2f},  a2 = {A2:+.2f}  (code convention)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = "alpha0_angular_theory.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")

    # numeric summary for quick checks
    print(f"W(0)   = {W(np.radians(0), A1, A2):.3f}")
    print(f"W(90)  = {W(np.radians(90), A1, A2):.3f}")
    print(f"W(180) = {W(np.radians(180), A1, A2):.3f}")
    print(f"W(90)/W(0) = {W(np.radians(90),A1,A2)/W(np.radians(0),A1,A2):.3f}"
          f"   (Becker 148 keV panel reads ~ 0.44)")


if __name__ == "__main__":
    main()
