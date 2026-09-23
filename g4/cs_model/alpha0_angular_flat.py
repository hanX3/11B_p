#!/usr/bin/env python3
"""Theoretical primary-alpha0 angular distribution (flat single panel).

    W(theta) = 1 + a1 P1(cos theta) + a2 P2(cos theta)

Flat (short) layout for a slide; legend sits in the empty space directly above
the 90-deg dip so it never overlaps the curve. Overlays the digitized Becker
Fig.10 (148 keV) points.
"""
from __future__ import annotations
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
import numpy as np
import matplotlib.pyplot as plt

# ---- coefficients (code convention) ----
A1 = 0.0
A2 = 0.6
# Becker Fig.10, 148 keV digitized points (theta_deg, relative yield)
BECKER_148 = [(0, 9.5), (30, 7.5), (60, 5.0), (90, 4.2),
              (120, 5.5), (150, 8.0), (180, 9.5)]

TEAL = "#0E7C6B"
INK = "#1A2E35"


def W(theta_rad, a1, a2):
    c = np.cos(theta_rad)
    return 1.0 + a1 * c + a2 * 0.5 * (3.0 * c**2 - 1.0)


def main() -> None:
    th_deg = np.linspace(0, 180, 361)
    w = W(np.radians(th_deg), A1, A2)

    fig, ax = plt.subplots(figsize=(7.2, 3.0))   # wide + short = flat
    ax.plot(th_deg, w, "-", color=TEAL, lw=2.4,
            label=r"$W(\theta)=1+%.2f\,P_1+%.2f\,P_2$" % (A1, A2))

    bt = np.array([p[0] for p in BECKER_148])
    by = np.array([p[1] for p in BECKER_148])
    ax.plot(bt, by * (w.mean() / by.mean()), "o", color=INK, ms=6,
            label="Becker Fig.10, 148 keV (digitized)")

    ax.set_yscale("log")
    ax.set_xlim(0, 180)
    ax.set_xticks([0, 30, 60, 90, 120, 150, 180])
    ax.set_xlabel(r"$\theta_{\mathrm{c.m.}}$ (deg)")
    ax.set_ylabel(r"$W(\theta)$  (rel. units)")
    ax.set_title(r"$^{11}$B(p,$\alpha_0$) primary angular distribution")
    ax.grid(alpha=0.25, which="both")

    # headroom, then place legend right above the dip (empty region)
    ax.set_ylim(top=w.max() * 3.5)
    ax.legend(fontsize=8, loc="center", bbox_to_anchor=(0.5, 0.58),
              framealpha=0.95, handlelength=1.5)

    fig.tight_layout()
    fig.savefig("alpha0_angular_flat.png", dpi=200)
    print("wrote alpha0_angular_flat.png")


if __name__ == "__main__":
    main()
