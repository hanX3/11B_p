#!/usr/bin/env python3
"""Secondary-alpha (alpha2 + alpha3) LAB energy spectra AND angular
distributions, split by primary branch, in a single 2x2 figure.
16:9 layout optimised for PowerPoint slides.

    branch_id == 0  ->  alpha0  (12C* -> a0 + 8Be g.s.),  8Be(0+) -> a + a
    branch_id == 1  ->  alpha1  (12C* -> a1 + 8Be* 2+),   8Be(2+) -> a + a

Panels:
    (top-left)     alpha0 secondary LAB energy spectrum
    (top-right)    alpha1 secondary LAB energy spectrum
    (bottom-left)  alpha0 secondary angular distribution   (expect flat)
    (bottom-right) alpha1 secondary angular distribution

Only the LAB frame is shown for the spectra; both spectrum panels share the
same x-range. Angular panels have y fixed to [0, 1].

Angular variables (per the reaction code):
    alpha1 branch : cos_theta_secondary_correlation  (vs primary-a1 in 12C c.m.)
    alpha0 branch : cos_chi_secondary_8be            (vs 8Be recoil axis)

Give two files (OFF then ON) to overlay the alpha1 angular baseline vs the
Treado A2/A4 anisotropy in the bottom-right panel. Spectra use the first file.

Usage:
    python3 plot_secondary.py FILE.root                # single run
    python3 plot_secondary.py OFF.root ON.root         # baseline vs physics
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import uproot

# ---- 16:9 PPT styling ----
plt.rcParams.update({
    "font.size": 15,
    "axes.titlesize": 16,
    "axes.labelsize": 15,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "axes.linewidth": 1.2,
    "figure.dpi": 100,
})

_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

if len(sys.argv) not in (2, 3):
    sys.exit("usage: python3 plot_secondary.py FILE.root [ON.root]")

FIELDS = [
    "e_alpha2", "e_alpha3",
    "cos_theta_secondary_correlation", "cos_chi_secondary_8be",
    "branch_id",
]

def load(path):
    with uproot.open(path) as f:
        return f["tr"].arrays(FIELDS, library="np")

def clean_cos(c):
    c = np.asarray(c, dtype=float)
    return c[np.isfinite(c) & (c >= -1.0) & (c <= 1.0)]

def both(arr, k2, k3, mask):
    return np.concatenate([np.asarray(arr[k2])[mask],
                           np.asarray(arr[k3])[mask]])

C_LAB = "#0E7C6B"   # teal = lab
C_OFF = "#C77D2E"   # orange = file 1
C_ON  = "#0E7C6B"   # teal   = file 2

a_off = load(sys.argv[1])
a_on  = load(sys.argv[2]) if len(sys.argv) == 3 else None

bid = a_off["branch_id"]
m0 = bid == 0
m1 = bid == 1

# common x-range for both spectrum panels
e0 = both(a_off, "e_alpha2", "e_alpha3", m0); e0 = e0[np.isfinite(e0)]
e1 = both(a_off, "e_alpha2", "e_alpha3", m1); e1 = e1[np.isfinite(e1)]
emax = max(e0.max() if e0.size else 0, e1.max() if e1.size else 0) * 1.05
xbins = np.linspace(0, emax, 150)

# 16:9 aspect ratio
fig, ax = plt.subplots(2, 2, figsize=(16, 8))

# ---------------- LAB energy spectra (shared x-range, no legend) -------------
def energy_panel(axis, elab, title):
    if elab.size == 0:
        axis.text(0.5, 0.5, "no events", ha="center", va="center",
                  transform=axis.transAxes)
        axis.set_title(title); return
    axis.hist(elab, bins=xbins, histtype="step", lw=2.2, color=C_LAB)
    axis.set_xlabel(r"secondary $\alpha$ energy (MeV)")
    axis.set_ylabel("counts")
    axis.set_xlim(0, emax)
    axis.set_title(title)

energy_panel(ax[0, 0], e0,
             r"$\alpha_0$ branch ($^8$Be g.s.) - secondary spectrum")
energy_panel(ax[0, 1], e1,
             r"$\alpha_1$ branch ($^8$Be$^*\,2^+$) - secondary spectrum")

# ---------------- angular distributions (y fixed to 0-1) ----------------
abins = np.linspace(-1, 1, 41)
actr  = 0.5 * (abins[:-1] + abins[1:])

def ang_hist(c):
    h, _ = np.histogram(c, bins=abins, density=True)
    return h

# alpha0: 8Be(0+) isotropic -> cos_chi_secondary_8be (no legend needed)
c0 = clean_cos(a_off["cos_chi_secondary_8be"][m0])
ax[1, 0].step(actr, ang_hist(c0), where="mid", lw=2.4, color=C_LAB)
ax[1, 0].axhline(0.5, color="#999999", lw=1.0, ls=":")
ax[1, 0].set_xlabel(r"$\cos\theta$  (secondary $\alpha$ vs $^8$Be axis, $\alpha_0$)")
ax[1, 0].set_ylabel("normalised yield")
ax[1, 0].set_xlim(-1, 1); ax[1, 0].set_ylim(0, 1)
ax[1, 0].set_title(r"$\alpha_0$ secondary angular ($^8$Be $0^+$: isotropic)")

# alpha1: 8Be(2+) -> cos_theta_secondary_correlation
c1_off = clean_cos(a_off["cos_theta_secondary_correlation"][m1])
ax[1, 1].step(actr, ang_hist(c1_off), where="mid", lw=2.4, color=C_OFF)

if a_on is not None:
    m1_on = a_on["branch_id"] == 1
    c1_on = clean_cos(a_on["cos_theta_secondary_correlation"][m1_on])
    ax[1, 1].step(actr, ang_hist(c1_on), where="mid", lw=2.4, color=C_ON)
    A2, A4 = -0.489, 0.647
    x = np.linspace(-1, 1, 400)
    P2 = 0.5 * (3 * x**2 - 1)
    P4 = (35 * x**4 - 30 * x**2 + 3) / 8
    W = 1 + A2 * P2 + A4 * P4
    W = W / (_trapz(W, x) / 2.0) * 0.5
    ax[1, 1].plot(x, W, "--", lw=2.0, color="#1A2E35")

ax[1, 1].axhline(0.5, color="#999999", lw=1.0, ls=":")
ax[1, 1].set_xlabel(r"$\cos\theta''$  (secondary $\alpha$, $\alpha_1$ branch)")
ax[1, 1].set_ylabel("normalised yield")
ax[1, 1].set_xlim(-1, 1); ax[1, 1].set_ylim(0, 1)
ax[1, 1].set_title(r"$\alpha_1$ secondary angular ($^8$Be $2^+$)")

plt.tight_layout(pad=1.5)
plt.savefig("secondary_alpha_overview.png", dpi=200, bbox_inches="tight")
print("wrote secondary_alpha_overview.png")
print(f"  alpha0 events: {int(m0.sum())},  alpha1 events: {int(m1.sum())}")
