#!/usr/bin/env python3
"""Validate the alpha0 primary angular distribution: off vs on.

Two figures are produced for the alpha0 channel (branch_id == 0):
  angular_alpha0_cm.png  : cos(theta) in the 12C c.m. frame, with the analytic
                           W(theta)=1+a2 P2 overlaid. This is the actual
                           validation (sampling frame == input frame).
  angular_alpha0_lab.png : cos(theta) in the lab frame, off vs on. Here the
                           distribution is boosted, so it does NOT equal the
                           input W(theta); shown to illustrate the boost.

alpha1 is not shown (isotropic by construction).

Usage:
    python3 plot_angular_alpha0.py OFF.root ON.root
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import uproot

# numpy 1.x/2.x compatibility (trapezoid was added in 2.0, was trapz before)
_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

if len(sys.argv) != 3:
    sys.exit("usage: python3 plot_angular_alpha0.py OFF.root ON.root")

# alpha0 Becker coefficients used in the ON run (for the reference curve)
A1, A2 = 0, 0.6

def load(path):
    with uproot.open(path) as f:
        a = f["tr"].arrays(
            ["cos_theta_primary_cm", "theta_lab_alpha1", "branch_id"],
            library="np")
    m = a["branch_id"] == 0
    cos_cm = a["cos_theta_primary_cm"][m]
    cos_lab = np.cos(a["theta_lab_alpha1"][m])
    def clean(c):
        return c[np.isfinite(c) & (c >= -1) & (c <= 1)]
    return clean(cos_cm), clean(cos_lab)

off_cm, off_lab = load(sys.argv[1])
on_cm,  on_lab  = load(sys.argv[2])

bins = np.linspace(-1, 1, 41)
ctr = 0.5 * (bins[:-1] + bins[1:])

def hist(c):
    h, _ = np.histogram(c, bins=bins, density=True)
    return h

# ---------- c.m. frame (validation, with analytic curve) ----------
x = np.linspace(-1, 1, 400)
W = 1 + A1 * x + A2 * 0.5 * (3 * x**2 - 1)
W = W / (_trapz(W, x) / 2.0) * 0.5

plt.figure(figsize=(6.0, 3.2))
plt.step(ctr, hist(off_cm), where="mid", lw=1.8, color="#C77D2E",
         label=f"off (isotropic), N={off_cm.size}")
plt.step(ctr, hist(on_cm), where="mid", lw=1.8, color="#0E7C6B",
         label=f"on ($a_1$={A1}, $a_2$={A2}), N={on_cm.size}")
plt.plot(x, W, "--", lw=1.4, color="#1A2E35",
         label=r"$1+a_1P_1+a_2P_2$")
plt.axhline(0.5, color="#999999", lw=0.8, ls=":")
plt.xlabel(r"$\cos\theta_{\mathrm{c.m.}}$  (primary $\alpha_0$)")
plt.ylabel("normalised yield")
plt.xlim(-1, 1); plt.ylim(bottom=0)
plt.legend(frameon=False, fontsize=9)
plt.tight_layout()
plt.savefig("angular_alpha0_cm.png", dpi=200)
print("wrote angular_alpha0_cm.png")

# ---------- lab frame (boosted; no analytic curve) ----------
plt.figure(figsize=(6.0, 3.2))
plt.step(ctr, hist(off_lab), where="mid", lw=1.8, color="#C77D2E",
         label=f"off (isotropic), N={off_lab.size}")
plt.step(ctr, hist(on_lab), where="mid", lw=1.8, color="#0E7C6B",
         label=f"on ($a_1$={A1}, $a_2$={A2}), N={on_lab.size}")
plt.axhline(0.5, color="#999999", lw=0.8, ls=":")
plt.xlabel(r"$\cos\theta_{\mathrm{lab}}$  (primary $\alpha_0$)")
plt.ylabel("normalised yield")
plt.xlim(-1, 1); plt.ylim(bottom=0)
plt.legend(frameon=False, fontsize=9)
plt.tight_layout()
plt.savefig("angular_alpha0_lab.png", dpi=200)
print("wrote angular_alpha0_lab.png")
