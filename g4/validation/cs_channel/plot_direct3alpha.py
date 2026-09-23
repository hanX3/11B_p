#!/usr/bin/env python3
"""Direct-decay (democratic phase-space) 3-alpha diagnostics.

Three panels side-by-side, 16:8 layout for a single PowerPoint slide:

    (left)   pooled LAB energy spectrum (all three alphas)
    (middle) alpha1 angular: c.m. (isotropic) vs LAB (boost-tilted) overlay
    (right)  pairwise c.m. opening-angle correlation (12/13/23 pooled)

Direct-decay events are selected via  branch_id == -1.

Middle panel uses alpha1 ONLY for a clean, like-for-like c.m.-vs-lab
comparison, because in the direct-decay path the 3-alpha-c.m. polar angle
(cos_theta_primary_cm) is stored for alpha1 only, while all three alphas
have a lab angle. Using the same particle for both curves isolates the
Lorentz-boost effect: c.m. is flat (isotropic), lab is forward-tilted.

Branches read from tree "tr":
    e_alpha1/2/3               : lab kinetic energies (MeV)
    theta_lab_alpha1           : alpha1 lab polar angle (rad)
    cos_theta_primary_cm       : alpha1 polar angle in 3-alpha c.m. (cos)
    opening_angle_alpha12/13/23_cm : pairwise c.m. opening angles (rad)
    branch_id                  : -1 for direct decay

Usage:
    python3 plot_direct3alpha_3panel.py FILE.root            # single run
    python3 plot_direct3alpha_3panel.py FILE.root  out_stem  # custom png stem
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import uproot

# ---- 16:8 PPT styling ----
plt.rcParams.update({
    "font.size": 15,
    "axes.titlesize": 16,
    "axes.labelsize": 15,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "axes.linewidth": 1.2,
    "legend.fontsize": 13,
    "figure.dpi": 100,
})

if len(sys.argv) not in (2, 3):
    sys.exit("usage: python3 plot_direct3alpha_3panel.py FILE.root [out_stem]")

INFILE   = sys.argv[1]
OUT_STEM = sys.argv[2] if len(sys.argv) == 3 else "direct3alpha_3panel"

FIELDS = [
    "e_alpha1", "e_alpha2", "e_alpha3",
    "theta_lab_alpha1",
    "cos_theta_primary_cm",
    "opening_angle_alpha12_cm", "opening_angle_alpha13_cm",
    "opening_angle_alpha23_cm",
    "branch_id",
]

C_LAB = "#0E7C6B"   # teal   = lab
C_CM  = "#C77D2E"   # orange = c.m. / opening angle


def load(path):
    with uproot.open(path) as f:
        return f["tr"].arrays(FIELDS, library="np")


def pool(arr, keys, mask):
    return np.concatenate([np.asarray(arr[k])[mask] for k in keys])


def finite(x):
    x = np.asarray(x, dtype=float)
    return x[np.isfinite(x)]


def clean_cos(c):
    c = np.asarray(c, dtype=float)
    return c[np.isfinite(c) & (c >= -1.0) & (c <= 1.0)]


a = load(INFILE)

bid = np.asarray(a["branch_id"])
m = bid == -1
if m.sum() == 0:
    sys.exit("no direct-decay events (branch_id == -1) found in this file.\n"
             "  Did you run validation_162/675_direct3alpha.mac "
             "(SequentialDecayFraction = 0)?")

# --- panel 1: pooled LAB energy spectrum ---
e_lab = finite(pool(a, ["e_alpha1", "e_alpha2", "e_alpha3"], m))
emax  = (e_lab.max() if e_lab.size else 0) * 1.05
ebins = np.linspace(0, emax, 150)

# --- panel 2: alpha1 c.m. vs lab cos(theta) ---
cos_cm  = clean_cos(np.asarray(a["cos_theta_primary_cm"])[m])
cos_lab = clean_cos(np.cos(np.asarray(a["theta_lab_alpha1"])[m]))

# --- panel 3: pooled pairwise c.m. opening angle ---
opening_deg = np.degrees(finite(pool(a, ["opening_angle_alpha12_cm",
                                         "opening_angle_alpha13_cm",
                                         "opening_angle_alpha23_cm"], m)))

# 16:8 aspect, three columns
fig, ax = plt.subplots(1, 3, figsize=(16, 8))

# ---------------- panel 1: energy ----------------
ax[0].hist(e_lab, bins=ebins, histtype="step", lw=2.2, color=C_LAB)
ax[0].set_xlabel(r"$\alpha$ energy (MeV)")
ax[0].set_ylabel("counts")
ax[0].set_xlim(0, emax)
ax[0].set_title(r"direct $3\alpha$ - pooled LAB spectrum")

# ---------------- panel 2: c.m. vs lab angular (alpha1) ----------------
abins = np.linspace(-1, 1, 41)
actr  = 0.5 * (abins[:-1] + abins[1:])
h_cm,  _ = np.histogram(cos_cm,  bins=abins, density=True)
h_lab, _ = np.histogram(cos_lab, bins=abins, density=True)
ax[1].step(actr, h_cm,  where="mid", lw=2.4, color=C_CM,
           label=r"c.m. ($3\alpha$ frame)")
ax[1].step(actr, h_lab, where="mid", lw=2.4, color=C_LAB,
           label="lab")
ax[1].axhline(0.5, color="#999999", lw=1.0, ls=":")
ax[1].set_xlabel(r"$\cos\theta$  ($\alpha_1$)")
ax[1].set_ylabel("normalised yield")
ax[1].set_xlim(-1, 1); ax[1].set_ylim(0, 1)
ax[1].legend(loc="upper left", frameon=False)
ax[1].set_title(r"$\alpha_1$ angular: c.m. vs lab")

# ---------------- panel 3: pairwise opening angle ----------------
obins = np.linspace(0, 180, 61)
octr  = 0.5 * (obins[:-1] + obins[1:])
h_op, _ = np.histogram(opening_deg, bins=obins, density=True)
ax[2].step(octr, h_op, where="mid", lw=2.4, color=C_CM)
ax[2].set_xlabel(r"pairwise opening angle $\theta_{ij}$ (deg, c.m.)")
ax[2].set_ylabel("normalised yield")
ax[2].set_xlim(0, 180)
ax[2].set_ylim(bottom=0)
ax[2].set_title(r"direct $3\alpha$ opening-angle correlation")

plt.tight_layout(pad=1.5)
plt.savefig(f"{OUT_STEM}.png", dpi=200, bbox_inches="tight")
print(f"wrote {OUT_STEM}.png")
print(f"  direct-decay events: {int(m.sum())}  "
      f"(pooled alphas: {e_lab.size})")
