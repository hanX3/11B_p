#!/usr/bin/env python3
"""Plot primary-alpha energy spectra in one panel.
Two colors distinguish the frame (c.m. vs lab); alpha0/alpha1 are marked
with text annotations on the peaks.

    branch_id == 0  ->  alpha0  (8Be g.s.)
    branch_id == 1  ->  alpha1  (8Be 2+)

Usage:
    python3 plot_primary_spectrum.py FILE.root
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import uproot

if len(sys.argv) != 2:
    sys.exit("usage: python3 plot_primary_spectrum.py FILE.root")

with uproot.open(sys.argv[1]) as f:
    a = f["tr"].arrays(["e_alpha1", "e_alpha1_cm", "branch_id"], library="np")

bid = a["branch_id"]
C_CM  = "#B23A48"   # red   = c.m. frame
C_LAB = "#0E7C6B"   # teal  = lab frame

emax = max(a["e_alpha1"].max(), a["e_alpha1_cm"].max()) * 1.05
bins = np.linspace(0, emax, 150)

plt.figure(figsize=(5.5, 4.5))

# both channels, split only by frame (color)
plt.hist(a["e_alpha1_cm"], bins=bins, histtype="step", lw=1.6, color=C_CM)
plt.hist(a["e_alpha1"],    bins=bins, histtype="step", lw=1.6, color=C_LAB)

# annotate alpha0 / alpha1 by peak position (median of each channel, c.m.)
e_cm = a["e_alpha1_cm"]
a0_x = np.median(e_cm[bid == 0])
a1_x = np.median(e_cm[bid == 1])
ymax = plt.ylim()[1]
plt.annotate(r"$\alpha_1$", xy=(a1_x, ymax * 0.9), ha="center", fontsize=13)
plt.annotate(r"$\alpha_0$", xy=(a0_x+0.3, ymax * 0.5), ha="center", fontsize=13)

# legend: colors = frames only
handles = [Line2D([0], [0], color=C_CM,  lw=1.6, label="c.m."),
           Line2D([0], [0], color=C_LAB, lw=1.6, label="lab")]
plt.legend(handles=handles, loc="upper left", frameon=False)

plt.xlabel(r"primary $\alpha$ energy (MeV)")
plt.ylabel("counts")
plt.tight_layout()
plt.savefig("primary_alpha_spectrum.png", dpi=200)
print("wrote primary_alpha_spectrum.png")
