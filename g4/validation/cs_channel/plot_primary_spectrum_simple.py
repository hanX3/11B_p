#!/usr/bin/env python3
"""Plot primary-alpha energy spectra for alpha0 and alpha1 in one panel.
Four distinct colors; lab = solid, c.m. = dashed. Unfilled step histograms.

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

# four distinct colors
C_A1_LAB = "#C77D2E"   # amber
C_A1_CM  = "#7A4E9E"   # purple
C_A0_LAB = "#0E7C6B"   # teal
C_A0_CM  = "#B23A48"   # red

emax = max(a["e_alpha1"].max(), a["e_alpha1_cm"].max()) * 1.05
bins = np.linspace(0, emax, 150)

plt.figure(figsize=(5.5, 4.5))

series = [
    (a["e_alpha1"][bid == 1],    C_A1_LAB, "-",  r"$\alpha_1$ lab"),
    (a["e_alpha1_cm"][bid == 1], C_A1_CM,  "-", r"$\alpha_1$ c.m."),
    (a["e_alpha1"][bid == 0],    C_A0_LAB, "-",  r"$\alpha_0$ lab"),
    (a["e_alpha1_cm"][bid == 0], C_A0_CM,  "-", r"$\alpha_0$ c.m."),
]
handles = []
for data, color, style, label in series:
    plt.hist(data, bins=bins, histtype="step", lw=1.8, ls=style, color=color)
    handles.append(Line2D([0], [0], color=color, lw=1.8, ls=style, label=label))

plt.xlabel(r"primary $\alpha$ energy (MeV)")
plt.ylabel("counts")
# legend as lines (no boxes), placed upper-left away from the peaks
plt.legend(handles=handles, loc="upper left", frameon=False)
plt.tight_layout()
plt.savefig("primary_alpha_spectrum.png", dpi=200)
print("wrote primary_alpha_spectrum.png")
