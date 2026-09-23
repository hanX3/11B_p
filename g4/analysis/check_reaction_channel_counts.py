#!/usr/bin/env python3
import sys
import numpy as np
import uproot

def arr(t, name):
    return t[name].array(library="np") if name in t.keys() else None

def counts(name, a):
    if a is None:
        print(f"{name}: branch not found")
        return
    v, c = np.unique(a, return_counts=True)
    print(f"{name}:")
    for x, n in zip(v, c):
        print(f"  {x}: {n}")
    print(f"  total: {len(a)}")

if len(sys.argv) != 2:
    print("Usage: python3 analysis/check_reaction_channel_counts.py reaction_merged.root")
    sys.exit(2)

f = uproot.open(sys.argv[1])
if "tr" in f:
    t = f["tr"]
elif "reaction" in f:
    t = f["reaction"]
else:
    print("ERROR: no reaction tree found. Expected 'tr' or 'reaction'.")
    print("Available keys:", f.keys())
    sys.exit(1)

ch = arr(t, "reaction_channel")
gres = arr(t, "gamma_resonance")
gbr = arr(t, "gamma_branch")

counts("reaction_channel", ch)

if ch is not None:
    m = ch == 3
    print(f"gamma events: {np.count_nonzero(m)}")
    if gres is not None:
        counts("gamma_resonance for gamma events", gres[m])
    if gbr is not None:
        counts("gamma_branch for gamma events", gbr[m])

for name in [
    "sigma_background_b",
    "sigma_background_sampling_b",
    "channel_probability_background",
    "background_bias_factor",
    "cross_section_bias_factor",
    "gamma_bias_factor",
    "sigma_gamma_total_b",
    "channel_probability_gamma",
]:
    a = arr(t, name)
    if a is None or len(a) == 0:
        print(f"{name}: branch not found or empty")
    else:
        finite = a[np.isfinite(a)]
        if len(finite):
            print(f"{name}: first={finite[0]:.8e}, min={finite.min():.8e}, max={finite.max():.8e}")

from_table = arr(t, "gamma_branch_from_relative_table")
if from_table is not None and ch is not None:
    counts("gamma_branch_from_relative_table for gamma events", from_table[ch == 3])
