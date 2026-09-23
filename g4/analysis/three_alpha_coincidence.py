#!/usr/bin/env python3
"""Three-alpha Si-segment coincidence analysis.

Quantifies, for the strip-segmented Si array, how often the three alpha
particles of an 11B(p,3alpha) event are recorded on three *distinct* Si
segments (a triple coincidence) -- the observable that the DSSD upgrade
enables and that downstream Dalitz / angular-correlation studies rely on.

Usage:
    python3 analysis/three_alpha_coincidence.py 'data/event_<stamp>_t*.root'

The argument is a glob for the per-thread event ROOT files of a single run.
Each Si hit row is one fired (module, segment); an alpha hit is identified by
pdg == 1000020040.  A segment is counted once per event even if several alphas
share it.  Reported quantities:

  * events with >=1 / >=2 / >=3 distinct alpha Si segments
  * the triple-coincidence fraction (>=3 distinct segments) -- the baseline
  * per-ring hit accounting (1=barrel, 2=backward annular, 3=forward annular)
"""
import sys
import glob

import numpy as np
import uproot

ALPHA_PDG = 1000020040
RING_NAMES = {1: "barrel", 2: "backward_annular", 3: "forward_annular"}


def load(pattern):
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit(f"no files match: {pattern}")
    branches = ["event", "detector_type", "ring_id", "copy_no", "pdg"]
    cols = {b: [] for b in branches}
    for fn in files:
        arr = uproot.open(fn)["tr"].arrays(branches, library="np")
        for b in branches:
            cols[b].append(arr[b])
    for b in branches:
        cols[b] = np.concatenate(cols[b])
    return files, cols


def main():
    pattern = sys.argv[1] if len(sys.argv) > 1 else "data/event_*_t*.root"
    files, d = load(pattern)

    si = d["detector_type"] == 1
    alpha = si & (d["pdg"] == ALPHA_PDG)

    n_events = int(len(np.unique(d["event"]))) if d["event"].size else 0
    n_events_with_si_alpha = int(len(np.unique(d["event"][alpha]))) if alpha.any() else 0

    # distinct alpha Si segments per event: unique copy_no per event
    ev = d["event"][alpha]
    cn = d["copy_no"][alpha].astype(np.int64)
    pair = ev.astype(np.int64) * (10 ** 12) + cn  # unique (event, copy_no) key
    uniq_pairs = np.unique(pair)
    ev_of_pair = uniq_pairs // (10 ** 12)
    seg_counts = np.bincount(np.unique(ev_of_pair, return_inverse=True)[1]) if ev_of_pair.size else np.array([])

    n_ge1 = int((seg_counts >= 1).sum())
    n_ge2 = int((seg_counts >= 2).sum())
    n_ge3 = int((seg_counts >= 3).sum())

    def frac(n):
        return n / n_events if n_events else 0.0

    print("=" * 64)
    print("Three-alpha Si-segment coincidence")
    print("=" * 64)
    print(f"input files            : {len(files)}")
    print(f"primary events         : {n_events}")
    print(f"total Si hit rows       : {int(si.sum())}")
    print(f"alpha Si hit rows       : {int(alpha.sum())}")
    print(f"events w/ >=1 alpha seg : {n_ge1}  ({frac(n_ge1):.4f})")
    print(f"events w/ >=2 alpha seg : {n_ge2}  ({frac(n_ge2):.4f})")
    print(f"events w/ >=3 alpha seg : {n_ge3}  ({frac(n_ge3):.4f})   <-- triple-coincidence baseline")
    print("-" * 64)
    print("per-ring alpha hit rows:")
    for ring, name in RING_NAMES.items():
        m = alpha & (d["ring_id"] == ring)
        print(f"  ring {ring} ({name:17s}): {int(m.sum())}")
    print("=" * 64)
    print(
        f"BASELINE triple-alpha Si-segment coincidence fraction = {frac(n_ge3):.4f} "
        f"({n_ge3}/{n_events})"
    )


if __name__ == "__main__":
    main()
