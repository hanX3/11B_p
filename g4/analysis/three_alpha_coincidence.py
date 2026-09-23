#!/usr/bin/env python3
"""Three-alpha DSSD strip-readout multiplicity analysis.

Usage:
    python3 analysis/three_alpha_coincidence.py 'data/event_<stamp>_t*.root'

The event tree stores the two DSSD faces as independent electronic signals:

  readout_side = 0: front side
    W1 -> local-x / azimuthal strip
    S3 -> angular-sector strip

  readout_side = 1: back side
    W1 -> local-y / polar strip
    S3 -> radial-ring strip

A physical deposit contributes the same charge magnitude to one front strip and
one back strip.  Deposits sharing a strip are summed before output.  Therefore
this script reports strip multiplicities; it does not assume a unique front/back
pixel pairing.  The condition >=3 front and >=3 back alpha strips is a necessary
readout-multiplicity condition for resolving three alpha hits, not a complete
pixel-reconstruction efficiency.
"""
import glob
import sys

import numpy as np
import uproot

ALPHA_PDG = 1000020040
FRONT = 0
BACK = 1
SUBARRAY_NAMES = {
    1: "drum",
    2: "backward_s3",
    3: "forward_s3",
    4: "forward_cap",
    5: "backward_cap",
}


def load(pattern):
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit(f"no files match: {pattern}")

    branches = [
        "event",
        "detector_type",
        "subarray_id",
        "module_id",
        "readout_side",
        "strip_id",
        "pdg",
    ]
    cols = {branch: [] for branch in branches}
    for filename in files:
        arrays = uproot.open(filename)["tr"].arrays(branches, library="np")
        for branch in branches:
            cols[branch].append(arrays[branch])

    for branch in branches:
        cols[branch] = np.concatenate(cols[branch]) if cols[branch] else np.array([])
    return files, cols


def per_event_unique_strip_count(data, mask):
    """Count unique (subarray,module,strip) channels for each selected event."""
    if not mask.any():
        return {}

    rows = np.rec.fromarrays(
        [
            data["event"][mask].astype(np.int64),
            data["subarray_id"][mask].astype(np.int32),
            data["module_id"][mask].astype(np.int32),
            data["strip_id"][mask].astype(np.int32),
        ],
        names="event,subarray,module,strip",
    )
    unique_rows = np.unique(rows)
    events, counts = np.unique(unique_rows.event, return_counts=True)
    return dict(zip(events.tolist(), counts.tolist()))


def count_at_least(counts, threshold):
    return sum(value >= threshold for value in counts.values())


def main():
    pattern = sys.argv[1] if len(sys.argv) > 1 else "data/event_*_t*.root"
    files, data = load(pattern)

    si = data["detector_type"] == 1
    alpha = si & (data["pdg"] == ALPHA_PDG)
    front_alpha = alpha & (data["readout_side"] == FRONT)
    back_alpha = alpha & (data["readout_side"] == BACK)

    front_counts = per_event_unique_strip_count(data, front_alpha)
    back_counts = per_event_unique_strip_count(data, back_alpha)

    event_ids = set(front_counts) | set(back_counts)
    n_ge1_both = sum(front_counts.get(event, 0) >= 1 and back_counts.get(event, 0) >= 1 for event in event_ids)
    n_ge2_both = sum(front_counts.get(event, 0) >= 2 and back_counts.get(event, 0) >= 2 for event in event_ids)
    n_ge3_both = sum(front_counts.get(event, 0) >= 3 and back_counts.get(event, 0) >= 3 for event in event_ids)

    # This remains the number of events represented in the hit tree, not the
    # total generated-event denominator.  Use the reaction tree for an absolute
    # efficiency denominator.
    n_saved_events = int(len(np.unique(data["event"]))) if data["event"].size else 0

    def fraction(value):
        return value / n_saved_events if n_saved_events else 0.0

    print("=" * 72)
    print("Three-alpha DSSD independent-strip readout")
    print("=" * 72)
    print(f"input files                         : {len(files)}")
    print(f"events represented in hit tree      : {n_saved_events}")
    print(f"Si alpha front-strip rows            : {int(front_alpha.sum())}")
    print(f"Si alpha back-strip rows             : {int(back_alpha.sum())}")
    print(f"events with >=1 strip on both sides  : {n_ge1_both}  ({fraction(n_ge1_both):.4f})")
    print(f"events with >=2 strips on both sides : {n_ge2_both}  ({fraction(n_ge2_both):.4f})")
    print(f"events with >=3 strips on both sides : {n_ge3_both}  ({fraction(n_ge3_both):.4f})")
    print("-" * 72)
    print("alpha strip rows by subarray and side:")
    for subarray_id, name in SUBARRAY_NAMES.items():
        front_rows = int((front_alpha & (data["subarray_id"] == subarray_id)).sum())
        back_rows = int((back_alpha & (data["subarray_id"] == subarray_id)).sum())
        print(f"  {subarray_id} ({name:13s})  front={front_rows:8d}  back={back_rows:8d}")
    print("=" * 72)
    print(
        ">=3 strips on both sides is a readout multiplicity condition. "
        "Front/back pairing and ghost-pixel rejection require a separate reconstruction step."
    )


if __name__ == "__main__":
    main()
