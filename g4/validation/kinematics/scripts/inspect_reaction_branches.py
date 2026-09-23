#!/usr/bin/env python3
"""Print every branch in the "tr" reaction tree of a merged ROOT file, plus
reaction_channel value counts, so the other validate_*.py scripts (and
anyone reading a report) can see exactly what was available without
guessing from source code alone.

Usage:
    python3 validation/kinematics/scripts/inspect_reaction_branches.py \
        data/reaction_merged_xxx.root
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402


def main() -> int:
    if len(sys.argv) < 2:
        print("Usage: inspect_reaction_branches.py <root_file> [root_file ...]")
        return 2

    for path_str in sys.argv[1:]:
        path = Path(path_str)
        print(f"==== {path} ====")
        tree = c.open_tree(path)
        n_entries = tree.num_entries
        print(f"tree 'tr': {n_entries} entries, {len(tree.keys())} branches")

        for name in sorted(tree.keys()):
            print(f"  {name}")

        if "reaction_channel" in tree.keys():
            values = tree["reaction_channel"].array(library="np")
            unique, counts = np.unique(values, return_counts=True)
            print("reaction_channel value counts:")
            for value, count in zip(unique, counts):
                print(f"  reaction_channel={value}: {count} ({count / n_entries:.4%})" if n_entries else f"  reaction_channel={value}: {count}")
        else:
            print("MISSING BRANCH: reaction_channel")

        for extra in ("resonance_id", "gamma_resonance", "gamma_branch", "branch_id", "h11b675_alpha_decay_model"):
            if extra in tree.keys():
                values = tree[extra].array(library="np")
                unique, counts = np.unique(values, return_counts=True)
                print(f"{extra} value counts: " + ", ".join(f"{u}:{n}" for u, n in zip(unique, counts)))
            else:
                print(f"MISSING BRANCH: {extra}")
        print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
