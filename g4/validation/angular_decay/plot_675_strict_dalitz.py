#!/usr/bin/env python3
"""Plot normalized-energy Dalitz coordinates for a 675 strict sample."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import _common as c


BRANCHES = [
    "reaction_channel", "resonance_id",
    "e_3alpha_cm_alpha1", "e_3alpha_cm_alpha2", "e_3alpha_cm_alpha3",
]


def dalitz(data, mask):
    e1 = data["e_3alpha_cm_alpha1"][mask]
    e2 = data["e_3alpha_cm_alpha2"][mask]
    e3 = data["e_3alpha_cm_alpha3"][mask]
    total = e1 + e2 + e3
    x = np.sqrt(3.0) * (e2 - e3) / total
    y = (2.0 * e1 - e2 - e3) / total
    return x[np.isfinite(x) & np.isfinite(y)], y[np.isfinite(x) & np.isfinite(y)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    data, missing, _ = c.load_arrays(args.input, BRANCHES)
    if missing:
        raise SystemExit(f"missing branches: {missing}")
    mask = (data["reaction_channel"] == 1) & (data["resonance_id"] == 675)
    x, y = dalitz(data, mask)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(5.8, 5.2))
    ax.hist2d(x, y, bins=60, range=[[-1.8, 1.8], [-1.8, 1.8]], cmap="viridis")
    ax.set_xlabel("sqrt(3) * (E2 - E3) / sum(E)")
    ax.set_ylabel("(2E1 - E2 - E3) / sum(E)")
    ax.set_title(f"675 strict Dalitz: {args.label} (n={x.size})")
    fig.tight_layout()
    fig.savefig(outdir / f"dalitz_675_strict_{args.label}.png", dpi=140)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

