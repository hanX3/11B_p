#!/usr/bin/env python3
"""Single-panel Dalitz plot, one decay mechanism per figure.

Each run of this script reads ONE reaction ROOT file and produces ONE
png, sized so that six of them tile a 16:8 PowerPoint slide as a
2-row x 3-column grid (each panel = 16/3 : 8/2 = 4:3 aspect).

Mechanism is chosen with --mode and selected via branch_id:
    alpha0  : branch_id ==  0   (12C* -> a0 + 8Be g.s. 0+)
    alpha1  : branch_id ==  1   (12C* -> a1 + 8Be* 2+)
    direct  : branch_id == -1   (democratic phase space)

Dalitz coordinates follow Kuhlwein 2022 (PLB 825, 136857):
    x = sqrt(3) (E2 - E3) / (E1 + E2 + E3)
    y = (2 E1 - E2 - E3) / (E1 + E2 + E3)
with E_i the alpha kinetic energies in the 3-alpha c.m. frame
(tree branches e_3alpha_cm_alpha1/2/3). Each event is filled with all
six permutations of (E1, E2, E3) -- identical bosons, exactly as an
experiment that cannot tag the primary alpha would populate the plot.
The kinematic boundary (unit circle) is drawn as a dashed reference.

Usage:
    python3 plot_dalitz_panel.py alpha0.root --mode alpha0
    python3 plot_dalitz_panel.py alpha1.root --mode alpha1 --title "my title"
    python3 plot_dalitz_panel.py direct.root --mode direct --out dalitz_direct
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import uproot

# ---- styling: panel will sit at 1/3 slide width, keep fonts readable ----
plt.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 15,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "axes.linewidth": 1.1,
    "figure.dpi": 100,
})

FIELDS = [
    "e_3alpha_cm_alpha1", "e_3alpha_cm_alpha2", "e_3alpha_cm_alpha3",
    "branch_id",
]

MODES = {
    "alpha0": (0,  r"$\alpha_0$: $^8$Be g.s. ($0^+$)"),
    "alpha1": (1,  r"$\alpha_1$: $^8$Be$^*$ ($2^+$)"),
    "direct": (-1, r"direct $3\alpha$ (phase space)"),
}

LIM = 1.15   # axis range: unit circle + margin


def dalitz_xy_symmetrized(E):
    """E: (N, 3) c.m. kinetic energies -> x, y with all 6 permutations."""
    perms = [(0, 1, 2), (0, 2, 1), (1, 0, 2),
             (1, 2, 0), (2, 0, 1), (2, 1, 0)]
    s = E.sum(axis=1)
    ok = np.isfinite(s) & (s > 0)
    E, s = E[ok], s[ok]
    xs, ys = [], []
    for i, j, k in perms:
        xs.append(np.sqrt(3.0) * (E[:, j] - E[:, k]) / s)
        ys.append((2.0 * E[:, i] - E[:, j] - E[:, k]) / s)
    return np.concatenate(xs), np.concatenate(ys)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rootfile")
    ap.add_argument("--mode", required=True, choices=MODES)
    ap.add_argument("--title", default=None,
                    help="override the default panel title")
    ap.add_argument("--out", default=None,
                    help="output png stem (default: dalitz_<mode>)")
    ap.add_argument("--nbins", type=int, default=200,
                    help="2D bins per axis (default 200; try 120 if the "
                         "alpha0 lines look broken)")
    args = ap.parse_args()

    want_bid, default_title = MODES[args.mode]
    title = args.title if args.title is not None else default_title
    stem  = args.out if args.out is not None else f"dalitz_{args.mode}"

    with uproot.open(args.rootfile) as f:
        a = f["tr"].arrays(FIELDS, library="np")

    bid = np.asarray(a["branch_id"])
    m = bid == want_bid
    n = int(m.sum())
    if n == 0:
        sys.exit(f"no events with branch_id == {want_bid} "
                 f"({args.mode}) in {args.rootfile}")

    E = np.column_stack([np.asarray(a["e_3alpha_cm_alpha1"])[m],
                         np.asarray(a["e_3alpha_cm_alpha2"])[m],
                         np.asarray(a["e_3alpha_cm_alpha3"])[m]])
    x, y = dalitz_xy_symmetrized(E)

    # panel aspect for a 2x3 grid on a 16:8 slide -> (16/3) : (8/2) = 4:3
    fig, axis = plt.subplots(figsize=(16 / 3, 8 / 2))

    edges = np.linspace(-LIM, LIM, args.nbins + 1)
    h, _, _ = np.histogram2d(x, y, bins=[edges, edges])
    h = np.ma.masked_equal(h, 0)          # empty bins -> white
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("white")
    im = axis.pcolormesh(edges, edges, h.T, cmap=cmap,
                         rasterized=True)
    fig.colorbar(im, ax=axis, fraction=0.046, pad=0.03)

    t = np.linspace(0, 2 * np.pi, 400)
    axis.plot(np.cos(t), np.sin(t), color="#CCCCCC", lw=1.0, ls="--")

    axis.set_xlabel(r"$x=\sqrt{3}\,(E_2-E_3)/\sum E_i$")
    axis.set_ylabel(r"$y=(2E_1-E_2-E_3)/\sum E_i$")
    axis.set_xlim(-LIM, LIM)
    axis.set_ylim(-LIM, LIM)
    axis.set_aspect("equal")
    axis.set_title(title)

    plt.tight_layout(pad=0.8)
    plt.savefig(f"{stem}.png", dpi=200, bbox_inches="tight")
    print(f"wrote {stem}.png   ({args.mode}: {n} events)")


if __name__ == "__main__":
    main()
