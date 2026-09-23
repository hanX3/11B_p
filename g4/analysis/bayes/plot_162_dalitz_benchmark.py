#!/usr/bin/env python3
"""Dalitz plot from a pure single-channel 162-keV benchmark ROOT file.

Each benchmark file is a single reaction channel (one of the three below).
The channel is read from the reaction tree's h11b_reaction_channel label, so
the title and default output name are set automatically; no need to pass the
channel by hand.

    h11b_reaction_channel == 1  ->  Seq162BeGround   (8Be g.s. 0+)
    h11b_reaction_channel == 2  ->  Seq162Be2Plus    (8Be* 2+)
    h11b_reaction_channel == 3  ->  Direct162        (direct 3-alpha)

Dalitz coordinates (Kuhlwein 2022), with E_i the alpha kinetic energies in
the 3-alpha c.m. frame (tree branches e_3alpha_cm_alpha1/2/3):
    x = sqrt(3) (E2 - E3) / (E1 + E2 + E3)
    y = (2 E1 - E2 - E3) / (E1 + E2 + E3)
Every event is filled with all six permutations (identical bosons).

Dual use:
    terminal :  python3 plot_162_dalitz_benchmark.py benchmark_162seq_Be2plus.root
                python3 plot_162_dalitz_benchmark.py f1.root f2.root f3.root
                python3 plot_162_dalitz_benchmark.py f.root --out my_name.png
                (saves one PNG per input file)
    jupyter  :  from plot_162_dalitz_benchmark import plot_162_dalitz_benchmark
                fig = plot_162_dalitz_benchmark("benchmark_162seq_Be2plus.root")
                fig = plot_162_dalitz_benchmark("f.root", save="my_name.png")
                (figure displays inline automatically under %matplotlib inline)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import uproot

TREE = "reaction"

FIELDS = [
    "e_3alpha_cm_alpha1", "e_3alpha_cm_alpha2", "e_3alpha_cm_alpha3",
    "h11b_reaction_channel",
]

# h11b_reaction_channel -> (default output stem, panel title)
CHANNEL_INFO = {
    1: ("dalitz_162seq_Be_gs",   r"162 seq. $^{8}$Be g.s. ($0^{+}$)"),
    2: ("dalitz_162seq_Be2plus", r"162 seq. $^{8}$Be$^{*}$ ($2^{+}$)"),
    3: ("dalitz_162direct",      r"162 direct $3\alpha$"),
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


def plot_162_dalitz_benchmark(rootfile, save=False, nbins=180, title=None):
    """Dalitz plot for a pure single-channel 162-keV benchmark file.

    Parameters
    ----------
    rootfile : str or pathlib.Path
        Pure single-channel benchmark ROOT file (reaction tree).
    save : bool or str, optional
        False (default) -> do not write a file.
        True            -> write dalitz_<channel>.png (auto-named).
        str             -> write to that path (".png" appended if missing).
    nbins : int, optional
        2D bins per axis (default 180).
    title : str, optional
        Override the auto-detected panel title.

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    with uproot.open(rootfile) as f:
        a = f[TREE].arrays(FIELDS, library="np")

    channel = np.asarray(a["h11b_reaction_channel"], dtype=np.int64)
    present = np.unique(channel)

    known = [c for c in present if c in CHANNEL_INFO]
    if not known:
        raise ValueError(
            f"no 162 channel (1/2/3) found; h11b_reaction_channel={present}")
    if len(known) > 1:
        raise ValueError(
            "this file is not a pure single-channel benchmark "
            f"(channels present: {known}). Use plot_162_dalitz_euler.py for a "
            "mixed file, or pass a pure benchmark file here.")

    chan = known[0]
    default_stem, default_title = CHANNEL_INFO[chan]
    panel_title = default_title if title is None else title

    m = channel == chan
    E = np.column_stack([np.asarray(a["e_3alpha_cm_alpha1"], dtype=float)[m],
                         np.asarray(a["e_3alpha_cm_alpha2"], dtype=float)[m],
                         np.asarray(a["e_3alpha_cm_alpha3"], dtype=float)[m]])
    x, y = dalitz_xy_symmetrized(E)
    if x.size == 0:
        raise ValueError(f"no valid Dalitz entries for channel {chan}")

    limits = (-LIM, LIM)
    hist2d, x_edges, y_edges = np.histogram2d(
        x, y, bins=nbins, range=[limits, limits])
    positive = hist2d[hist2d > 0]

    fig, ax = plt.subplots(figsize=(6.2, 5.5))
    im = ax.pcolormesh(
        x_edges, y_edges, hist2d.T, shading="auto", cmap="viridis",
        norm=LogNorm(vmin=1.0, vmax=float(np.max(positive))))
    fig.colorbar(im, ax=ax, pad=0.02, label="Counts")

    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(t), np.sin(t), color="#CCCCCC", lw=1.0, ls="--")

    ax.set_xlabel(r"$x=\sqrt{3}\,(E_2-E_3)/\sum E_i$")
    ax.set_ylabel(r"$y=(2E_1-E_2-E_3)/\sum E_i$")
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(panel_title)
    fig.tight_layout()

    print(f"channel {chan} ({default_stem}): {int(m.sum())} events")

    if save:
        path = f"{default_stem}.png" if save is True else str(save)
        if not path.lower().endswith(".png"):
            path += ".png"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        print(f"wrote {path}")

    return fig


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(
        description="Dalitz plot(s) from pure single-channel 162-keV "
                    "benchmark ROOT file(s).")
    ap.add_argument("rootfiles", nargs="+",
                    help="one or more pure benchmark ROOT files")
    ap.add_argument("--out", default=None,
                    help="output png name (only valid with a single input; "
                         "otherwise names are auto-detected per channel)")
    ap.add_argument("--nbins", type=int, default=180)
    args = ap.parse_args()

    if args.out is not None and len(args.rootfiles) > 1:
        ap.error("--out cannot be used with multiple input files")

    for rf in args.rootfiles:
        out = True if args.out is None else args.out
        plot_162_dalitz_benchmark(rf, save=out, nbins=args.nbins)
