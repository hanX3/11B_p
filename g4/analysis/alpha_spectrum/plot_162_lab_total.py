#!/usr/bin/env python3
"""Plot the 5/95/0 inclusive alpha LAB spectrum at reaction generation.

Pass either one merged file or the disjoint worker files, never both.
Terminal: python3 plot_162_lab_total.py [FILE.root] [--out OUTPUT_STEM]
Total curve only: python3 plot_162_lab_total.py --total-only
Jupyter:  fig = plot_spectrum()  # or plot_spectrum("FILE.root", save=True)

Defaults are relative to this script, independent of the working directory.
Terminal runs save PNG and PDF next to this script. Importing
does not change the Matplotlib backend. The function returns an open Figure
and only writes files when save=True or an output stem is supplied.
Each alpha contributes one count; the bin counts sum to three times the
number of reactions. No normalization or detector cuts.
"""

import argparse
from pathlib import Path

import numpy as np
import uproot

HERE = Path(__file__).resolve().parent
DEFAULT_ROOT = HERE.parents[1] / "data/spectrum_162_lab_05_95_00.root"
DEFAULT_OUT = HERE / "alpha_lab_total"


def plot_spectrum(rootfiles=None, out=None, *, save=False, show_components=True):
    """Return an editable Figure for terminal or Jupyter use.

    rootfiles: one path, a list of paths, or None for the existing 5/95/0 run.
    out: optional output path stem; supplying it also enables saving.
    save: write PNG/PDF, default False for notebook use.
    show_components: False draws only the black total curve.
    """
    import matplotlib.pyplot as plt

    if rootfiles is None:
        rootfiles = [DEFAULT_ROOT]
    elif isinstance(rootfiles, (str, Path)):
        rootfiles = [rootfiles]
    paths = [Path(p).resolve() for p in rootfiles]
    if not paths:
        raise ValueError("No input files")
    if len(set(paths)) != len(paths):
        raise ValueError("Duplicate input files")
    fields = [
        "e_alpha1", "e_alpha2", "e_alpha3", "projectile_kinetic_lab",
        "resonance_id", "h11b_reaction_channel", "event_weight",
        "configured_162_alpha0_branching_fraction",
        "sequential_decay_fraction_162", "direct_decay_fraction_162",
    ]
    pieces = []
    for path in paths:
        # These are local files. Explicit memory mapping avoids the installed
        # fsspec asynchronous reader, which stalls on basket reads here.
        with uproot.open(path, handler=uproot.source.file.MemmapSource) as f:
            pieces.append(f["reaction"].arrays(fields, library="np"))
    a = {k: np.concatenate([p[k] for p in pieces]) for k in fields}
    n = len(a["resonance_id"])
    if n == 0:
        raise ValueError("No reaction entries")
    channel = a["h11b_reaction_channel"]
    if not (np.all(a["resonance_id"] == 162) and np.all(np.isin(channel, [1, 2]))):
        raise ValueError("Input must contain only 162-keV gs/excited reactions")
    for key, expected in [
        ("configured_162_alpha0_branching_fraction", 0.05),
        ("sequential_decay_fraction_162", 1.0),
        ("direct_decay_fraction_162", 0.0),
    ]:
        if not np.allclose(a[key], expected, rtol=0, atol=1e-10):
            raise ValueError(f"Incorrect run configuration: {key}")

    energies = np.column_stack([a[f"e_alpha{i}"] for i in (1, 2, 3)])
    weights = a["event_weight"]
    proton = a["projectile_kinetic_lab"]
    if not (np.all(np.isfinite(energies)) and np.all(energies > 0)
            and np.all(np.isfinite(weights)) and np.all(weights > 0)
            and np.all(np.isfinite(proton))):
        raise ValueError("Nonfinite or nonpositive energy/weight")
    if not np.allclose(weights, weights[0], rtol=1e-10, atol=0):
        raise ValueError("This fixed-energy spectrum expects a common event weight")
    q_values = energies.sum(axis=1) - proton
    if np.ptp(q_values) > 1e-7:
        raise ValueError("LAB energy conservation failed: sum(T_alpha) - T_p varies")

    # Raw simulated particle counts: each alpha fills one bin once.
    edges = np.linspace(0, 6.5, 261)  # 25 keV per bin
    total, _ = np.histogram(energies.ravel(), edges)
    gs, _ = np.histogram(energies[channel == 1].ravel(), edges)
    exc, _ = np.histogram(energies[channel == 2].ravel(), edges)
    if total.sum() != 3 * n or not np.array_equal(total, gs + exc):
        raise ValueError("Histogram lost entries or channel decomposition failed")

    fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
    ax.stairs(total, edges, color="#222222" if show_components else "black",
              lw=1.7, label="Total")
    if show_components:
        ax.stairs(exc, edges, color="#167D9A", lw=1.2,
                  label=r"via $^8$Be($2^+$)")
        ax.stairs(gs, edges, color="#D48728", lw=1.3,
                  label=r"via $^8$Be(g.s.)")
    ax.set(xlabel=r"$\alpha$ kinetic energy, LAB (MeV)",
           ylabel="Counts",
           xlim=(0, 6.5), ylim=(0, None),
           title=r"$p+{}^{11}$B, $E_p^{\mathrm{lab}}=162$ keV")
    ax.text(0.03, 0.95, "At generation; all three alphas; full solid angle\n"
            f"Input fractions 5% / 95% / 0%; {n:,} reactions",
            transform=ax.transAxes, va="top", fontsize=9)
    if show_components:
        ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(1, 0.80))
    ax.grid(axis="y", alpha=0.15)

    if save or out is not None:
        default_out = DEFAULT_OUT if show_components else DEFAULT_OUT.with_name(DEFAULT_OUT.name + "_total_only")
        out = default_out if out is None else Path(out)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out.with_suffix(".png"), dpi=200)
        fig.savefig(out.with_suffix(".pdf"))
        print(f"Saved {out}.png/.pdf ({n:,} reactions)")
    return fig


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rootfiles", nargs="*", help="Default: the existing 5/95/0 ROOT")
    parser.add_argument("--out", help="Output path stem")
    parser.add_argument("--total-only", action="store_true", help="Only draw the black total curve")
    args = parser.parse_args(argv)
    try:
        from IPython import get_ipython
        in_ipython = get_ipython() is not None
    except ImportError:
        in_ipython = False
    if not in_ipython:
        import matplotlib
        matplotlib.use("Agg")
    fig = plot_spectrum(args.rootfiles or None, args.out, save=True,
                        show_components=not args.total_only)
    if in_ipython:
        # Also support `%run plot_162_lab_total.py` in a notebook.
        from IPython.display import display
        display(fig)
    import matplotlib.pyplot as plt
    plt.close(fig)


if __name__ == "__main__":
    main()
