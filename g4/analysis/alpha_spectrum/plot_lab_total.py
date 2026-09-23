#!/usr/bin/env python3
"""Count all three alphas at generation in the LAB frame (25 keV bins).

Terminal: python3 analysis/alpha_spectrum/plot_lab_total.py [FILE.root ...]
Total curve only: python3 analysis/alpha_spectrum/plot_lab_total.py --total-only
Separate totals and overview: python3 analysis/alpha_spectrum/plot_lab_total.py --energies 50 200 10 --separate --overview
Jupyter:  spectra = load_spectra(); fig = plot_spectrum(spectra=spectra)
Pass the combined ROOT or the individual energy-point files, never both.
Only PNG/PDF are saved. Imports preserve the notebook's Matplotlib backend.
"""

import argparse
from pathlib import Path

import numpy as np
import uproot

HERE = Path(__file__).resolve().parent
DEFAULT_ROOT = HERE.parents[1] / "data/spectrum_150_300_step10_lab_sum.root"
DEFAULT_OUT = HERE / "alpha_lab_150_300_sum"
EDGES = np.linspace(0, 7, 281)
CHANNELS = (1, 2, 4)
COLORS = ("#D48728", "#167D9A", "#8D5AA5")
LABELS = (r"162 resonance: $^8$Be(g.s.)",
          r"162 resonance: $^8$Be($2^+$)",
          r"675 resonance: $^8$Be($2^+$)")


def load_spectra(rootfiles=None):
    """Read one point at a time; return histograms and per-point statistics.

    ScanInfo gives the incident energy and disjoint entry range of each run.
    Legacy single-energy ROOT files without ScanInfo are also supported.
    """
    if rootfiles is None:
        rootfiles = [DEFAULT_ROOT]
    elif isinstance(rootfiles, (str, Path)):
        rootfiles = [rootfiles]
    paths = [Path(p).resolve() for p in rootfiles]
    if not paths or len(set(paths)) != len(paths):
        raise ValueError("Provide distinct input ROOT files")
    fields = ["event", "e_alpha1", "e_alpha2", "e_alpha3", "projectile_kinetic_lab",
              "resonance_id", "h11b_reaction_channel", "event_weight",
              "configured_162_alpha0_branching_fraction", "enable_direct_decay",
              "sigma_162_total_barn", "sigma_675_total_barn",
              "h11b675_strict_sampling_attempts", "h11b675_strict_weight",
              "h11b675_strict_weight_max"]
    points = []
    seen = set()
    q_min, q_max = np.inf, -np.inf
    common_weight = None
    for path in paths:
        # Avoid the installed fsspec async reader, which stalls on local baskets.
        with uproot.open(path, handler=uproot.source.file.MemmapSource) as f:
            tree = f["reaction"]
            if "ScanInfo" in f:
                info = f["ScanInfo"].arrays(library="np")
                ranges = list(zip(info["energy_keV"], info["reaction_entries"],
                                  info["incident_events"], info["random_seed"]))
                if sum(int(r[1]) for r in ranges) != tree.num_entries:
                    raise ValueError("ScanInfo does not cover the reaction tree")
            else:
                seed = np.unique(f["RunInfo"]["random_seed"].array(library="np"))
                if len(seed) != 1:
                    raise ValueError("Multiple runs require ScanInfo entry ranges")
                ranges = [(None, tree.num_entries, None, seed[0])]
            start = 0
            for energy, count, incident, seed in ranges:
                count = int(count)
                if count <= 0:
                    raise ValueError("Empty reaction sample")
                key = (int(seed), None if energy is None else int(energy))
                if key in seen:
                    raise ValueError("The same run appears in multiple inputs")
                seen.add(key)
                a = tree.arrays(fields, entry_start=start, entry_stop=start + count, library="np")
                start += count
                events = a["event"]
                if len(np.unique(events)) != count:
                    raise ValueError("Duplicate event IDs within an energy point")
                if incident is not None and (count > int(incident) or events.min() < 0
                                              or events.max() >= int(incident)):
                    raise ValueError("Reaction IDs exceed the incident event count")
                alpha = np.column_stack([a[f"e_alpha{i}"] for i in (1, 2, 3)])
                proton = a["projectile_kinetic_lab"]
                weights = a["event_weight"]
                if (not np.all(np.isfinite(alpha)) or np.any(alpha <= 0)
                        or not np.all(np.isfinite(proton)) or np.any(proton <= 0)
                        or not np.all(np.isfinite(weights)) or np.any(weights <= 0)):
                    raise ValueError("Invalid energy or event weight")
                if common_weight is None:
                    common_weight = weights[0]
                if not np.allclose(weights, common_weight, rtol=1e-10, atol=0):
                    raise ValueError("Raw-count comparison expects a common event bias")
                if energy is None:
                    energy = int(round(proton.mean() * 1000))
                energy = int(energy)
                if proton.max() * 1000 > energy + 1e-6:
                    raise ValueError("Reaction energies disagree with the beam-energy label")
                q = alpha.sum(axis=1) - proton
                q_min = min(q_min, float(q.min()))
                q_max = max(q_max, float(q.max()))
                channel = a["h11b_reaction_channel"]
                if (not np.all(np.isin(channel, CHANNELS))
                        or not np.all(a["enable_direct_decay"] == 0)
                        or not np.allclose(a["configured_162_alpha0_branching_fraction"], .05)):
                    raise ValueError("Expected sequential 162/675 channels with 162 gs fraction 5%")
                if not np.all(a["resonance_id"] == np.where(channel == 4, 675, 162)):
                    raise ValueError("Resonance and channel labels disagree")
                strict = channel == 4
                if (np.any(a["h11b675_strict_sampling_attempts"][strict] >= 10000)
                        or np.any(a["h11b675_strict_weight"][strict] > a["h11b675_strict_weight_max"][strict])):
                    raise ValueError("675 rejection sampler exceeded its limits")
                histogram = np.array([np.histogram(alpha[channel == c].ravel(), EDGES)[0]
                                      for c in CHANNELS])
                if histogram.sum() != count * 3:
                    raise ValueError("Histogram energy range lost alpha entries")
                s162, s675 = a["sigma_162_total_barn"], a["sigma_675_total_barn"]
                if np.any(s162 + s675 <= 0):
                    raise ValueError("Invalid total cross section")
                p162 = s162 / (s162 + s675)
                expected = np.array([.05 * p162.mean(), .95 * p162.mean(), 1 - p162.mean()])
                points.append(dict(energy=energy, reactions=count, incident=incident,
                                   seed=int(seed), histogram=histogram,
                                   channels=np.array([np.count_nonzero(channel == c) for c in CHANNELS]),
                                   expected=expected, proton_mean=proton.mean() * 1000,
                                   proton_min=proton.min() * 1000,
                                   mean_alpha=float(alpha.mean()),
                                   sigma_mb=np.array([s162.mean(), s675.mean()]) * 1000,
                                   alpha_below_1=int(np.count_nonzero(alpha < 1)),
                                   alpha_above_5=int(np.count_nonzero(alpha > 5))))
    if q_max - q_min > 1e-7:
        raise ValueError("LAB energy conservation failed across input runs")
    points.sort(key=lambda p: p["energy"])
    if len({p["energy"] for p in points}) != len(points):
        raise ValueError("Multiple samples for one energy would change its relative weight")
    return dict(points=points, edges=EDGES.copy(), q_range=(q_min, q_max), weight=common_weight)


def save_figure(fig, out):
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".png"), dpi=200)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"Saved {out}.png/.pdf")


def plot_spectrum(rootfiles=None, out=None, *, save=False, spectra=None, show_components=True):
    """Return the summed spectrum; show_components=False draws only the total."""
    import matplotlib.pyplot as plt
    spectra = load_spectra(rootfiles) if spectra is None else spectra
    points = spectra["points"]
    histogram = sum(p["histogram"] for p in points)
    fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
    ax.stairs(histogram.sum(axis=0), EDGES, color="black", lw=1.6, label="Total")
    if show_components:
        for counts, color, label in zip(histogram, COLORS, LABELS):
            if counts.sum():
                ax.stairs(counts, EDGES, color=color, lw=1.1, label=label)
    if len(points) == 1:
        title = f"$p+^{{11}}$B, $E_p^{{lab}}={points[0]['energy']}$ keV"
    else:
        energies = [p["energy"] for p in points]
        title = f"$p+^{{11}}$B, LAB: {energies[0]}–{energies[-1]} keV, {len(points)} energy points"
    ax.set(xlabel=r"$\alpha$ kinetic energy, LAB (MeV)", ylabel="Counts",
           title=title, xlim=(0, 7), ylim=(0, None))
    if show_components:
        ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(1, .80), fontsize=9)
    ax.grid(axis="y", alpha=.15)
    if save or out is not None:
        default_out = DEFAULT_OUT if show_components else DEFAULT_OUT.with_name(DEFAULT_OUT.name + "_total_only")
        save_figure(fig, default_out if out is None else out)
    return fig


def plot_individual_spectra(spectra, out_dir=None, *, save=False):
    """Return {energy_keV: Figure}, one black total per point, each with its own y scale."""
    figures = {}
    out_dir = HERE if out_dir is None else Path(out_dir)
    for point in spectra["points"]:
        fig = plot_spectrum(spectra={**spectra, "points": [point]}, show_components=False)
        if save:
            save_figure(fig, out_dir / f"alpha_lab_{point['energy']}keV_total")
        figures[point["energy"]] = fig
    return figures


def plot_energy_points(spectra, *, show_components=True, share_y=True):
    """Return an energy-point grid; disable components and shared y for separate totals."""
    import matplotlib.pyplot as plt
    points = spectra["points"]
    columns = min(4, len(points))
    rows = (len(points) + columns - 1) // columns
    fig, axes = plt.subplots(rows, columns, figsize=(14, 2.6 * rows),
                             sharex=True, sharey=share_y, squeeze=False, layout="constrained")
    ymax = 1.08 * max(p["histogram"].sum(axis=0).max() for p in points)
    for ax, point in zip(axes.ravel(), points):
        hist = point["histogram"]
        ax.stairs(hist.sum(axis=0), EDGES, color="black", lw=1.1)
        if show_components:
            for counts, color in zip(hist, COLORS):
                ax.stairs(counts, EDGES, color=color, lw=.8, alpha=.85)
        local_ymax = ymax if share_y else 1.08 * hist.sum(axis=0).max()
        ax.set(title=f"{point['energy']} keV", xlim=(0, 7), ylim=(0, local_ymax))
        ax.ticklabel_format(axis="y", style="plain", useOffset=False)
        ax.tick_params(labelsize=9)
        ax.grid(axis="y", alpha=.15)
    for ax in axes.ravel()[len(points):]:
        ax.set_visible(False)
    fig.supxlabel(r"$\alpha$ kinetic energy, LAB (MeV)")
    fig.supylabel("Counts")
    fig.suptitle(r"$p+{}^{11}$B, " + f"{points[0]['energy']}–{points[-1]['energy']} keV")
    return fig


def plot_channel_fractions(spectra):
    """Observed reaction fractions and expected fractions from the runtime model."""
    import matplotlib.pyplot as plt
    points = spectra["points"]
    energy = np.array([p["energy"] for p in points])
    observed = np.array([p["channels"] / p["reactions"] for p in points])
    expected = np.array([p["expected"] for p in points])
    fig, ax = plt.subplots(figsize=(8, 4.5), layout="constrained")
    for i, (color, label) in enumerate(zip(COLORS, LABELS)):
        ax.plot(energy, 100 * expected[:, i], color=color, lw=1.2, label=label)
        ax.scatter(energy, 100 * observed[:, i], color=color, s=18, zorder=3)
    ax.set(xlabel="Incident proton energy, LAB (keV)", ylabel="Reaction fraction (%)",
           ylim=(0, 100), title="Channel composition: lines = model; points = simulation")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=.15)
    return fig


def print_summary(spectra):
    print("Ep (keV)  Reactions     162-gs    162-exc    675-exc   675 (%)   Mean alpha (MeV)")
    for p in spectra["points"]:
        gs, exc162, exc675 = p["channels"]
        print(f"{p['energy']:8d} {p['reactions']:10,d} {gs:10,d} {exc162:10,d} {exc675:10,d} "
              f"{exc675 / p['reactions'] * 100:9.4f} {p['mean_alpha']:18.6f}")
    total = sum(p["reactions"] for p in spectra["points"])
    print(f"Total: {total:,} reactions, {3 * total:,} alpha entries")
    print(f"Q-value spread: {np.ptp(spectra['q_range']):.3g} MeV")
    print("Counts are the unweighted sum of the simulated samples, not absolute physical yields.")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rootfiles", nargs="*")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--total-only", action="store_true", help="Only draw the summed alpha curve")
    parser.add_argument("--energies", type=int, nargs=3, metavar=("START", "STOP", "STEP"),
                        help="Read individual data/spectrum_*keV_lab_162_675.root files")
    parser.add_argument("--separate", action="store_true", help="Only draw individual black total spectra")
    parser.add_argument("--overview", action="store_true", help="Draw a grid of black total spectra")
    parser.add_argument("--out-dir", type=Path, help="Directory for --separate/--overview PNG/PDF outputs")
    args = parser.parse_args(argv)
    if args.energies:
        if args.rootfiles:
            parser.error("Choose ROOT files or --energies, not both")
        first, last, step = args.energies
        if first <= 0 or last < first or step <= 0 or (last - first) % step:
            parser.error("Use positive energies/step and a stop reachable from start in whole steps")
        args.rootfiles = [HERE.parents[1] / "data" / f"spectrum_{energy}keV_lab_162_675.root"
                          for energy in range(first, last + 1, step)]
    if (args.separate or args.overview) and args.out is not None:
        parser.error("Use --out-dir for separate plots or the overview")
    if args.out_dir is not None and not (args.separate or args.overview):
        parser.error("--out-dir requires --separate or --overview")
    try:
        from IPython import get_ipython
        in_ipython = get_ipython() is not None
    except ImportError:
        in_ipython = False
    if not in_ipython:
        import matplotlib
        matplotlib.use("Agg")
    spectra = load_spectra(args.rootfiles or None)
    print_summary(spectra)
    if args.separate or args.overview:
        figures = list(plot_individual_spectra(spectra, args.out_dir, save=True).values()) if args.separate else []
        if args.overview:
            fig = plot_energy_points(spectra, show_components=False, share_y=False)
            first, last = spectra["points"][0]["energy"], spectra["points"][-1]["energy"]
            out_dir = HERE if args.out_dir is None else args.out_dir
            save_figure(fig, out_dir / f"alpha_lab_{first}_{last}_overview")
            figures.append(fig)
        import matplotlib.pyplot as plt
        for fig in figures:
            if in_ipython:
                from IPython.display import display
                display(fig)
            plt.close(fig)
        return
    out = args.out or (DEFAULT_OUT.with_name(DEFAULT_OUT.name + "_total_only")
                       if args.total_only else DEFAULT_OUT)
    figures = [plot_spectrum(spectra=spectra, out=out, show_components=not args.total_only)]
    if not args.total_only:
        figures.append(plot_spectrum(spectra=spectra, show_components=False,
                                     out=out.with_name(out.name + "_total_only")))
    if len(spectra["points"]) > 1 and not args.total_only:
        for suffix, plot in [("_by_energy", plot_energy_points), ("_channels", plot_channel_fractions)]:
            fig = plot(spectra)
            save_figure(fig, out.with_name(out.name + suffix))
            figures.append(fig)
    import matplotlib.pyplot as plt
    for fig in figures:
        if in_ipython:
            from IPython.display import display
            display(fig)
        plt.close(fig)


if __name__ == "__main__":
    main()
