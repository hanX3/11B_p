#!/usr/bin/env python3
"""Optimize two forward-backward-symmetric additional polar-angle gaps.

The fixed baseline blind regions are, by default,

    [0, 20], [80, 100], [160, 180] degrees.

For each requested additional missing width M, k is fixed to 2. The two new
blind gaps are constrained to be mirror images about 90 degrees, and each has
width

    additional gap width = M / 2.

For a forward gap [a, a+w], the backward gap is

    [180-a-w, 180-a].

The additional gaps may touch but may not overlap the fixed baseline blind
regions. Every candidate forward position on the requested angular grid is
evaluated, so the result is the exact optimum on that discrete grid. Event
loss is evaluated from complete three-alpha events: an event is accepted only
when all three alphas avoid every baseline and additional blind region.

The default robust objective maximizes the minimum efficiency among the three
pure channels. A weighted objective is also available.

Outputs:
  * terminal table with optimized blind and covered intervals;
  * efficiency-versus-additional-loss PNG;
  * optimal angular-coverage interval PNG.

This is a standalone script. No companion analysis module is required.
No JSON file or Python bytecode cache is written.

Usage
-----
Terminal:
    python3 optimize_162_theta_coverage.py \
        gs.root be2plus.root direct.root
    python3 optimize_162_theta_coverage.py \
        gs.root be2plus.root direct.root --objective weighted --show

Jupyter / IPython (figures render inline, nothing is closed).  The three ROOT
files may be passed positionally, in the order gs / 2+ / direct:
    from optimize_162_theta_coverage import optimize_162_theta_coverage

    fig_eff, fig_cov = optimize_162_theta_coverage(
        "data/benchmark_162seq_Be_gs.root",
        "data/benchmark_162seq_Be2plus.root",
        "data/benchmark_162direct.root",
        output_dir="analysis/efficiency",
    )
    # any CLI flag is a keyword (dashes -> underscores):
    fig_eff, fig_cov = optimize_162_theta_coverage(
        seq_be_gs="gs.root", seq_be2plus="be2plus.root", direct="direct.root",
        angle_step=0.5, objective="weighted", output_dir="figures",
    )

Figure names default to
    angular_gap_efficiency_162.png
    angular_gap_coverage_162.png
and are written into --output-dir / output_dir= (default: the current working
directory), which is created if it does not exist.
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path
from typing import Any

sys.dont_write_bytecode = True
_MPL_TEMP_DIR = tempfile.TemporaryDirectory(prefix="optimize-k2-matplotlib-")
os.environ.setdefault("MPLCONFIGDIR", _MPL_TEMP_DIR.name)

import matplotlib


def _in_notebook() -> bool:
    """True only inside a Jupyter/IPython ZMQ kernel."""
    try:
        from IPython import get_ipython

        shell = get_ipython()
        return shell is not None and shell.__class__.__name__ == "ZMQInteractiveShell"
    except Exception:
        return False


IN_NOTEBOOK = _in_notebook()

# Backend selection (works in a head-less terminal AND in Jupyter):
#   * Jupyter / IPython        -> leave the inline (or user) backend untouched;
#   * terminal with a display  -> keep the interactive backend so --show works;
#   * head-less terminal       -> force Agg so savefig still works.
# savefig() works on every backend, so the PNGs are always written.
if (
    not IN_NOTEBOOK
    and os.name != "nt"
    and not os.environ.get("DISPLAY")
    and not os.environ.get("MPLBACKEND")
):
    matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402  (import after backend choice)
import numpy as np  # noqa: E402


def _require_uproot():
    """Import uproot on demand.

    Deferred so that this module can be imported (and the geometry helpers
    used) without uproot installed, and so that a notebook import does not
    die with SystemExit.
    """
    try:
        import uproot
    except ImportError as exc:  # pragma: no cover - environment check
        raise ImportError(
            "Missing dependency 'uproot'. Install with "
            "'python3 -m pip install uproot numpy matplotlib'."
        ) from exc
    return uproot


ALPHA_PDG = 1000020040
ALPHA_MASS_MEV = 3727.379378
H11B_REACTION_PRODUCT = 1

CHANNELS = (
    ("162seq_Be_gs", r"162 seq. $^{8}$Be(g.s.)", "#0072B2"),
    ("162seq_Be2plus", r"162 seq. $^{8}$Be($2^{+}$)", "#D55E00"),
    ("162direct", "162 direct", "#009E73"),
)

SPHERE_BRANCHES = (
    "event_id",
    "pdg",
    "kinetic_energy_MeV",
    "px_MeV_c",
    "py_MeV_c",
    "pz_MeV_c",
    "h11b_particle_source",
    "generator_particle_index",
)


# Fixed axes rectangle shared by BOTH figures.
#
# tight_layout() sizes the margins from each figure's own labels, so the two
# plots ended up with different axes boxes (the coverage plot's "44.8%" tick
# labels are wider than the efficiency plot's "20/40/60").  Pinning the same
# rectangle in both, and saving without bbox_inches="tight", guarantees the
# plot frames are pixel-identical and line up when stacked on a slide.
# Font size of the angle labels drawn inside the coverage bars.  Kept small
# enough that the narrowest bar (e.g. "100-130", 30 deg wide) still fits.
BAR_LABEL_FONTSIZE = 10.0

AXES_RECT = {
    "left": 0.150,
    "right": 0.930,
    "top": 0.930,
    "bottom": 0.230,
}


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 14,
            "axes.labelsize": 15,
            "axes.titlesize": 15,
            "legend.fontsize": 14,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "axes.linewidth": 1.1,
            "lines.linewidth": 1.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )


def object_names(root_file: Any) -> set[str]:
    return {str(name).split(";")[0] for name in root_file.keys()}


def branch_names(tree: Any) -> set[str]:
    return {str(name).split(";")[0] for name in tree.keys()}


def complete_three_alpha_rows(
    event: np.ndarray,
    generator_index: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, dict[int, int], int]:
    """Return rows for events containing generator indices 0, 1 and 2."""
    order = np.lexsort((generator_index, event))
    event_sorted = event[order]
    index_sorted = generator_index[order]

    unique_event, first, multiplicity = np.unique(
        event_sorted,
        return_index=True,
        return_counts=True,
    )
    mult_values, mult_counts = np.unique(
        multiplicity,
        return_counts=True,
    )
    multiplicity_distribution = {
        int(value): int(count)
        for value, count in zip(mult_values, mult_counts)
    }

    three = multiplicity == 3
    starts = first[three]
    candidate_events = unique_event[three]
    if starts.size == 0:
        return (
            np.empty(0, dtype=np.int64),
            np.empty((0, 3), dtype=np.int64),
            multiplicity_distribution,
            0,
        )

    offsets = np.arange(3, dtype=np.int64)
    sorted_positions = starts[:, None] + offsets[None, :]
    expected_indices = np.array(
        [0, 1, 2],
        dtype=index_sorted.dtype,
    )
    correct_indices = np.all(
        index_sorted[sorted_positions] == expected_indices,
        axis=1,
    )
    wrong_indices = int(np.count_nonzero(~correct_indices))
    rows = order[sorted_positions[correct_indices]]

    return (
        np.asarray(
            candidate_events[correct_indices],
            dtype=np.int64,
        ),
        np.asarray(rows, dtype=np.int64),
        multiplicity_distribution,
        wrong_indices,
    )


def load_channel_cosines(
    root_path: Path,
    mass_shell_tolerance: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Load quality-selected complete three-alpha momentum directions."""
    if not root_path.is_file():
        raise RuntimeError(f"Input file not found: {root_path}")

    uproot = _require_uproot()
    with uproot.open(root_path) as root_file:
        if "virtual_sphere" not in object_names(root_file):
            raise RuntimeError(
                f"{root_path}: missing 'virtual_sphere' tree"
            )

        tree = root_file["virtual_sphere"]
        missing = sorted(
            set(SPHERE_BRANCHES) - branch_names(tree)
        )
        if missing:
            raise RuntimeError(
                f"{root_path}: virtual_sphere is missing branches: "
                f"{', '.join(missing)}"
            )

        sphere = tree.arrays(
            SPHERE_BRANCHES,
            library="np",
        )

    selected = (
        (np.asarray(sphere["pdg"]) == ALPHA_PDG)
        & (
            np.asarray(sphere["h11b_particle_source"])
            == H11B_REACTION_PRODUCT
        )
        & (
            np.asarray(sphere["generator_particle_index"])
            >= 0
        )
        & (
            np.asarray(sphere["generator_particle_index"])
            <= 2
        )
    )
    if not np.any(selected):
        raise RuntimeError(
            f"{root_path}: no selected H11B alpha records"
        )

    event = np.asarray(
        sphere["event_id"][selected],
        dtype=np.int64,
    )
    generator_index = np.asarray(
        sphere["generator_particle_index"][selected],
        dtype=np.int32,
    )
    kinetic = np.asarray(
        sphere["kinetic_energy_MeV"][selected],
        dtype=float,
    )
    momentum = np.column_stack(
        [
            np.asarray(
                sphere["px_MeV_c"][selected],
                dtype=float,
            ),
            np.asarray(
                sphere["py_MeV_c"][selected],
                dtype=float,
            ),
            np.asarray(
                sphere["pz_MeV_c"][selected],
                dtype=float,
            ),
        ]
    )

    (
        complete_event,
        rows,
        multiplicities,
        wrong_indices,
    ) = complete_three_alpha_rows(
        event,
        generator_index,
    )
    if rows.size == 0:
        raise RuntimeError(
            f"{root_path}: no complete three-alpha events"
        )

    event_kinetic = kinetic[rows]
    event_momentum = momentum[rows]

    finite = np.all(
        np.isfinite(event_kinetic),
        axis=1,
    ) & np.all(
        np.isfinite(event_momentum),
        axis=(1, 2),
    )
    nonnegative = np.all(
        event_kinetic >= 0.0,
        axis=1,
    )

    momentum2 = np.sum(
        event_momentum**2,
        axis=2,
    )
    expected_momentum2 = event_kinetic * (
        event_kinetic + 2.0 * ALPHA_MASS_MEV
    )
    mass_shell_residual = np.abs(
        momentum2 - expected_momentum2
    ) / np.maximum(
        expected_momentum2,
        1.0,
    )
    mass_shell = (
        np.max(
            mass_shell_residual,
            axis=1,
        )
        <= mass_shell_tolerance
    )

    momentum_norm = np.sqrt(momentum2)
    angular = np.all(
        momentum_norm > 0.0,
        axis=1,
    )
    accepted = (
        finite
        & nonnegative
        & mass_shell
        & angular
    )
    if not np.any(accepted):
        raise RuntimeError(
            f"{root_path}: no events survive quality requirements"
        )

    accepted_momentum = event_momentum[accepted]
    accepted_norm = momentum_norm[accepted]
    accepted_kinetic = event_kinetic[accepted]

    cos_theta = np.clip(
        accepted_momentum[:, :, 2] / accepted_norm,
        -1.0,
        1.0,
    )

    diagnostics = {
        "selected_records": int(np.count_nonzero(selected)),
        "multiplicities": multiplicities,
        "wrong_indices": wrong_indices,
        "n_complete": int(complete_event.size),
        "n_finite": int(np.count_nonzero(finite)),
        "n_nonnegative": int(
            np.count_nonzero(finite & nonnegative)
        ),
        "n_mass_shell": int(
            np.count_nonzero(
                finite & nonnegative & mass_shell
            )
        ),
        "n_angular": int(np.count_nonzero(accepted)),
    }
    return cos_theta, accepted_kinetic, diagnostics


POPCOUNT = np.array(
    [int(value).bit_count() for value in range(256)],
    dtype=np.uint8,
)


def normalize_intervals(
    values: list[float] | tuple[float, ...],
) -> list[tuple[float, float]]:
    """Validate, sort and merge LOW HIGH polar-angle pairs in [0, 180]."""
    if len(values) % 2 != 0:
        raise ValueError(
            "angle intervals must be supplied as LOW HIGH pairs"
        )

    intervals: list[tuple[float, float]] = []
    for index in range(0, len(values), 2):
        start = float(values[index])
        stop = float(values[index + 1])
        if not np.isfinite(start) or not np.isfinite(stop):
            raise ValueError("angle interval bounds must be finite")
        if start < 0.0 or stop > 180.0 or stop <= start:
            raise ValueError(
                "each angle interval must satisfy "
                "0 <= LOW < HIGH <= 180 degrees"
            )
        intervals.append((start, stop))

    intervals.sort()
    merged: list[tuple[float, float]] = []
    tolerance = 1.0e-10
    for start, stop in intervals:
        if not merged or start > merged[-1][1] + tolerance:
            merged.append((start, stop))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], stop))
    return merged


def validate_forward_backward_symmetry(
    intervals: list[tuple[float, float]],
) -> None:
    """Require the complete fixed blind layout to be symmetric about 90 deg."""
    mirrored = sorted(
        (180.0 - stop, 180.0 - start)
        for start, stop in intervals
    )
    if len(intervals) != len(mirrored):
        raise ValueError("internal interval-symmetry validation failure")

    for original, mirror in zip(intervals, mirrored):
        if not np.allclose(original, mirror, atol=1.0e-9, rtol=0.0):
            raise ValueError(
                "--base-blind-theta must be forward-backward symmetric "
                "about 90 degrees"
            )


def merge_intervals(
    intervals: list[tuple[float, float]],
) -> list[tuple[float, float]]:
    flat = [value for interval in intervals for value in interval]
    return normalize_intervals(flat)


def complement_intervals(
    blind_intervals: list[tuple[float, float]],
) -> list[tuple[float, float]]:
    """Return covered polar-angle intervals inside [0, 180] degrees."""
    covered: list[tuple[float, float]] = []
    cursor = 0.0
    tolerance = 1.0e-10
    for start, stop in blind_intervals:
        if start > cursor + tolerance:
            covered.append((cursor, start))
        cursor = max(cursor, stop)
    if cursor < 180.0 - tolerance:
        covered.append((cursor, 180.0))
    return covered


def interval_width(
    intervals: list[tuple[float, float]],
) -> float:
    return float(sum(stop - start for start, stop in intervals))


def intervals_overlap(
    first: tuple[float, float],
    second: tuple[float, float],
) -> bool:
    tolerance = 1.0e-9
    return (
        first[0] < second[1] - tolerance
        and second[0] < first[1] - tolerance
    )


def candidate_pair(
    forward_start: float,
    width_deg: float,
) -> list[tuple[float, float]]:
    """Return a forward gap and its mirror about 90 degrees."""
    forward = (forward_start, forward_start + width_deg)
    backward = (
        180.0 - forward_start - width_deg,
        180.0 - forward_start,
    )
    return [forward, backward]


def candidate_forward_starts(
    width_deg: float,
    step_deg: float,
    base_blind: list[tuple[float, float]],
) -> np.ndarray:
    """Return all valid unique forward-gap starts on the discrete grid."""
    maximum = 90.0 - width_deg
    if maximum < -1.0e-10:
        return np.empty(0, dtype=float)

    starts = np.arange(
        0.0,
        max(maximum, 0.0) + 1.0e-10,
        step_deg,
        dtype=float,
    )
    if starts.size == 0 or maximum - starts[-1] > 1.0e-8:
        starts = np.append(starts, maximum)

    # Include positions that exactly touch fixed blind-region boundaries,
    # even when a user-selected grid does not pass through those positions.
    special: list[float] = [0.0, maximum]
    for start, stop in base_blind:
        special.extend(
            [
                start,
                stop,
                start - width_deg,
                stop - width_deg,
                180.0 - start - width_deg,
                180.0 - stop - width_deg,
            ]
        )

    starts = np.unique(
        np.round(
            np.concatenate(
                [
                    starts,
                    np.clip(np.asarray(special, dtype=float), 0.0, maximum),
                ]
            ),
            decimals=10,
        )
    )

    valid: list[float] = []
    for start in starts:
        pair = candidate_pair(float(start), width_deg)
        if intervals_overlap(pair[0], pair[1]):
            continue
        if any(
            intervals_overlap(extra, fixed)
            for extra in pair
            for fixed in base_blind
        ):
            continue
        valid.append(float(start))

    return np.asarray(valid, dtype=float)


def event_hits_intervals(
    theta: np.ndarray,
    intervals: list[tuple[float, float]],
) -> np.ndarray:
    """Return one Boolean value per event: any alpha enters any interval."""
    event_hit = np.zeros(theta.shape[0], dtype=bool)
    for start, stop in intervals:
        if stop >= 180.0 - 1.0e-10:
            inside = (theta >= start) & (theta <= stop)
        else:
            inside = (theta >= start) & (theta < stop)
        event_hit |= np.any(inside, axis=1)
    return event_hit


def packed_event_hits(
    theta: np.ndarray,
    intervals: list[tuple[float, float]],
) -> np.ndarray:
    return np.packbits(
        event_hits_intervals(theta, intervals),
        bitorder="little",
    )


def objective_from_efficiencies(
    efficiencies: np.ndarray,
    objective: str,
    weights: np.ndarray,
) -> float:
    if objective == "robust":
        return float(np.min(efficiencies))
    return float(np.dot(weights, efficiencies))


def efficiencies_from_packed_unions(
    unions_by_channel: list[np.ndarray],
    event_counts: np.ndarray,
) -> np.ndarray:
    efficiencies = np.empty(3, dtype=float)
    for channel, (packed, n_events) in enumerate(
        zip(unions_by_channel, event_counts)
    ):
        lost = int(np.sum(POPCOUNT[packed], dtype=np.int64))
        efficiencies[channel] = 1.0 - lost / float(n_events)
    return efficiencies


def baseline_result(
    base_packed_by_channel: list[np.ndarray],
    event_counts: np.ndarray,
    base_blind: list[tuple[float, float]],
    objective: str,
    weights: np.ndarray,
) -> dict[str, object]:
    efficiencies = efficiencies_from_packed_unions(
        base_packed_by_channel,
        event_counts,
    )
    return {
        "additional_missing_deg": 0.0,
        "additional_gap_width_deg": 0.0,
        "forward_start_deg": None,
        "additional_gaps": [],
        "blind_intervals": base_blind,
        "covered_intervals": complement_intervals(base_blind),
        "efficiencies": efficiencies,
        "score": objective_from_efficiencies(
            efficiencies,
            objective,
            weights,
        ),
        "candidate_count": 0,
    }


def optimize_additional_missing(
    additional_missing_deg: float,
    angle_step_deg: float,
    theta_by_channel: list[np.ndarray],
    event_counts: np.ndarray,
    base_blind: list[tuple[float, float]],
    base_packed_by_channel: list[np.ndarray],
    objective: str,
    weights: np.ndarray,
) -> dict[str, object]:
    """Optimize one symmetric k=2 additional-gap configuration."""
    width_deg = additional_missing_deg / 2.0
    starts = candidate_forward_starts(
        width_deg,
        angle_step_deg,
        base_blind,
    )
    if starts.size == 0:
        raise RuntimeError(
            f"no valid symmetric k=2 placement exists for "
            f"{additional_missing_deg:g} deg additional missing width"
        )

    best_score = -np.inf
    best_start: float | None = None
    best_efficiencies: np.ndarray | None = None

    for start in starts:
        extra_gaps = candidate_pair(float(start), width_deg)
        candidate_packed_by_channel = [
            packed_event_hits(theta, extra_gaps)
            for theta in theta_by_channel
        ]
        unions = [
            np.bitwise_or(base, extra)
            for base, extra in zip(
                base_packed_by_channel,
                candidate_packed_by_channel,
            )
        ]
        efficiencies = efficiencies_from_packed_unions(
            unions,
            event_counts,
        )
        score = objective_from_efficiencies(
            efficiencies,
            objective,
            weights,
        )

        # Deterministic tie handling: keep the smaller forward start.
        if score > best_score + 1.0e-12:
            best_score = score
            best_start = float(start)
            best_efficiencies = efficiencies

    if best_start is None or best_efficiencies is None:
        raise RuntimeError("symmetric k=2 optimization failed")

    additional_gaps = candidate_pair(best_start, width_deg)
    blind_intervals = merge_intervals(base_blind + additional_gaps)
    covered_intervals = complement_intervals(blind_intervals)

    return {
        "additional_missing_deg": additional_missing_deg,
        "additional_gap_width_deg": width_deg,
        "forward_start_deg": best_start,
        "additional_gaps": additional_gaps,
        "blind_intervals": blind_intervals,
        "covered_intervals": covered_intervals,
        "efficiencies": best_efficiencies,
        "score": best_score,
        "candidate_count": int(starts.size),
    }


def interval_text(
    intervals: list[tuple[float, float]],
    decimals: int = 2,
) -> str:
    if not intervals:
        return "none"
    return ", ".join(
        f"[{start:.{decimals}f},{stop:.{decimals}f}]"
        for start, stop in intervals
    )


def plot_efficiency(
    results: list[dict[str, object]],
    output: Path | None = None,
) -> "plt.Figure":
    additional = np.asarray(
        [result["additional_missing_deg"] for result in results],
        dtype=float,
    )
    efficiencies = np.stack(
        [np.asarray(result["efficiencies"]) for result in results]
    )

    fig, axis = plt.subplots(figsize=(6, 3.6))
    for channel, (_, label, color) in enumerate(CHANNELS):
        axis.plot(
            additional,
            100.0 * efficiencies[:, channel],
            marker="o",
            color=color,
            label=label,
        )

    axis.set_xticks(additional)
    axis.set_xlabel("Additional symmetric missing polar width (deg)")
    axis.set_ylabel(r"Triple-$\alpha$ efficiency (%)")
    axis.grid(alpha=0.25)
    axis.legend(frameon=False, loc="best")
    fig.subplots_adjust(**AXES_RECT)
    if output is not None:
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=300)
    return fig


def plot_coverage_intervals(
    results: list[dict[str, object]],
    output: Path | None = None,
) -> "plt.Figure":
    """Plot retained angular coverage directly, one row per loss scenario."""
    y_positions = np.arange(len(results), dtype=float)
    row_labels: list[str] = []

    fig, axis = plt.subplots(
        figsize=(6, 3.6)
    )

    for row, result in enumerate(results):
        covered = result["covered_intervals"]
        efficiencies = np.asarray(result["efficiencies"])
        minimum_efficiency = float(np.min(efficiencies))

        for start_deg, stop_deg in covered:
            axis.hlines(
                y=row,
                xmin=start_deg,
                xmax=stop_deg,
                linewidth=22.0,
            )
            interval_width_deg = stop_deg - start_deg
            if interval_width_deg >= 12.0:
                # Label drawn inside the bar, in white for contrast.
                axis.text(
                    0.5 * (start_deg + stop_deg),
                    row,
                    f"{start_deg:g}–{stop_deg:g}°",
                    ha="center",
                    va="center",
                    fontsize=BAR_LABEL_FONTSIZE,
                    color="white",
                    fontweight="bold",
                    zorder=3,
                )

        row_labels.append(f"{100.0 * minimum_efficiency:.1f}%")

    axis.set_xlim(0.0, 180.0)
    axis.set_ylim(-0.75, len(results) - 0.25)
    axis.set_xticks(np.arange(0.0, 181.0, 20.0))
    axis.set_yticks(y_positions, row_labels)
    axis.invert_yaxis()
    axis.set_xlabel(
        r"Retained laboratory polar-angle coverage "
        r"$\theta^{\mathrm{lab}}$ (deg)"
    )
    axis.grid(axis="x", alpha=0.25)
    axis.tick_params(axis="y", length=0)
    fig.subplots_adjust(**AXES_RECT)
    if output is not None:
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=300)
    return fig


def resolve_output_path(
    output_dir: Path | str,
    filename: Path | str,
) -> Path:
    """Join filename onto output_dir, unless filename is already absolute."""
    path = Path(filename)
    if path.is_absolute():
        return path
    return Path(output_dir) / path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Optimize fixed-k=2 forward-backward-symmetric additional "
            "polar-angle gaps for 162-keV triple-alpha efficiency."
        )
    )
    # The three ROOT files are positional on the command line but optional in
    # the parser, so that a notebook can call parse_args([]) and supply them
    # as keywords instead.  They are validated in the main function.
    parser.add_argument(
        "seq_be_gs",
        nargs="?",
        default=None,
        type=Path,
        help="pure 162seq_Be_gs ROOT file",
    )
    parser.add_argument(
        "seq_be2plus",
        nargs="?",
        default=None,
        type=Path,
        help="pure 162seq_Be2plus ROOT file",
    )
    parser.add_argument(
        "direct",
        nargs="?",
        default=None,
        type=Path,
        help="pure 162direct ROOT file",
    )
    parser.add_argument(
        "--base-blind-theta",
        type=float,
        nargs="+",
        default=(0.0, 20.0, 80.0, 100.0, 160.0, 180.0),
        metavar="ANGLE",
        help=(
            "fixed forward-backward-symmetric blind LOW HIGH pairs; "
            "default: 0 20 80 100 160 180"
        ),
    )
    parser.add_argument(
        "--additional-missing",
        type=float,
        nargs="+",
        default=(10.0, 20.0, 30.0, 40.0, 50.0, 60.0),
        metavar="DEG",
        help=(
            "additional total missing widths. k is fixed to 2, so each "
            "new symmetric gap has half this width; default: "
            "10 20 30 40 50 60"
        ),
    )
    parser.add_argument(
        "--angle-step",
        type=float,
        default=1.0,
        help=(
            "forward-gap start grid in degrees; every valid grid position "
            "is evaluated exactly (default: 1)"
        ),
    )
    parser.add_argument(
        "--objective",
        choices=("robust", "weighted"),
        default="robust",
        help=(
            "robust=maximize minimum channel efficiency; "
            "weighted=maximize assumed-mixture efficiency"
        ),
    )
    parser.add_argument(
        "--weights",
        type=float,
        nargs=3,
        default=(0.2, 0.6, 0.2),
        metavar=("GS", "BE2PLUS", "DIRECT"),
        help=(
            "channel weights for --objective weighted; "
            "default: 0.2 0.6 0.2"
        ),
    )
    parser.add_argument(
        "--mass-shell-tolerance",
        type=float,
        default=1.0e-4,
    )
    parser.add_argument(
        "--efficiency-output",
        type=Path,
        default=Path(
            "angular_gap_efficiency_162.png"
        ),
    )
    parser.add_argument(
        "--coverage-output",
        "--layout-output",
        dest="coverage_output",
        type=Path,
        default=Path(
            "angular_gap_coverage_162.png"
        ),
        help=(
            "output PNG for the retained angular-coverage interval plot; "
            "--layout-output is retained as a compatibility alias"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        metavar="DIR",
        help=(
            "directory the figures are written to; created if missing "
            "(default: current working directory). Ignored for any output "
            "given as an absolute path."
        ),
    )
    parser.add_argument(
        "--no-save",
        action="store_true",
        help="skip writing the PNGs (useful in notebooks)",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="open the figures in a window (interactive backends only)",
    )

    # In a notebook, sys.argv carries the kernel's own flags (e.g. -f ...),
    # so fall back to defaults there and ignore any stray/unknown flags.
    if argv is None:
        argv = [] if IN_NOTEBOOK else sys.argv[1:]
    args, _unknown = parser.parse_known_args(argv)
    return args


def optimize_162_theta_coverage(
    seq_be_gs: Path | str | None = None,
    seq_be2plus: Path | str | None = None,
    direct: Path | str | None = None,
    *,
    argv: list[str] | None = None,
    **overrides,
) -> tuple["plt.Figure", "plt.Figure"]:
    """Run the optimization and return (efficiency_figure, coverage_figure).

    Terminal:
        python3 optimize_162_theta_coverage.py \
            gs.root be2plus.root direct.root --show

    Jupyter / IPython (figures are shown inline).  The three ROOT files may be
    given positionally, in the order gs / 2+ / direct:

        from optimize_162_theta_coverage import optimize_162_theta_coverage

        fig_eff, fig_cov = optimize_162_theta_coverage(
            "data/benchmark_162seq_Be_gs.root",
            "data/benchmark_162seq_Be2plus.root",
            "data/benchmark_162direct.root",
            output_dir="analysis/efficiency",
        )

    or as keywords:

        fig_eff, fig_cov = optimize_162_theta_coverage(
            seq_be_gs="gs.root",
            seq_be2plus="be2plus.root",
            direct="direct.root",
        )

    Any command-line option is available as a keyword with dashes replaced by
    underscores, e.g. --angle-step becomes angle_step, --objective becomes
    objective.  By default both PNGs are written into output_dir (default: the
    current working directory, created if missing) under their standard names;
    pass no_save=True to skip writing them.
    """
    args = parse_args(argv)

    # The three inputs may arrive positionally or as keywords; Python itself
    # rejects supplying the same one both ways.
    explicit: dict[str, object] = {}
    for key, value in (
        ("seq_be_gs", seq_be_gs),
        ("seq_be2plus", seq_be2plus),
        ("direct", direct),
    ):
        if value is not None:
            explicit[key] = value
    explicit.update(overrides)

    path_keys = {"seq_be_gs", "seq_be2plus", "direct",
                 "efficiency_output", "coverage_output", "output_dir"}
    for key, value in explicit.items():
        if not hasattr(args, key):
            raise TypeError(f"unknown option: {key!r}")
        if key in path_keys and value is not None:
            value = Path(value)
        setattr(args, key, value)

    missing_inputs = [
        name
        for name in ("seq_be_gs", "seq_be2plus", "direct")
        if getattr(args, name) is None
    ]
    if missing_inputs:
        raise ValueError(
            "the three ROOT inputs are required; missing: "
            + ", ".join(missing_inputs)
            + ".  Pass them positionally in the order gs / 2+ / direct, "
            "e.g. optimize_162_theta_coverage('gs.root', 'be2plus.root', "
            "'direct.root'), or as keywords seq_be_gs=, seq_be2plus=, "
            "direct=."
        )

    for path in (args.seq_be_gs, args.seq_be2plus, args.direct):
        if not path.is_file():
            raise ValueError(f"Input file not found: {path}")

    if not np.isfinite(args.angle_step) or args.angle_step <= 0.0:
        raise ValueError("--angle-step must be finite and positive")

    additional_missing = np.unique(
        np.asarray(args.additional_missing, dtype=float)
    )
    if (
        additional_missing.size == 0
        or np.any(~np.isfinite(additional_missing))
        or np.any(additional_missing <= 0.0)
        or np.any(additional_missing >= 180.0)
    ):
        raise ValueError(
            "--additional-missing values must be finite and between "
            "0 and 180 degrees"
        )

    weights = np.asarray(args.weights, dtype=float)
    if (
        np.any(~np.isfinite(weights))
        or np.any(weights < 0.0)
        or np.sum(weights) <= 0.0
    ):
        raise ValueError(
            "--weights must be finite, non-negative and have a positive sum"
        )
    weights /= np.sum(weights)

    base_blind = normalize_intervals(
        list(args.base_blind_theta)
    )
    validate_forward_backward_symmetry(base_blind)
    base_width = interval_width(base_blind)

    if np.any(base_width + additional_missing >= 180.0 - 1.0e-10):
        raise ValueError(
            "fixed baseline width plus every additional missing width "
            "must remain below 180 degrees"
        )

    configure_style()
    paths = (
        args.seq_be_gs,
        args.seq_be2plus,
        args.direct,
    )
    theta_by_channel: list[np.ndarray] = []
    event_counts: list[int] = []

    print("Loading quality-selected complete three-alpha events...")
    for path, (key, _, _) in zip(paths, CHANNELS):
        cos_theta, _, diagnostics = load_channel_cosines(
            path,
            args.mass_shell_tolerance,
        )
        theta = np.rad2deg(
            np.arccos(np.clip(cos_theta, -1.0, 1.0))
        )
        theta_by_channel.append(theta)
        event_counts.append(theta.shape[0])
        print(
            f"  {key:<22} accepted={theta.shape[0]:,}, "
            f"complete={diagnostics['n_complete']:,}"
        )

    event_counts_array = np.asarray(
        event_counts,
        dtype=np.int64,
    )
    base_packed_by_channel = [
        packed_event_hits(theta, base_blind)
        for theta in theta_by_channel
    ]

    print("\nFixed-k symmetric gap optimization")
    print("=" * 160)
    print("k=2, forward-backward symmetry required")
    print(
        f"fixed baseline blind={interval_text(base_blind)}, "
        f"baseline missing width={base_width:g} deg"
    )
    print(
        f"additional missing widths="
        f"{', '.join(f'{value:g}' for value in additional_missing)} deg, "
        f"objective={args.objective}, "
        f"grid step={args.angle_step:g} deg"
    )
    if args.objective == "weighted":
        print(f"normalized weights={weights}")
    print("-" * 160)

    results: list[dict[str, object]] = [
        baseline_result(
            base_packed_by_channel,
            event_counts_array,
            base_blind,
            args.objective,
            weights,
        )
    ]

    for missing in additional_missing:
        results.append(
            optimize_additional_missing(
                float(missing),
                args.angle_step,
                theta_by_channel,
                event_counts_array,
                base_blind,
                base_packed_by_channel,
                args.objective,
                weights,
            )
        )

    print(
        f"{'extra':>7} {'gap width':>10} "
        f"{'g.s. eff.':>12} {'2+ eff.':>12} "
        f"{'direct eff.':>12} {'objective':>12} "
        f"{'candidates':>11}"
    )
    print("-" * 160)

    for result in results:
        efficiencies = np.asarray(result["efficiencies"])
        print(
            f"{result['additional_missing_deg']:7.1f} "
            f"{result['additional_gap_width_deg']:10.1f} "
            f"{efficiencies[0]:12.6%} "
            f"{efficiencies[1]:12.6%} "
            f"{efficiencies[2]:12.6%} "
            f"{result['score']:12.6%} "
            f"{result['candidate_count']:11d}"
        )
        print(
            "    optimized extra blind: "
            f"{interval_text(result['additional_gaps'])}"
        )
        print(
            "    total blind:           "
            f"{interval_text(result['blind_intervals'])}"
        )
        print(
            "    optimal coverage:      "
            f"{interval_text(result['covered_intervals'])}"
        )

    print("=" * 160)

    efficiency_path = resolve_output_path(
        args.output_dir,
        args.efficiency_output,
    )
    coverage_path = resolve_output_path(
        args.output_dir,
        args.coverage_output,
    )

    efficiency_figure = plot_efficiency(
        results,
        None if args.no_save else efficiency_path,
    )
    coverage_figure = plot_coverage_intervals(
        results,
        None if args.no_save else coverage_path,
    )

    if args.no_save:
        print("Figures not written to disk (--no-save / no_save=True).")
    else:
        print("Wrote:")
        print(f"  {efficiency_path}")
        print(f"  {coverage_path}")
    print(
        "The reported result is the exact optimum on the specified "
        "one-dimensional forward-start grid."
    )

    # --show works on every backend: interactive backends open a window; on
    # Agg it is a harmless no-op, so say so rather than failing silently.
    interactive = matplotlib.get_backend().lower() != "agg"
    if args.show and interactive:
        plt.show()
    elif args.show and not interactive:
        print(
            "[info] --show requested but a non-interactive (Agg) backend is "
            "active; the figures were saved to disk instead."
        )

    return efficiency_figure, coverage_figure


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point."""
    try:
        figures = optimize_162_theta_coverage(argv=argv)
    except (ValueError, RuntimeError, ImportError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    for figure in figures:
        plt.close(figure)
    return 0


if __name__ == "__main__":
    if IN_NOTEBOOK:
        # `%run` inside Jupyter: show inline, keep the kernel alive.
        optimize_162_theta_coverage()
        plt.show()
    else:
        raise SystemExit(main())
