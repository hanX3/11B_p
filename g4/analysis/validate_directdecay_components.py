#!/usr/bin/env python3
"""Validate direct-decay and scaled cross-section diagnostics in reaction ROOT files."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


REQUIRED_BRANCHES = [
    "sigma_eval_b",
    "sigma_162_model_b",
    "sigma_162_total_b",
    "sigma_162_used_b",
    "sigma_675_total_b",
    "sigma_675_used_b",
    "sigma_total_used_b",
    "sigma_directdecay_b",
    "sigma_directdecay_sampling_b",
    "direct_decay_fraction",
    "direct_decay_fraction_162",
    "direct_decay_fraction_675",
    "sequential_decay_fraction_162",
    "sequential_decay_fraction_675",
    "enable_direct_decay",
    "scale_factor_162",
    "scale_factor_675",
    "channel_probability_directdecay",
    "reaction_channel",
    "event_weight",
    "event_sampling_weight",
    "reaction",
]

OPTIONAL_ALIAS_BRANCHES = [
    "sigma_background_b",
    "sigma_background_sampling_b",
    "background_bias_factor",
    "channel_probability_background",
]

TOL = 1.0e-10


def load_uproot():
    try:
        import uproot  # type: ignore
    except ImportError as exc:
        raise SystemExit("uproot is not installed; cannot validate ROOT files.") from exc
    return uproot


def max_abs(values: np.ndarray) -> float:
    if values.size == 0:
        return float("nan")
    return float(np.nanmax(np.abs(values)))


def branch_array(tree, name: str) -> np.ndarray:
    return tree[name].array(library="np")


def decode_reactions(values: np.ndarray) -> np.ndarray:
    decoded = []
    for value in values:
        if isinstance(value, bytes):
            decoded.append(value.decode("utf-8", errors="replace").rstrip("\x00"))
        else:
            decoded.append(str(value).rstrip("\x00"))
    return np.asarray(decoded, dtype=object)


def summarize_channels(channels: np.ndarray) -> dict[int, int]:
    return {channel: int(np.count_nonzero(channels == channel)) for channel in (0, 1, 2, 3)}


def print_metric(name: str, value: float, tol: float | None = None) -> None:
    suffix = ""
    if tol is not None:
        suffix = " OK" if np.isfinite(value) and value < tol else " WARN"
    print(f"  {name} = {value:.12g}{suffix}")


def validate_file(path: Path) -> int:
    uproot = load_uproot()
    status = 0

    with uproot.open(path) as root_file:
        if "tr" not in root_file:
            print(f"\n{path}: missing tree 'tr'")
            return 1
        tree = root_file["tr"]
        branches = set(tree.keys())

        print(f"\n{path}")
        print(f"  entries = {tree.num_entries}")

        missing = [name for name in REQUIRED_BRANCHES if name not in branches]
        if missing:
            print("  missing required branches: " + ", ".join(missing))
            return 1

        alias_present = [name for name in OPTIONAL_ALIAS_BRANCHES if name in branches]
        alias_missing = [name for name in OPTIONAL_ALIAS_BRANCHES if name not in branches]
        print("  alias branches present: " + (", ".join(alias_present) if alias_present else "none"))
        print("  alias branches absent: " + (", ".join(alias_missing) if alias_missing else "none"))

        if tree.num_entries == 0:
            print("  no reaction entries; numerical checks skipped")
            return 1

        data = {name: branch_array(tree, name) for name in REQUIRED_BRANCHES if name != "reaction"}
        reactions = decode_reactions(branch_array(tree, "reaction"))

        for name in (
            "sigma_eval_b",
            "sigma_162_model_b",
            "sigma_162_total_b",
            "sigma_162_used_b",
            "sigma_675_total_b",
            "sigma_675_used_b",
            "sigma_total_used_b",
            "sigma_directdecay_b",
            "sigma_directdecay_sampling_b",
            "direct_decay_fraction",
            "direct_decay_fraction_162",
            "direct_decay_fraction_675",
            "sequential_decay_fraction_162",
            "sequential_decay_fraction_675",
            "scale_factor_162",
            "scale_factor_675",
            "channel_probability_directdecay",
            "event_weight",
            "event_sampling_weight",
        ):
            values = data[name]
            print(
                f"  mean {name} = {float(np.nanmean(values)):.12g}, "
                f"min = {float(np.nanmin(values)):.12g}, max = {float(np.nanmax(values)):.12g}"
            )

        if np.any(data["sigma_directdecay_b"] < -TOL):
            print("  sigma_directdecay_b has negative values WARN")
            status = 1
        if np.any(data["sigma_directdecay_sampling_b"] < -TOL):
            print("  sigma_directdecay_sampling_b has negative values WARN")
            status = 1

        sigma_162_expected = data["sequential_decay_fraction_162"] * data["sigma_162_total_b"]
        diff_162 = max_abs(data["sigma_162_used_b"] - sigma_162_expected)
        print_metric("max |sigma_162_used - expected|", diff_162, TOL)
        if diff_162 >= TOL:
            status = 1

        sigma_675_expected = data["sequential_decay_fraction_675"] * data["sigma_675_total_b"]
        diff_675 = max_abs(data["sigma_675_used_b"] - sigma_675_expected)
        print_metric("max |sigma_675_used - expected|", diff_675, TOL)
        if diff_675 >= TOL:
            status = 1

        sigma_direct_expected = np.where(
            data["enable_direct_decay"] != 0,
            (1.0 - data["sequential_decay_fraction_162"]) * data["sigma_162_total_b"]
            + (1.0 - data["sequential_decay_fraction_675"]) * data["sigma_675_total_b"],
            0.0,
        )
        direct_fraction_expected = np.divide(
            sigma_direct_expected,
            data["sigma_total_used_b"],
            out=np.zeros_like(sigma_direct_expected),
            where=data["sigma_total_used_b"] > 0,
        )
        diff_direct = max_abs(data["sigma_directdecay_b"] - sigma_direct_expected)
        print_metric("max |sigma_directdecay - expected|", diff_direct, TOL)

        sigma_sum = (
            data["sigma_162_used_b"]
            + data["sigma_675_used_b"]
            + data["sigma_directdecay_b"]
        )
        diff_total = max_abs(data["sigma_total_used_b"] - sigma_sum)
        print_metric("max |sigma_total_used - sum|", diff_total, TOL)
        if diff_total >= TOL:
            status = 1

        channels = data["reaction_channel"]
        counts = summarize_channels(channels)
        total = int(channels.size)
        print(f"  n_total_events = {total}")
        for channel, label in (
            (0, "n_162"),
            (1, "n_675"),
            (2, "n_directdecay"),
            (3, "n_gamma"),
        ):
            fraction = counts[channel] / total if total else float("nan")
            print(f"  {label} = {counts[channel]}, fraction = {fraction:.8g}")

        direct_reactions = reactions[channels == 2]
        unique_direct = sorted(set(direct_reactions.tolist()))
        if unique_direct:
            print("  directdecay reaction strings: " + ", ".join(unique_direct))
        else:
            print("  directdecay reaction strings: none")
        legacy_strings = [text for text in unique_direct if "background3alpha" in text]
        if legacy_strings:
            print("  legacy background reaction string still present WARN")
            status = 1

        diff_direct_fraction = max_abs(data["direct_decay_fraction"] - direct_fraction_expected)
        print_metric("max |direct_decay_fraction - expected|", diff_direct_fraction, TOL)
        if diff_direct_fraction >= TOL:
            status = 1
        diff_direct_fraction_162 = max_abs(
            data["direct_decay_fraction_162"] - (1.0 - data["sequential_decay_fraction_162"])
        )
        print_metric("max |direct_decay_fraction_162 - expected|", diff_direct_fraction_162, TOL)
        if diff_direct_fraction_162 >= TOL:
            status = 1
        diff_direct_fraction_675 = max_abs(
            data["direct_decay_fraction_675"] - (1.0 - data["sequential_decay_fraction_675"])
        )
        print_metric("max |direct_decay_fraction_675 - expected|", diff_direct_fraction_675, TOL)
        if diff_direct_fraction_675 >= TOL:
            status = 1

        optional_data = {name: branch_array(tree, name) for name in alias_present}
        if "sigma_background_b" in optional_data:
            diff = max_abs(optional_data["sigma_background_b"] - data["sigma_directdecay_b"])
            print_metric("max |sigma_background - sigma_directdecay|", diff, TOL)
            if diff >= TOL:
                status = 1
        if "sigma_background_sampling_b" in optional_data:
            diff = max_abs(optional_data["sigma_background_sampling_b"] - data["sigma_directdecay_sampling_b"])
            print_metric("max |sigma_background_sampling - sigma_directdecay_sampling|", diff, TOL)
            if diff >= TOL:
                status = 1
        if "channel_probability_background" in optional_data:
            diff = max_abs(optional_data["channel_probability_background"] - data["channel_probability_directdecay"])
            print_metric("max |prob_background - prob_directdecay|", diff, TOL)
            if diff >= TOL:
                status = 1
        if "background_bias_factor" in optional_data:
            diff = max_abs(optional_data["background_bias_factor"] - data["direct_decay_fraction"])
            print_metric("max |background_bias_factor - direct_decay_fraction|", diff, TOL)
            if diff >= TOL:
                status = 1

    return status


def pyroot_string(view) -> str:
    return bytes(view[:128]).split(b"\0", 1)[0].decode("utf-8", errors="replace")


def validate_file_pyroot(path: Path) -> int:
    try:
        import ROOT  # type: ignore
    except ImportError as exc:
        raise SystemExit("PyROOT is not installed; cannot use --backend pyroot.") from exc

    status = 0
    root_file = ROOT.TFile.Open(str(path))
    if not root_file or root_file.IsZombie():
        print(f"\n{path}: could not open ROOT file")
        return 1
    tree = root_file.Get("tr")
    if not tree:
        print(f"\n{path}: missing tree 'tr'")
        return 1

    branches = {branch.GetName() for branch in tree.GetListOfBranches()}
    print(f"\n{path}")
    print(f"  entries = {tree.GetEntries()}")

    missing = [name for name in REQUIRED_BRANCHES if name not in branches]
    if missing:
        print("  missing required branches: " + ", ".join(missing))
        return 1

    alias_present = [name for name in OPTIONAL_ALIAS_BRANCHES if name in branches]
    alias_missing = [name for name in OPTIONAL_ALIAS_BRANCHES if name not in branches]
    print("  alias branches present: " + (", ".join(alias_present) if alias_present else "none"))
    print("  alias branches absent: " + (", ".join(alias_missing) if alias_missing else "none"))

    if tree.GetEntries() == 0:
        print("  no reaction entries; numerical checks skipped")
        return 1

    numeric_names = [name for name in REQUIRED_BRANCHES if name != "reaction"]
    data_lists = {name: [] for name in numeric_names}
    optional_lists = {name: [] for name in alias_present}

    for index in range(tree.GetEntries()):
        tree.GetEntry(index)
        for name in numeric_names:
            data_lists[name].append(getattr(tree, name))
        for name in alias_present:
            optional_lists[name].append(getattr(tree, name))

    data = {name: np.asarray(values) for name, values in data_lists.items()}
    optional_data = {name: np.asarray(values) for name, values in optional_lists.items()}

    for name in (
        "sigma_eval_b",
        "sigma_162_model_b",
        "sigma_162_total_b",
        "sigma_162_used_b",
        "sigma_675_total_b",
        "sigma_675_used_b",
        "sigma_total_used_b",
        "sigma_directdecay_b",
        "sigma_directdecay_sampling_b",
        "direct_decay_fraction",
        "direct_decay_fraction_162",
        "direct_decay_fraction_675",
        "sequential_decay_fraction_162",
        "sequential_decay_fraction_675",
        "scale_factor_162",
        "scale_factor_675",
        "channel_probability_directdecay",
        "event_weight",
        "event_sampling_weight",
    ):
        values = data[name]
        print(
            f"  mean {name} = {float(np.nanmean(values)):.12g}, "
            f"min = {float(np.nanmin(values)):.12g}, max = {float(np.nanmax(values)):.12g}"
        )

    if np.any(data["sigma_directdecay_b"] < -TOL):
        print("  sigma_directdecay_b has negative values WARN")
        status = 1
    if np.any(data["sigma_directdecay_sampling_b"] < -TOL):
        print("  sigma_directdecay_sampling_b has negative values WARN")
        status = 1

    sigma_162_expected = data["sequential_decay_fraction_162"] * data["sigma_162_total_b"]
    diff_162 = max_abs(data["sigma_162_used_b"] - sigma_162_expected)
    print_metric("max |sigma_162_used - expected|", diff_162, TOL)
    if diff_162 >= TOL:
        status = 1

    sigma_675_expected = data["sequential_decay_fraction_675"] * data["sigma_675_total_b"]
    diff_675 = max_abs(data["sigma_675_used_b"] - sigma_675_expected)
    print_metric("max |sigma_675_used - expected|", diff_675, TOL)
    if diff_675 >= TOL:
        status = 1

    sigma_direct_expected = np.where(
        data["enable_direct_decay"] != 0,
        (1.0 - data["sequential_decay_fraction_162"]) * data["sigma_162_total_b"]
        + (1.0 - data["sequential_decay_fraction_675"]) * data["sigma_675_total_b"],
        0.0,
    )
    direct_fraction_expected = np.divide(
        sigma_direct_expected,
        data["sigma_total_used_b"],
        out=np.zeros_like(sigma_direct_expected),
        where=data["sigma_total_used_b"] > 0,
    )
    diff_direct = max_abs(data["sigma_directdecay_b"] - sigma_direct_expected)
    print_metric("max |sigma_directdecay - expected|", diff_direct, TOL)
    if diff_direct >= TOL:
        status = 1

    diff_direct_fraction = max_abs(data["direct_decay_fraction"] - direct_fraction_expected)
    print_metric("max |direct_decay_fraction - expected|", diff_direct_fraction, TOL)
    if diff_direct_fraction >= TOL:
        status = 1

    diff_direct_fraction_162 = max_abs(
        data["direct_decay_fraction_162"] - (1.0 - data["sequential_decay_fraction_162"])
    )
    print_metric("max |direct_decay_fraction_162 - expected|", diff_direct_fraction_162, TOL)
    if diff_direct_fraction_162 >= TOL:
        status = 1

    diff_direct_fraction_675 = max_abs(
        data["direct_decay_fraction_675"] - (1.0 - data["sequential_decay_fraction_675"])
    )
    print_metric("max |direct_decay_fraction_675 - expected|", diff_direct_fraction_675, TOL)
    if diff_direct_fraction_675 >= TOL:
        status = 1

    sigma_sum = data["sigma_162_used_b"] + data["sigma_675_used_b"] + data["sigma_directdecay_b"]
    diff_total = max_abs(data["sigma_total_used_b"] - sigma_sum)
    print_metric("max |sigma_total_used - sum|", diff_total, TOL)
    if diff_total >= TOL:
        status = 1

    channels = data["reaction_channel"]
    counts = summarize_channels(channels)
    total = int(channels.size)
    print(f"  n_total_events = {total}")
    for channel, label in (
        (0, "n_162"),
        (1, "n_675"),
        (2, "n_directdecay"),
        (3, "n_gamma"),
    ):
        fraction = counts[channel] / total if total else float("nan")
        print(f"  {label} = {counts[channel]}, fraction = {fraction:.8g}")

    print("  reaction string branch present; PyROOT backend skips char-array contents")

    if "sigma_background_b" in optional_data:
        diff = max_abs(optional_data["sigma_background_b"] - data["sigma_directdecay_b"])
        print_metric("max |sigma_background - sigma_directdecay|", diff, TOL)
        if diff >= TOL:
            status = 1
    if "sigma_background_sampling_b" in optional_data:
        diff = max_abs(optional_data["sigma_background_sampling_b"] - data["sigma_directdecay_sampling_b"])
        print_metric("max |sigma_background_sampling - sigma_directdecay_sampling|", diff, TOL)
        if diff >= TOL:
            status = 1
    if "channel_probability_background" in optional_data:
        diff = max_abs(optional_data["channel_probability_background"] - data["channel_probability_directdecay"])
        print_metric("max |prob_background - prob_directdecay|", diff, TOL)
        if diff >= TOL:
            status = 1
    if "background_bias_factor" in optional_data:
        diff = max_abs(optional_data["background_bias_factor"] - data["direct_decay_fraction"])
        print_metric("max |background_bias_factor - direct_decay_fraction|", diff, TOL)
        if diff >= TOL:
            status = 1

    root_file.Close()
    return status


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--backend",
        choices=("uproot", "pyroot"),
        default="uproot",
        help="ROOT reader backend. Default uses uproot as requested by the validation task.",
    )
    parser.add_argument("root_files", nargs="+", type=Path)
    args = parser.parse_args()

    status = 0
    for path in args.root_files:
        if not path.exists():
            print(f"{path}: file not found")
            status = 1
            continue
        if args.backend == "pyroot":
            status = max(status, validate_file_pyroot(path))
        else:
            status = max(status, validate_file(path))
    return status


if __name__ == "__main__":
    sys.exit(main())
