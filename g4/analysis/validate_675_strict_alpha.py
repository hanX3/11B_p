#!/usr/bin/env python3
"""Validate 675-keV alpha strict-model diagnostics in reaction ROOT files."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np


REQUIRED_FIELDS = [
    "reaction_channel",
    "reaction",
    "resonance_id",
    "branch_id",
    "e_alpha1",
    "e_alpha2",
    "e_alpha3",
    "e_3alpha_cm_alpha1",
    "e_3alpha_cm_alpha2",
    "e_3alpha_cm_alpha3",
    "eaa_8Be",
    "e_alpha8Be",
    "ex_8Be",
    "e_8be_excitation",
    "cos_theta_primary_cm",
    "cos_chi_exit",
    "cos_chi_secondary_8be",
    "opening_angle_alpha12_cm",
    "opening_angle_alpha13_cm",
    "opening_angle_alpha23_cm",
    "h11b675_alpha_decay_model",
    "h11b675_strict_coherent_l13",
    "h11b675_strict_permutation_symmetrized",
    "h11b675_strict_l1_fraction",
    "h11b675_strict_l13_phase",
    "h11b675_strict_8be_lambda_energy_keV",
    "h11b675_strict_8be_reduced_width_squared_keV",
    "h11b675_strict_weight",
    "h11b675_strict_weight_max",
    "h11b675_strict_sampling_attempts",
    "sigma_eval_b",
    "sigma_162_used_b",
    "sigma_675_used_b",
    "sigma_directdecay_b",
    "sigma_total_used_b",
]

MODEL_NAMES = {
    0: "legacyLegendreA2A4",
    1: "symmetrizedCoherentL1L3",
}

CHANNEL_NAMES = {
    0: "n_162_alpha",
    1: "n_675_alpha",
    2: "n_directdecay",
    3: "n_gamma",
}

TOL = 1.0e-10


def load_uproot():
    try:
        import uproot  # type: ignore
    except ImportError as exc:
        raise SystemExit("uproot is not installed; cannot validate ROOT files.") from exc
    return uproot


def load_pyroot():
    try:
        import ROOT  # type: ignore
    except ImportError as exc:
        raise SystemExit("PyROOT is not installed; cannot validate ROOT files with --backend pyroot.") from exc
    return ROOT


def arr(tree, name: str) -> np.ndarray:
    return tree[name].array(library="np")


def finite_stats(values: np.ndarray) -> tuple[float, float, float]:
    values = np.asarray(values)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return (float("nan"), float("nan"), float("nan"))
    return (float(np.min(finite)), float(np.mean(finite)), float(np.max(finite)))


def print_stats(name: str, values: np.ndarray) -> None:
    vmin, vmean, vmax = finite_stats(values)
    print(f"  {name}: min={vmin:.12g} mean={vmean:.12g} max={vmax:.12g}")


def unique_summary(values: np.ndarray) -> str:
    unique, counts = np.unique(values, return_counts=True)
    parts = []
    for value, count in zip(unique, counts):
        if isinstance(value, np.integer):
            label = MODEL_NAMES.get(int(value), "")
            suffix = f" ({label})" if label else ""
            parts.append(f"{int(value)}{suffix}:{int(count)}")
        elif isinstance(value, np.floating):
            parts.append(f"{float(value):.12g}:{int(count)}")
        else:
            parts.append(f"{value}:{int(count)}")
    return ", ".join(parts) if parts else "none"


def count_outside(values: np.ndarray, lo: float, hi: float) -> int:
    values = np.asarray(values)
    finite = np.isfinite(values)
    return int(np.count_nonzero(finite & ((values < lo - TOL) | (values > hi + TOL))))


def load_file(path: Path, backend: str):
    if backend == "uproot":
        uproot = load_uproot()
        with uproot.open(path) as root_file:
            if "tr" not in root_file:
                return None, list(root_file.keys()), 0, {}
            tree = root_file["tr"]
            branches = set(tree.keys())
            data = {
                field: arr(tree, field)
                for field in REQUIRED_FIELDS
                if field != "reaction" and field in branches
            }
            return branches, list(root_file.keys()), tree.num_entries, data

    ROOT = load_pyroot()
    root_file = ROOT.TFile.Open(str(path))
    if not root_file or root_file.IsZombie():
        return None, [], 0, {}
    tree = root_file.Get("tr")
    if not tree:
        keys = [key.GetName() for key in root_file.GetListOfKeys()]
        root_file.Close()
        return None, keys, 0, {}

    branches = {branch.GetName() for branch in tree.GetListOfBranches()}
    data_fields = [field for field in REQUIRED_FIELDS if field != "reaction" and field in branches]
    values = {field: [] for field in data_fields}
    entries = int(tree.GetEntries())
    for index in range(entries):
        tree.GetEntry(index)
        for field in data_fields:
            values[field].append(getattr(tree, field))
    data = {field: np.asarray(field_values) for field, field_values in values.items()}
    root_file.Close()
    return branches, ["tr"], entries, data


def validate_file(path: Path, backend: str) -> int:
    status = 0
    print(f"\n{path}")

    branches, keys, entries, data = load_file(path, backend)
    if branches is None:
        print("  missing tree 'tr'")
        print("  keys: " + ", ".join(keys))
        return 1

    missing = [field for field in REQUIRED_FIELDS if field not in branches]
    print(f"  total entries = {entries}")
    print("  missing fields: " + (", ".join(missing) if missing else "none"))
    if missing:
        return 1
    if entries == 0:
        print("  no reaction entries; checks skipped")
        return 1

    channels = data["reaction_channel"]

    for channel, label in CHANNEL_NAMES.items():
        print(f"  {label} = {int(np.count_nonzero(channels == channel))}")

    print("  reaction_channel values: " + unique_summary(channels))
    print("  resonance_id values: " + unique_summary(data["resonance_id"]))
    print("  branch_id values: " + unique_summary(data["branch_id"]))

    alpha675 = channels == 1
    if not np.any(alpha675):
        print("  no reaction_channel == 1 events; 675-alpha checks skipped")
        return 1

    print("  675-alpha entries = " + str(int(np.count_nonzero(alpha675))))
    print("  model enum mapping: 0=legacyLegendreA2A4, 1=symmetrizedCoherentL1L3")
    for field in (
        "h11b675_alpha_decay_model",
        "h11b675_strict_coherent_l13",
        "h11b675_strict_permutation_symmetrized",
    ):
        print(f"  unique {field}: {unique_summary(data[field][alpha675])}")

    for field in (
        "h11b675_strict_l1_fraction",
        "h11b675_strict_l13_phase",
        "h11b675_strict_8be_lambda_energy_keV",
        "h11b675_strict_8be_reduced_width_squared_keV",
    ):
        print_stats(field, data[field][alpha675])

    strict = alpha675 & (data["h11b675_alpha_decay_model"] == 1)
    if np.any(strict):
        weight = data["h11b675_strict_weight"][strict]
        weight_max = data["h11b675_strict_weight_max"][strict]
        attempts = data["h11b675_strict_sampling_attempts"][strict]
        print_stats("h11b675_strict_weight", weight)
        print_stats("h11b675_strict_weight_max", weight_max)
        print_stats("h11b675_strict_sampling_attempts", attempts)
        print(f"  fraction sampling_attempts > 100 = {np.mean(attempts > 100):.12g}")
        print(f"  fraction sampling_attempts > 1000 = {np.mean(attempts > 1000):.12g}")
        print(f"  fraction sampling_attempts >= 10000 = {np.mean(attempts >= 10000):.12g}")
        print(f"  fraction sampling_attempts >= 3000 = {np.mean(attempts >= 3000):.12g}")
        valid_ratio = np.isfinite(weight) & np.isfinite(weight_max) & (weight_max > 0)
        ratio = np.divide(weight, weight_max, out=np.full_like(weight, np.nan, dtype=float), where=valid_ratio)
        print(f"  fraction weight > weight_max = {np.mean(ratio > 1.0 + TOL):.12g}")
        print_stats("strict weight/weight_max", ratio)
        if np.any(~np.isfinite(weight)) or np.any(weight < -TOL):
            print("  WARN strict weight has non-finite or negative values")
            status = 1
        if np.any(~np.isfinite(weight_max)) or np.any(weight_max <= 0):
            print("  WARN strict weight_max has non-finite or non-positive values")
            status = 1
    else:
        print("  no strict-model 675 events; strict weight checks skipped")

    for field in (
        "e_alpha1",
        "e_alpha2",
        "e_alpha3",
        "e_3alpha_cm_alpha1",
        "e_3alpha_cm_alpha2",
        "e_3alpha_cm_alpha3",
    ):
        values = data[field][alpha675]
        print(f"  {field} negative count = {int(np.count_nonzero(values < -TOL))}")

    e_sum_cm = (
        data["e_3alpha_cm_alpha1"][alpha675]
        + data["e_3alpha_cm_alpha2"][alpha675]
        + data["e_3alpha_cm_alpha3"][alpha675]
    )
    print_stats("e_sum_cm", e_sum_cm)

    for field in ("cos_theta_primary_cm", "cos_chi_exit", "cos_chi_secondary_8be"):
        print(f"  {field} outside [-1,1] = {count_outside(data[field][alpha675], -1.0, 1.0)}")
    for field in (
        "opening_angle_alpha12_cm",
        "opening_angle_alpha13_cm",
        "opening_angle_alpha23_cm",
    ):
        print(f"  {field} outside [0,pi] = {count_outside(data[field][alpha675], 0.0, math.pi)}")

    sigma_sum = data["sigma_162_used_b"] + data["sigma_675_used_b"] + data["sigma_directdecay_b"]
    diff = np.abs(data["sigma_total_used_b"] - sigma_sum)
    print(f"  max |sigma_total_used - (sigma_162 + sigma_675 + sigma_directdecay)| = {float(np.nanmax(diff)):.12g}")
    if np.nanmax(diff) > 1.0e-9:
        status = 1

    return status


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=("uproot", "pyroot"), default="uproot")
    parser.add_argument("root_files", nargs="+", type=Path)
    args = parser.parse_args()

    status = 0
    for root_file in args.root_files:
        status |= validate_file(root_file, args.backend)
    return status


if __name__ == "__main__":
    sys.exit(main())
