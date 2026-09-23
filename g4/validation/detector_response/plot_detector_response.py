#!/usr/bin/env python3
"""Make detector-response validation plots from ROOT outputs."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from validate_detector_response import RootTree, DETECTOR_NAMES


def load(path: str | None, branches: list[str]) -> tuple[RootTree | None, dict[str, np.ndarray]]:
    if not path:
        return None, {}
    tree = RootTree(Path(path))
    return tree, tree.arrays(branches)


def savefig(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def map_channels(reaction: dict[str, np.ndarray], events: np.ndarray) -> np.ndarray:
    if "event" not in reaction or "reaction_channel" not in reaction:
        return np.full(events.shape, -1, dtype=int)
    lookup = {int(e): int(c) for e, c in zip(reaction["event"], reaction["reaction_channel"])}
    return np.asarray([lookup.get(int(e), -1) for e in events], dtype=int)


def plot_event_energy(outdir: Path, label: str, event: dict[str, np.ndarray], reaction: dict[str, np.ndarray]) -> list[Path]:
    paths = []
    if not event or "e" not in event or "detector_type" not in event:
        print(f"SKIP {label}: event energy plot needs e and detector_type")
        return paths
    plt.figure(figsize=(7.5, 5))
    for det_type in sorted(set(event["detector_type"].astype(int).tolist())):
        mask = event["detector_type"].astype(int) == det_type
        values = event["e"][mask]
        if values.size:
            plt.hist(values, bins=80, histtype="step", label=f"{det_type} {DETECTOR_NAMES.get(det_type, 'unknown')}")
    plt.xlabel("Detector hit energy (MeV)")
    plt.ylabel("Hits")
    plt.yscale("log")
    plt.title(f"{label}: detector hit energy by detector type")
    plt.legend()
    path = outdir / f"{label}_edep_by_detector.png"
    savefig(path)
    paths.append(path)

    if "event" in event and reaction:
        channels = map_channels(reaction, event["event"])
        plt.figure(figsize=(7.5, 5))
        for channel in sorted(set(channels.tolist())):
            if channel < 0:
                continue
            values = event["e"][channels == channel]
            if values.size:
                plt.hist(values, bins=80, histtype="step", label=f"channel {channel}")
        plt.xlabel("Detector hit energy (MeV)")
        plt.ylabel("Hits")
        plt.yscale("log")
        plt.title(f"{label}: detector hit energy by reaction channel")
        plt.legend()
        path = outdir / f"{label}_edep_by_channel.png"
        savefig(path)
        paths.append(path)
    return paths


def plot_hit_multiplicity(outdir: Path, label: str, event: dict[str, np.ndarray], reaction: dict[str, np.ndarray]) -> list[Path]:
    paths = []
    if not event or "event" not in event:
        print(f"SKIP {label}: hit multiplicity plot needs event branch")
        return paths
    counts = Counter(event["event"].astype(int).tolist())
    values = np.asarray(list(counts.values()), dtype=float)
    plt.figure(figsize=(7.5, 5))
    bins = np.arange(0.5, max(values.max() if values.size else 1, 1) + 1.5, 1)
    plt.hist(values, bins=bins, histtype="stepfilled", alpha=0.75)
    plt.xlabel("Above-threshold detector hits per hit event")
    plt.ylabel("Events")
    plt.title(f"{label}: total hit multiplicity")
    path = outdir / f"{label}_hit_multiplicity_total.png"
    savefig(path)
    paths.append(path)

    if "detector_type" in event:
        plt.figure(figsize=(7.5, 5))
        for det_type in sorted(set(event["detector_type"].astype(int).tolist())):
            mask = event["detector_type"].astype(int) == det_type
            det_counts = Counter(event["event"][mask].astype(int).tolist())
            det_values = np.asarray(list(det_counts.values()), dtype=float)
            if det_values.size:
                bins = np.arange(0.5, max(det_values.max(), 1) + 1.5, 1)
                plt.hist(det_values, bins=bins, histtype="step", label=f"{det_type} {DETECTOR_NAMES.get(det_type, 'unknown')}")
        plt.xlabel("Hits per hit event")
        plt.ylabel("Events")
        plt.title(f"{label}: hit multiplicity by detector type")
        plt.legend()
        path = outdir / f"{label}_hit_multiplicity_by_detector.png"
        savefig(path)
        paths.append(path)
    return paths


def plot_positions(outdir: Path, label: str, event: dict[str, np.ndarray]) -> list[Path]:
    paths = []
    required = {"x", "y", "detector_type"}
    if not event or not required.issubset(event):
        print(f"SKIP {label}: first-hit XY plot needs x/y/detector_type")
        return paths
    det = event["detector_type"].astype(int)
    for name, mask in [("si", det == 1), ("gamma_detector", np.isin(det, [2, 3]))]:
        if not mask.any():
            continue
        plt.figure(figsize=(6, 6))
        plt.scatter(event["x"][mask], event["y"][mask], s=6, alpha=0.45)
        plt.xlabel("x (Geant4 length unit)")
        plt.ylabel("y (Geant4 length unit)")
        plt.title(f"{label}: hit positions, {name}")
        plt.axis("equal")
        path = outdir / f"{label}_first_hit_xy_{name}.png"
        savefig(path)
        paths.append(path)
    return paths


def plot_step(outdir: Path, label: str, step_tree: RootTree | None, step: dict[str, np.ndarray]) -> list[Path]:
    paths = []
    if not step:
        print(f"SKIP {label}: no step tree")
        return paths
    if "de" in step:
        positive = step["de"][step["de"] > 0]
        plt.figure(figsize=(7.5, 5))
        if positive.size:
            plt.hist(positive, bins=100, histtype="step")
            plt.yscale("log")
            plt.xscale("log")
        plt.xlabel("Step energy deposit (MeV)")
        plt.ylabel("Steps")
        plt.title(f"{label}: positive step edep")
        path = outdir / f"{label}_step_edep.png"
        savefig(path)
        paths.append(path)
    if "length" in step:
        positive = step["length"][step["length"] > 0]
        plt.figure(figsize=(7.5, 5))
        if positive.size:
            plt.hist(positive, bins=100, histtype="step")
            plt.yscale("log")
            plt.xscale("log")
        plt.xlabel("Step length (Geant4 length unit)")
        plt.ylabel("Steps")
        plt.title(f"{label}: step length distribution")
        path = outdir / f"{label}_step_length.png"
        savefig(path)
        paths.append(path)
    if step_tree is not None and "volume" in step_tree.branches:
        volumes = step_tree.strings("volume", limit=250000)
        counts = Counter(volumes.tolist()).most_common(20)
        if counts:
            labels = [k for k, _ in counts][::-1]
            values = [v for _, v in counts][::-1]
            plt.figure(figsize=(8, 6))
            plt.barh(labels, values)
            plt.xlabel("Sampled step count")
            plt.title(f"{label}: most common step volumes")
            path = outdir / f"{label}_volume_step_counts.png"
            savefig(path)
            paths.append(path)
    return paths


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", required=True)
    parser.add_argument("--reaction", required=True)
    parser.add_argument("--event", required=True)
    parser.add_argument("--track", default=None)
    parser.add_argument("--step", default=None)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    _, reaction = load(args.reaction, ["event", "reaction_channel"])
    _, event = load(args.event, ["event", "detector_type", "copy_no", "e", "x", "y", "z"])
    step_tree, step = load(args.step, ["event", "track", "de", "length"]) if args.step else (None, {})

    paths = []
    paths.extend(plot_event_energy(outdir, args.label, event, reaction))
    paths.extend(plot_hit_multiplicity(outdir, args.label, event, reaction))
    paths.extend(plot_positions(outdir, args.label, event))
    paths.extend(plot_step(outdir, args.label, step_tree, step))
    for path in paths:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
