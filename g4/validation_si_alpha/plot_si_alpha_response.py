#!/usr/bin/env python3
"""Plot Si alpha response for 11B(p,3alpha) resonance simulations.

Inputs:
  - reaction ROOT tree: generator-level truth, one row per reaction event.
  - event ROOT tree: detector hits, one row per detector channel hit.

The script joins the two trees by event id and produces:
  - generated alpha lab spectra, split by branch and primary/cascade labels;
  - Si deposited-energy spectra, split by branch and approximate alpha track labels;
  - sector/channel response maps;
  - summary text files.

Assumptions for detected alpha labels:
  Geant4 usually assigns the three alpha secondaries track_id 2, 3, 4 in the
  order they are added in H11BReaction.cc: alpha1(primary), alpha2, alpha3.
  Therefore track_id==2 is treated as the primary alpha and track_id==3/4
  as the two cascade/secondary alphas. The script also writes a track-id
  diagnostic so this can be checked from the output.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot

ALPHA_PDG = 1000020040
SI_DETECTOR_TYPE = 1

REACTION_REQUIRED = [
    "event",
    "reaction_channel",
    "resonance_id",
    "branch_id",
    "e_alpha1",
    "e_alpha2",
    "e_alpha3",
]

EVENT_REQUIRED = [
    "event",
    "detector_type",
    "ring_id",
    "sector",
    "e",
    "pdg",
    "track_id",
    "parent_id",
]


@dataclass(frozen=True)
class Component:
    key: str
    label: str
    query: str


def open_tree(path: Path, tree_name: str | None = None):
    with uproot.open(path) as f:
        names = {k.split(";")[0] for k in f.keys()}
        if tree_name is not None:
            if tree_name not in names:
                raise RuntimeError(f"Tree {tree_name!r} not found in {path}. Available: {sorted(names)}")
            return tree_name
        if "tr" in names:
            return "tr"
        if "reaction" in names:
            return "reaction"
        raise RuntimeError(f"No tr/reaction tree found in {path}. Available: {sorted(names)}")


def load_reaction(path: Path, tree_name: str | None, channel: int) -> pd.DataFrame:
    actual_tree = open_tree(path, tree_name)
    with uproot.open(path) as f:
        tree = f[actual_tree]
        available = set(tree.keys())
        missing = [b for b in REACTION_REQUIRED if b not in available]
        if missing:
            raise RuntimeError(f"Reaction tree missing branches: {missing}")
        branches = list(REACTION_REQUIRED)
        for optional in [
            "event_sampling_weight",
            "event_weight",
            "cos_theta_primary_cm",
            "cos_theta_secondary_correlation",
            "primary_angular_mode",
            "secondary_angular_model",
            "primary_a1_used",
            "primary_a2_used",
            "secondary_a2_used",
            "secondary_a4_used",
        ]:
            if optional in available:
                branches.append(optional)
        arr = pd.DataFrame(tree.arrays(branches, library="np"))

    df = arr[arr["reaction_channel"] == channel].copy()
    if df.empty:
        raise RuntimeError(
            f"No reaction entries with reaction_channel=={channel} in {path}. "
            "Check that you passed the correct reaction ROOT file."
        )
    if "event_sampling_weight" not in df.columns:
        df["event_sampling_weight"] = 1.0
    if "event_weight" not in df.columns:
        df["event_weight"] = 1.0
    return df


def load_event(path: Path, tree_name: str | None) -> pd.DataFrame:
    actual_tree = open_tree(path, tree_name)
    with uproot.open(path) as f:
        tree = f[actual_tree]
        available = set(tree.keys())
        missing = [b for b in EVENT_REQUIRED if b not in available]
        if missing:
            raise RuntimeError(f"Event tree missing branches: {missing}")
        branches = list(EVENT_REQUIRED)
        for optional in ["array_id", "module_id", "segment_id", "copy_no", "detector"]:
            if optional in available:
                branches.append(optional)
        arr = pd.DataFrame(tree.arrays(branches, library="np"))
    return arr


def branch_name(branch_id: int) -> str:
    if branch_id == 0:
        return "branch0_alpha0_to_8Be_gs"
    if branch_id == 1:
        return "branch1_alpha1_to_8Be_exc"
    return f"branch{branch_id}"


def generated_component_frame(reaction: pd.DataFrame, use_weights: bool) -> pd.DataFrame:
    rows = []
    for _, r in reaction.iterrows():
        weight = float(r["event_sampling_weight"] if use_weights else 1.0)
        branch = int(r["branch_id"])
        if branch == 0:
            primary_label = "branch0 primary alpha0"
            secondary_label = "branch0 cascade alpha"
        elif branch == 1:
            primary_label = "branch1 primary alpha1"
            secondary_label = "branch1 cascade alpha"
        else:
            primary_label = f"branch{branch} primary alpha"
            secondary_label = f"branch{branch} cascade alpha"
        rows.append((float(r["e_alpha1"]), branch, branch_name(branch), primary_label, weight))
        rows.append((float(r["e_alpha2"]), branch, branch_name(branch), secondary_label, weight))
        rows.append((float(r["e_alpha3"]), branch, branch_name(branch), secondary_label, weight))
    return pd.DataFrame(rows, columns=["energy", "branch_id", "branch_name", "component", "weight"])


def label_detected_alpha(row: pd.Series) -> str:
    # Track ID convention from H11BReaction.cc secondary insertion order.
    # track_id==2: alpha1 primary; track_id==3/4: two secondary/cascade alphas.
    track_id = int(row["track_id"])
    branch = int(row["branch_id"])
    if branch == 0:
        if track_id == 2:
            return "branch0 primary alpha0"
        if track_id in (3, 4):
            return "branch0 cascade alpha"
        return "branch0 other alpha track"
    if branch == 1:
        if track_id == 2:
            return "branch1 primary alpha1"
        if track_id in (3, 4):
            return "branch1 cascade alpha"
        return "branch1 other alpha track"
    return f"branch{branch} alpha track {track_id}"


def select_si_alpha_hits(event: pd.DataFrame, reaction: pd.DataFrame, ring: int | None) -> pd.DataFrame:
    hit = event[
        (event["detector_type"] == SI_DETECTOR_TYPE)
        & (event["pdg"] == ALPHA_PDG)
        & np.isfinite(event["e"])
    ].copy()
    if ring is not None:
        hit = hit[hit["ring_id"] == ring].copy()

    event_to_reaction = reaction[
        ["event", "reaction_channel", "resonance_id", "branch_id", "event_sampling_weight", "event_weight"]
    ].drop_duplicates("event")
    joined = hit.merge(event_to_reaction, on="event", how="inner")
    if joined.empty:
        raise RuntimeError(
            "No Si alpha hits matched to selected reaction events. "
            "Check event ROOT file, ring selection, saveEvent=true, and detector geometry."
        )
    joined["branch_name"] = joined["branch_id"].astype(int).map(branch_name)
    joined["component"] = joined.apply(label_detected_alpha, axis=1)
    return joined


def hist_step(data: np.ndarray, bins: np.ndarray, weights: np.ndarray | None, label: str):
    if len(data) == 0:
        return
    plt.hist(data, bins=bins, weights=weights, histtype="step", label=label)


def save_hist_components(
    df: pd.DataFrame,
    value_col: str,
    group_col: str,
    weight_col: str,
    bins: np.ndarray,
    title: str,
    xlabel: str,
    output: Path,
    density: bool = False,
):
    plt.figure(figsize=(8, 5))
    for label, sub in sorted(df.groupby(group_col), key=lambda x: str(x[0])):
        weights = sub[weight_col].to_numpy() if weight_col in sub.columns else None
        if density and weights is not None:
            # Normalize each component to unit area for shape comparison.
            wsum = np.sum(weights)
            if wsum > 0:
                weights = weights / wsum
        elif density:
            weights = np.ones(len(sub)) / max(len(sub), 1)
        hist_step(sub[value_col].to_numpy(), bins, weights, str(label))
    plt.xlabel(xlabel)
    plt.ylabel("normalized counts" if density else "counts")
    plt.title(title)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(output, dpi=170)
    plt.close()


def save_total_hist(
    df: pd.DataFrame,
    value_col: str,
    weight_col: str,
    bins: np.ndarray,
    title: str,
    xlabel: str,
    output: Path,
):
    weights = df[weight_col].to_numpy() if weight_col in df.columns else None
    plt.figure(figsize=(8, 5))
    plt.hist(df[value_col].to_numpy(), bins=bins, weights=weights, histtype="step")
    plt.xlabel(xlabel)
    plt.ylabel("counts")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output, dpi=170)
    plt.close()


def save_sector_map(hit: pd.DataFrame, bins_e: np.ndarray, output: Path, title: str):
    plt.figure(figsize=(8, 5))
    sectors = hit["sector"].to_numpy()
    energy = hit["e"].to_numpy()
    if len(sectors) == 0:
        return
    smin = int(np.nanmin(sectors))
    smax = int(np.nanmax(sectors))
    sector_bins = np.arange(smin - 0.5, smax + 1.5, 1.0)
    plt.hist2d(sectors, energy, bins=[sector_bins, bins_e])
    plt.xlabel("Si sector")
    plt.ylabel("Si deposited energy [MeV]")
    plt.title(title)
    plt.colorbar(label="counts")
    plt.tight_layout()
    plt.savefig(output, dpi=170)
    plt.close()


def save_multiplicity(hit: pd.DataFrame, reaction: pd.DataFrame, output: Path, title: str):
    # Count number of Si alpha hits per selected reaction event.
    counts = hit.groupby("event").size().rename("n_si_alpha_hits").reset_index()
    base = reaction[["event", "branch_id"]].drop_duplicates("event").merge(counts, on="event", how="left")
    base["n_si_alpha_hits"] = base["n_si_alpha_hits"].fillna(0).astype(int)
    plt.figure(figsize=(8, 5))
    bins = np.arange(-0.5, max(6, base["n_si_alpha_hits"].max() + 1.5), 1.0)
    for bid, sub in sorted(base.groupby("branch_id")):
        plt.hist(sub["n_si_alpha_hits"], bins=bins, histtype="step", label=branch_name(int(bid)))
    plt.xlabel("Number of Si alpha hits per reaction event")
    plt.ylabel("events")
    plt.title(title)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(output, dpi=170)
    plt.close()


def write_summary(
    reaction: pd.DataFrame,
    hit: pd.DataFrame,
    gen: pd.DataFrame,
    out: Path,
    label: str,
    ring: int | None,
    weighted: bool,
):
    with open(out, "w") as w:
        w.write(f"Si alpha response summary: {label}\n")
        w.write("========================================\n\n")
        w.write(f"Selected Si ring: {'all Si rings' if ring is None else ring}\n")
        w.write(f"Weighted spectra: {weighted}\n\n")
        w.write("Reaction entries by channel / branch:\n")
        w.write(reaction.groupby(["reaction_channel", "resonance_id", "branch_id"]).size().to_string())
        w.write("\n\n")
        w.write("Generated alpha entries by component:\n")
        w.write(gen.groupby(["branch_id", "component"]).size().to_string())
        w.write("\n\n")
        w.write("Matched Si alpha hits by branch / component:\n")
        w.write(hit.groupby(["branch_id", "component"]).size().to_string())
        w.write("\n\n")
        w.write("Matched Si alpha hits by ring / sector:\n")
        w.write(hit.groupby(["ring_id", "sector"]).size().to_string())
        w.write("\n\n")
        w.write("Alpha hit track_id counts:\n")
        w.write(hit.groupby(["track_id", "branch_id", "component"]).size().to_string())
        w.write("\n\n")
        event_hit_counts = hit.groupby("event").size()
        n_reaction_events = reaction["event"].nunique()
        n_hit_events = event_hit_counts.index.nunique()
        w.write(f"Selected reaction events: {n_reaction_events}\n")
        w.write(f"Events with >=1 selected Si alpha hit: {n_hit_events}\n")
        w.write(f"Event-level selected Si alpha efficiency: {n_hit_events / max(n_reaction_events, 1):.6g}\n")
        w.write(f"Total selected Si alpha hits: {len(hit)}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reaction-root", required=True, type=Path)
    ap.add_argument("--event-root", required=True, type=Path)
    ap.add_argument("--channel", required=True, type=int, choices=[0, 1], help="0=165 resonance, 1=675 resonance")
    ap.add_argument("--label", required=True, help="label used in plot titles and filenames")
    ap.add_argument("--ring", default="all", help="Si ring_id to select; use 'all' for all Si rings. Default: all")
    ap.add_argument("--emin", type=float, default=0.0)
    ap.add_argument("--emax", type=float, default=8.0)
    ap.add_argument("--bins", type=int, default=400)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--weighted", action="store_true", help="Use event_sampling_weight for spectra")
    ap.add_argument("--reaction-tree", default=None)
    ap.add_argument("--event-tree", default=None)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    ring = None if str(args.ring).lower() == "all" else int(args.ring)

    reaction = load_reaction(args.reaction_root, args.reaction_tree, args.channel)
    event = load_event(args.event_root, args.event_tree)
    hit = select_si_alpha_hits(event, reaction, ring)

    gen = generated_component_frame(reaction, args.weighted)
    gen_weight_col = "weight"
    hit["weight"] = hit["event_sampling_weight"] if args.weighted else 1.0

    bins_e = np.linspace(args.emin, args.emax, args.bins + 1)

    save_total_hist(
        gen,
        "energy",
        gen_weight_col,
        bins_e,
        f"Generated lab alpha spectrum: {args.label}",
        "Generated alpha kinetic energy in lab [MeV]",
        args.outdir / f"{args.label}_generated_total.png",
    )
    save_hist_components(
        gen,
        "energy",
        "branch_name",
        gen_weight_col,
        bins_e,
        f"Generated lab alpha spectrum by branch: {args.label}",
        "Generated alpha kinetic energy in lab [MeV]",
        args.outdir / f"{args.label}_generated_by_branch_counts.png",
    )
    save_hist_components(
        gen,
        "energy",
        "component",
        gen_weight_col,
        bins_e,
        f"Generated lab alpha spectrum by component: {args.label}",
        "Generated alpha kinetic energy in lab [MeV]",
        args.outdir / f"{args.label}_generated_by_component_counts.png",
    )
    save_hist_components(
        gen,
        "energy",
        "component",
        gen_weight_col,
        bins_e,
        f"Generated lab alpha spectrum shape by component: {args.label}",
        "Generated alpha kinetic energy in lab [MeV]",
        args.outdir / f"{args.label}_generated_by_component_shape.png",
        density=True,
    )

    save_total_hist(
        hit,
        "e",
        "weight",
        bins_e,
        f"Si alpha deposited-energy spectrum: {args.label}",
        "Si deposited energy [MeV]",
        args.outdir / f"{args.label}_si_total.png",
    )
    save_hist_components(
        hit,
        "e",
        "branch_name",
        "weight",
        bins_e,
        f"Si alpha deposited energy by branch: {args.label}",
        "Si deposited energy [MeV]",
        args.outdir / f"{args.label}_si_by_branch_counts.png",
    )
    save_hist_components(
        hit,
        "e",
        "component",
        "weight",
        bins_e,
        f"Si alpha deposited energy by component: {args.label}",
        "Si deposited energy [MeV]",
        args.outdir / f"{args.label}_si_by_component_counts.png",
    )
    save_hist_components(
        hit,
        "e",
        "component",
        "weight",
        bins_e,
        f"Si alpha deposited-energy shape by component: {args.label}",
        "Si deposited energy [MeV]",
        args.outdir / f"{args.label}_si_by_component_shape.png",
        density=True,
    )
    save_sector_map(
        hit,
        bins_e,
        args.outdir / f"{args.label}_si_sector_energy_map.png",
        f"Si sector vs alpha deposited energy: {args.label}",
    )
    save_multiplicity(
        hit,
        reaction,
        args.outdir / f"{args.label}_si_alpha_hit_multiplicity.png",
        f"Si alpha hit multiplicity: {args.label}",
    )

    write_summary(
        reaction,
        hit,
        gen,
        args.outdir / f"{args.label}_summary.txt",
        args.label,
        ring,
        args.weighted,
    )

    print("Wrote plots and summary to:", args.outdir)
    print("Key files:")
    for p in [
        args.outdir / f"{args.label}_summary.txt",
        args.outdir / f"{args.label}_si_total.png",
        args.outdir / f"{args.label}_si_by_component_counts.png",
        args.outdir / f"{args.label}_si_sector_energy_map.png",
        args.outdir / f"{args.label}_generated_by_component_counts.png",
    ]:
        print("  ", p)


if __name__ == "__main__":
    main()
