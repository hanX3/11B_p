#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot

DETECTOR_LABELS = {
    2: "LaBr3",
    3: "HPGe",
}

REACTION_CHANNEL_GAMMA_CAPTURE = 3


def read_tree_arrays(root_path, branches):
    with uproot.open(root_path) as f:
        if "tr" not in f:
            raise KeyError(f"ROOT file {root_path} does not contain tree 'tr'")
        tree = f["tr"]
        missing = [b for b in branches if b not in tree.keys()]
        if missing:
            raise KeyError(f"ROOT file {root_path} missing branches: {missing}")
        return tree.arrays(branches, library="np")


def detector_name(detector_type):
    return DETECTOR_LABELS.get(int(detector_type), f"detector_type={int(detector_type)}")


def format_counts_by_key(keys, values):
    lines = []
    for key in sorted(np.unique(keys)):
        count = int(np.sum(keys == key))
        lines.append(f"  {key}: {count}")
    if not lines:
        lines.append("  none")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Plot LaBr3 and HPGe gamma deposited-energy spectra on the same log-scale figure."
    )
    parser.add_argument("--event-root", required=True, help="event_merged_*.root")
    parser.add_argument("--reaction-root", required=True, help="reaction_merged_*.root")
    parser.add_argument("--resonance", type=int, choices=[165, 675], required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--emin", type=float, default=0.0)
    parser.add_argument("--emax", type=float, default=18.0)
    parser.add_argument("--bins", type=int, default=900)
    args = parser.parse_args()

    event_root = Path(args.event_root)
    reaction_root = Path(args.reaction_root)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    reaction_branches = [
        "event",
        "reaction_channel",
        "gamma_resonance",
        "gamma_branch",
        "n_prompt_gammas",
        "gamma1_energy",
        "gamma2_energy",
        "gamma_event_weight",
    ]
    event_branches = [
        "event",
        "detector_type",
        "ring_id",
        "sector",
        "module_id",
        "e",
        "pdg",
        "track_id",
        "parent_id",
    ]

    r = read_tree_arrays(reaction_root, reaction_branches)
    ev = read_tree_arrays(event_root, event_branches)

    selected_reaction_mask = (
        (r["reaction_channel"] == REACTION_CHANNEL_GAMMA_CAPTURE)
        & (r["gamma_resonance"] == args.resonance)
    )
    selected_events = r["event"][selected_reaction_mask].astype(np.int64)

    if selected_events.size == 0:
        raise RuntimeError(
            f"No gamma-capture reaction entries found for gamma_resonance={args.resonance} "
            f"in {reaction_root}"
        )

    detector_mask = (
        np.isin(ev["event"].astype(np.int64), selected_events)
        & np.isin(ev["detector_type"], np.array([2, 3], dtype=ev["detector_type"].dtype))
        & np.isfinite(ev["e"])
        & (ev["e"] > 0.0)
    )

    det_type = ev["detector_type"][detector_mask]
    energy = ev["e"][detector_mask]
    hit_event = ev["event"][detector_mask].astype(np.int64)

    bins = np.linspace(args.emin, args.emax, args.bins + 1)

    plt.figure(figsize=(8.0, 5.2))
    max_count = 0
    plotted_any = False
    for dtype in [2, 3]:
        e_det = energy[det_type == dtype]
        counts, edges = np.histogram(e_det, bins=bins)
        if counts.size > 0:
            max_count = max(max_count, int(counts.max()))
        label = f"{detector_name(dtype)} ({len(e_det)} hits)"
        plt.stairs(counts, edges, label=label, linewidth=1.2)
        plotted_any = plotted_any or len(e_det) > 0

    if not plotted_any:
        raise RuntimeError(
            f"No LaBr3/HPGe detector hits found for gamma_resonance={args.resonance}. "
            "Check /output/saveEvent, detector geometry, and statistics."
        )

    plt.yscale("log")
    plt.ylim(bottom=0.8, top=max(10, max_count * 1.8))
    plt.xlabel("Deposited energy [MeV]")
    plt.ylabel("Counts / bin")
    plt.title(f"{args.resonance}-keV resonance gamma response")
    plt.legend(frameon=False)
    plt.tight_layout()

    spectrum_path = out_prefix.with_name(out_prefix.name + "_labr3_hpge_gamma_spectrum_log.png")
    plt.savefig(spectrum_path, dpi=180)
    plt.close()

    branch_values = r["gamma_branch"][selected_reaction_mask]
    n_prompt_values = r["n_prompt_gammas"][selected_reaction_mask]

    summary_path = out_prefix.with_name(out_prefix.name + "_summary.txt")
    with summary_path.open("w", encoding="utf-8") as f:
        f.write(f"event ROOT       : {event_root}\n")
        f.write(f"reaction ROOT    : {reaction_root}\n")
        f.write(f"resonance        : {args.resonance} keV\n")
        f.write(f"selected reactions: {selected_events.size}\n")
        f.write(f"events with LaBr3/HPGe hit: {len(np.unique(hit_event))}\n")
        f.write(f"total LaBr3/HPGe hits: {len(energy)}\n")
        f.write("\nDetector hit counts:\n")
        for dtype in [2, 3]:
            e_det = energy[det_type == dtype]
            ev_det = hit_event[det_type == dtype]
            f.write(f"  {detector_name(dtype)} hits: {len(e_det)}\n")
            f.write(f"  {detector_name(dtype)} events with hit: {len(np.unique(ev_det))}\n")
            if len(e_det) > 0:
                f.write(f"  {detector_name(dtype)} energy mean [MeV]: {np.mean(e_det):.8g}\n")
                f.write(f"  {detector_name(dtype)} energy max  [MeV]: {np.max(e_det):.8g}\n")
        f.write("\nGamma branch counts in selected reaction tree:\n")
        f.write(format_counts_by_key(branch_values, branch_values))
        f.write("\n\nNumber of prompt gammas per selected reaction:\n")
        f.write(format_counts_by_key(n_prompt_values, n_prompt_values))
        f.write("\n\nEvent-tree detector_type counts after selection:\n")
        f.write(format_counts_by_key(det_type, det_type))
        f.write("\n")

    print(f"Wrote {spectrum_path}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
