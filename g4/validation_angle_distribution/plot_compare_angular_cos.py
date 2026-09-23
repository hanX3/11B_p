#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


def resolve_tree(root_file, requested_tree):
    with uproot.open(root_file) as f:
        keys = {key.split(';')[0]: key for key in f.keys()}
        if requested_tree is not None:
            if requested_tree in keys:
                return requested_tree
            raise KeyError(f"Requested tree '{requested_tree}' not found. Available keys: {list(keys)}")
        if "tr" in keys:
            return "tr"
        if "reaction" in keys:
            return "reaction"
        raise KeyError(f"No supported reaction tree found. Available keys: {list(keys)}")


def load_cos_values(root_file, branch_name, channel, gamma_resonance=None, gamma_branch=None, tree_name=None):
    root_path = Path(root_file)
    if not root_path.exists():
        raise FileNotFoundError(f"ROOT file not found: {root_path}")

    tree_name = resolve_tree(root_path, tree_name)
    with uproot.open(root_path) as f:
        tree = f[tree_name]
        branches = set(tree.keys())

        required = {"reaction_channel", branch_name}
        if gamma_resonance is not None:
            required.add("gamma_resonance")
        if gamma_branch is not None:
            required.add("gamma_branch")

        missing = sorted(required - branches)
        if missing:
            available = "\n".join(sorted(branches))
            raise KeyError(
                f"Missing branch(es) in {root_path}: {missing}\n"
                f"Available branches:\n{available}"
            )

        arrays = tree.arrays(list(required), library="np")

    ch = np.asarray(arrays["reaction_channel"])
    values = np.asarray(arrays[branch_name], dtype=float)
    mask = ch == channel

    if gamma_resonance is not None:
        gres = np.asarray(arrays["gamma_resonance"])
        mask &= gres == gamma_resonance

    if gamma_branch is not None:
        gbr = np.asarray(arrays["gamma_branch"])
        mask &= gbr == gamma_branch

    values = values[mask]
    values = values[np.isfinite(values)]
    values = values[(values >= -1.0) & (values <= 1.0)]
    return values, tree_name


def summarize(label, values):
    print(f"{label}:")
    print(f"  selected entries: {len(values)}")
    if len(values) > 0:
        print(f"  mean: {np.mean(values):.8e}")
        print(f"  std : {np.std(values):.8e}")
        print(f"  min : {np.min(values):.8e}")
        print(f"  max : {np.max(values):.8e}")


def main():
    parser = argparse.ArgumentParser(description="Compare off/on angular cosine distributions from ROOT files.")
    parser.add_argument("--off", required=True, help="ROOT file from angular-off run")
    parser.add_argument("--on", required=True, help="ROOT file from angular-on run")
    parser.add_argument("--channel", required=True, type=int, help="reaction_channel selection")
    parser.add_argument("--branch", required=True, help="cosine branch to plot")
    parser.add_argument("--gamma-resonance", type=int, default=None, help="optional gamma_resonance filter")
    parser.add_argument("--gamma-branch", type=int, default=None, help="optional gamma_branch filter")
    parser.add_argument("--output", required=True, help="output PNG path")
    parser.add_argument("--bins", type=int, default=50, help="histogram bins")
    parser.add_argument("--tree", default=None, help="TTree name. Default: auto-detect tr or reaction")
    args = parser.parse_args()

    off, tree_off = load_cos_values(args.off, args.branch, args.channel, args.gamma_resonance, args.gamma_branch, args.tree)
    on, tree_on = load_cos_values(args.on, args.branch, args.channel, args.gamma_resonance, args.gamma_branch, args.tree)

    print(f"off tree: {tree_off}")
    print(f"on  tree: {tree_on}")
    summarize("off", off)
    summarize("on", on)

    if len(off) == 0:
        print("ERROR: off sample has zero selected entries.", file=sys.stderr)
        return 2
    if len(on) == 0:
        print("ERROR: on sample has zero selected entries.", file=sys.stderr)
        return 2

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    bins = np.linspace(-1.0, 1.0, args.bins + 1)
    plt.figure(figsize=(8, 6))
    plt.hist(off, bins=bins, density=True, histtype="step", linewidth=1.5, label="off")
    plt.hist(on, bins=bins, density=True, histtype="step", linewidth=1.5, label="on")
    plt.xlabel(args.branch)
    plt.ylabel("Normalized counts")
    title = f"channel={args.channel}"
    if args.gamma_resonance is not None:
        title += f", gamma_resonance={args.gamma_resonance}"
    if args.gamma_branch is not None:
        title += f", gamma_branch={args.gamma_branch}"
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output, dpi=200)
    print(f"saved: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
