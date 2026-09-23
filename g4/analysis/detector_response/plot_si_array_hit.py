#!/usr/bin/env python3
"""
Plot Si hit positions from the event tree.

Default behaviour:
  - alpha particles only;
  - use entry position;
  - keep axes visible;
  - display rotated so the beam direction runs from screen left to right.

Input:
  ROOT file containing the event tree with schema version 6 branches:
    si_hit_pdg[]
    si_hit_edep_MeV[]
    si_hit_x_entry_mm[], si_hit_y_entry_mm[], si_hit_z_entry_mm[]
    si_hit_x_edep_mm[],  si_hit_y_edep_mm[],  si_hit_z_edep_mm[]

Dual use:
    terminal :  python3 plot_si_array_hit.py <file.root>
                python3 plot_si_array_hit.py <file.root> --output my_name.png
                (saves si_array_hit.png)
    jupyter  :  from plot_si_array_hit import plot_si_array_hit
                fig = plot_si_array_hit("file.root")
                fig = plot_si_array_hit("file.root", save="my_name.png",
                                        position="edep", all_particles=True)
                (figure displays inline automatically under %matplotlib inline)

Dependencies:
  python3 -m pip install uproot awkward numpy matplotlib
"""

from __future__ import annotations

from pathlib import Path

import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
import uproot

ALPHA_PDG = 1000020040
TARGET_Z_MM = 170.0
DEFAULT_OUTPUT = "si_array_hit.png"


def set_equal_axes(ax, x, y, z):
    span = max(np.ptp(x), np.ptp(y), np.ptp(z), 1.0)
    half = 0.55 * span
    xmid = 0.5 * (np.min(x) + np.max(x))
    ymid = 0.5 * (np.min(y) + np.max(y))
    zmid = 0.5 * (np.min(z) + np.max(z))
    ax.set_xlim(xmid - half, xmid + half)
    ax.set_ylim(ymid - half, ymid + half)
    ax.set_zlim(zmid - half, zmid + half)


def plot_si_array_hit(root_file, save=False, position="entry",
                      all_particles=False, transparent=False):
    """Plot Si hit positions from the event tree.

    Parameters
    ----------
    root_file : str or pathlib.Path
        ROOT file containing the ``event`` tree (schema version 6).
    save : bool or str, optional
        False (default) -> do not write a file.
        True            -> write si_array_hit.png.
        str             -> write to that path (".png" appended if missing).
    position : {"entry", "edep"}, optional
        Use the entry position or the energy-weighted position.
    all_particles : bool, optional
        Plot all particles instead of alpha particles only.
    transparent : bool, optional
        Save with a transparent background.

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    if position not in ("entry", "edep"):
        raise ValueError("position must be 'entry' or 'edep'")

    root_path = Path(root_file).expanduser().resolve()

    x_branch = f"si_hit_x_{position}_mm"
    y_branch = f"si_hit_y_{position}_mm"
    z_branch = f"si_hit_z_{position}_mm"
    required = ["si_hit_pdg", "si_hit_edep_MeV",
                x_branch, y_branch, z_branch]

    with uproot.open(root_path) as root:
        if "event" not in root:
            raise ValueError("ROOT file does not contain event.")
        tree = root["event"]
        missing = [name for name in required if name not in tree.keys()]
        if missing:
            raise ValueError(
                f"Missing branches: {missing}. "
                "Recompile with schema version 6 and rerun the simulation.")
        arrays = tree.arrays(required, library="ak")

    pdg = ak.to_numpy(ak.flatten(arrays["si_hit_pdg"]))
    energy = ak.to_numpy(ak.flatten(arrays["si_hit_edep_MeV"]))
    x = ak.to_numpy(ak.flatten(arrays[x_branch]))
    y = ak.to_numpy(ak.flatten(arrays[y_branch]))
    z = ak.to_numpy(ak.flatten(arrays[z_branch]))

    mask = (
        np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
        & np.isfinite(energy) & (energy > 0.0)
    )
    if not all_particles:
        mask &= np.abs(pdg) == ALPHA_PDG

    # Coordinates relative to target centre (0, 0, 170 mm).
    x = x[mask]
    y = y[mask]
    z = z[mask] - TARGET_Z_MM

    if len(x) == 0:
        raise ValueError("No matching Si hits found.")

    # Rotate display coordinates so the beam direction (+z) runs left -> right.
    x_disp = z
    y_disp = x
    z_disp = y

    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(
        x_disp, y_disp, z_disp,
        s=1.2, alpha=0.42, depthshade=False, linewidths=0,
    )

    set_equal_axes(ax, x_disp, y_disp, z_disp)
    ax.view_init(elev=12.0, azim=-90.0)
    ax.grid(False)

    ax.set_xlabel(r"$z-z_{\rm target}$ [mm]  (beam direction)")
    ax.set_ylabel(r"$x-x_{\rm target}$ [mm]")
    ax.set_zlabel(r"$y-y_{\rm target}$ [mm]")

    fig.tight_layout()

    print(f"plotted hits: {len(x_disp)}")
    print(f"position: {position}")

    if save:
        out = DEFAULT_OUTPUT if save is True else str(save)
        if not out.lower().endswith(".png"):
            out += ".png"
        out_path = Path(out).expanduser().resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=300, bbox_inches="tight",
                    transparent=transparent)
        print(f"saved: {out_path}")

    return fig


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Plot Si hit positions from the event tree.")
    parser.add_argument("root_file", type=Path)
    parser.add_argument("--output", type=Path, default=None,
                        help="output png (default: si_array_hit.png)")
    parser.add_argument("--position", choices=("entry", "edep"),
                        default="entry",
                        help="entry or energy-weighted position")
    parser.add_argument("--all-particles", action="store_true",
                        help="plot all particles, not just alphas")
    parser.add_argument("--transparent", action="store_true",
                        help="save with transparent background")
    args = parser.parse_args()

    out = True if args.output is None else str(args.output)
    plot_si_array_hit(args.root_file, save=out, position=args.position,
                      all_particles=args.all_particles,
                      transparent=args.transparent)
