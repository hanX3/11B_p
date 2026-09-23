#!/usr/bin/env python3
"""
Plot alpha-particle crossings on the virtual sphere.

Fixed settings:
    - alpha particles only;
    - virtual-sphere centre at (0, 0, 170 mm);
    - coordinate axes retained;
    - no figure title;
    - beam direction displayed from screen left to right.

Dual use:
    terminal :  python3 plot_virtual_sphere.py <merged.root>
                (saves virtual_sphere.png)
    jupyter  :  from plot_virtual_sphere import plot_virtual_sphere
                fig = plot_virtual_sphere("merged.root")
                fig = plot_virtual_sphere("merged.root", save="my_name.png")
                (figure displays inline automatically under %matplotlib inline)
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


ALPHA_PDG = 1000020040
TARGET_Z_MM = 170.0
MAX_POINTS = 40000
DEFAULT_OUTPUT = "virtual_sphere.png"


def downsample_indices(number: int, limit: int) -> np.ndarray:
    if number <= limit:
        return np.arange(number)

    rng = np.random.default_rng(12345)
    return np.sort(
        rng.choice(number, size=limit, replace=False)
    )


def set_equal_axes(axis, x, y, z):
    span = max(np.ptp(x), np.ptp(y), np.ptp(z), 1.0)
    half = 0.6 * span

    x_mid = 0.5 * (np.min(x) + np.max(x))
    y_mid = 0.5 * (np.min(y) + np.max(y))
    z_mid = 0.5 * (np.min(z) + np.max(z))

    axis.set_xlim(x_mid - half, x_mid + half)
    axis.set_ylim(y_mid - half, y_mid + half)
    axis.set_zlim(z_mid - half, z_mid + half)


def plot_virtual_sphere(root_path, save=False, max_points=MAX_POINTS):
    """Plot alpha-particle crossings on the virtual sphere.

    Parameters
    ----------
    root_path : str or pathlib.Path
        ROOT file containing the ``virtual_sphere`` tree.
    save : bool or str, optional
        False (default) -> do not write a file.
        True            -> write virtual_sphere.png.
        str             -> write to that path (".png" appended if missing).
    max_points : int, optional
        Cap on plotted points (random downsample above this; default 40000).

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    root_path = Path(root_path).expanduser().resolve()

    with uproot.open(root_path) as root_file:
        if "virtual_sphere" not in root_file:
            raise ValueError(
                "ROOT file does not contain the virtual_sphere tree.")

        tree = root_file["virtual_sphere"]
        required_branches = ["pdg", "x_mm", "y_mm", "z_mm"]
        missing = [b for b in required_branches if b not in tree.keys()]
        if missing:
            raise ValueError(f"Missing branches: {missing}")

        data = tree.arrays(required_branches, library="np")

    alpha_mask = np.abs(data["pdg"]) == ALPHA_PDG
    x = data["x_mm"][alpha_mask]
    y = data["y_mm"][alpha_mask]
    z = data["z_mm"][alpha_mask] - TARGET_Z_MM

    if len(x) == 0:
        raise ValueError(
            "No alpha-particle virtual-sphere crossings were found.")

    selected = downsample_indices(len(x), max_points)
    x = x[selected]
    y = y[selected]
    z = z[selected]

    radius = np.median(np.sqrt(x * x + y * y + z * z))

    # Display-coordinate mapping:
    # physical beam direction +z -> horizontal direction in the figure.
    display_x = z
    display_y = x
    display_z = y

    figure = plt.figure(figsize=(7, 7))
    axis = figure.add_subplot(111, projection="3d")

    axis.scatter(
        display_x, display_y, display_z,
        s=0.65, alpha=0.38, depthshade=False, linewidths=0,
    )

    # Thin wireframe sphere for geometric reference.
    u = np.linspace(0.0, 2.0 * np.pi, 48)
    v = np.linspace(0.0, np.pi, 24)
    sphere_x = radius * np.outer(np.ones_like(u), np.cos(v))
    sphere_y = radius * np.outer(np.cos(u), np.sin(v))
    sphere_z = radius * np.outer(np.sin(u), np.sin(v))

    axis.plot_wireframe(
        sphere_x, sphere_y, sphere_z,
        rstride=4, cstride=4, linewidth=0.30, alpha=0.15,
    )

    set_equal_axes(axis, display_x, display_y, display_z)

    axis.view_init(elev=10.0, azim=-90.0)
    axis.grid(False)

    axis.set_xlabel(
        r"$z-z_{\rm target}$ [mm]  (beam direction)", labelpad=10)
    axis.set_ylabel(r"$x-x_{\rm target}$ [mm]", labelpad=10)
    axis.set_zlabel(r"$y-y_{\rm target}$ [mm]", labelpad=10)

    figure.tight_layout()

    print(f"alpha-particle crossings: {len(x)}")
    print(f"median radius: {radius:.6g} mm")

    if save:
        out = DEFAULT_OUTPUT if save is True else str(save)
        if not out.lower().endswith(".png"):
            out += ".png"
        figure.savefig(out, dpi=300, bbox_inches="tight")
        print(f"saved: {Path(out).resolve()}")

    return figure


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(
        description="Plot alpha-particle crossings on the virtual sphere.")
    ap.add_argument("rootfile", help="merged ROOT file (virtual_sphere tree)")
    ap.add_argument("--out", default=None,
                    help="output png name (default: virtual_sphere.png)")
    ap.add_argument("--max-points", type=int, default=MAX_POINTS,
                    help="cap on plotted points (default 40000)")
    args = ap.parse_args()

    out = True if args.out is None else args.out
    plot_virtual_sphere(args.rootfile, save=out, max_points=args.max_points)
