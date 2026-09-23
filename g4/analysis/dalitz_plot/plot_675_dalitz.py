#!/usr/bin/env python3

import argparse
from itertools import permutations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot

FIXED_MARGINS = dict(left=0.18, right=0.86, bottom=0.15, top=0.90)


DALITZ_LIMIT = 1.08
DALITZ_BINS = 200


def _get_branch_names(tree):
    """Return ROOT branch names without cycle numbers."""
    return {
        str(name).split(";")[0]
        for name in tree.keys()
    }


def _symmetrized_dalitz_coordinates(energies, event_weights=None):
    """
    Calculate symmetrized Dalitz coordinates for three identical alpha
    particles.

    The coordinates are

        x = sqrt(3) * (E2 - E3) / (E1 + E2 + E3)

        y = (2E1 - E2 - E3) / (E1 + E2 + E3)

    All six permutations of E1, E2 and E3 are filled for each event.

    Parameters
    ----------
    energies : ndarray, shape (N, 3)
        Three-alpha CM kinetic energies.

    event_weights : ndarray, shape (N,), optional
        Event weights.

    Returns
    -------
    x : ndarray
        Symmetrized Dalitz x coordinates.

    y : ndarray
        Symmetrized Dalitz y coordinates.

    weights : ndarray or None
        Repeated weights corresponding to the six permutations.
    """

    energies = np.asarray(
        energies,
        dtype=float,
    )

    energy_sum = np.sum(
        energies,
        axis=1,
    )

    valid = (
        np.all(np.isfinite(energies), axis=1)
        & np.all(energies >= 0.0, axis=1)
        & np.isfinite(energy_sum)
        & (energy_sum > 0.0)
    )

    energies = energies[valid]
    energy_sum = energy_sum[valid]

    if event_weights is not None:
        event_weights = np.asarray(
            event_weights,
            dtype=float,
        )[valid]

    if energies.shape[0] == 0:
        raise RuntimeError(
            "No valid three-alpha CM energies were found."
        )

    x_values = []
    y_values = []

    alpha_permutations = list(
        permutations([0, 1, 2])
    )

    for first, second, third in alpha_permutations:
        x = (
            np.sqrt(3.0)
            * (
                energies[:, second]
                - energies[:, third]
            )
            / energy_sum
        )

        y = (
            2.0 * energies[:, first]
            - energies[:, second]
            - energies[:, third]
        ) / energy_sum

        x_values.append(x)
        y_values.append(y)

    x_values = np.concatenate(
        x_values
    )

    y_values = np.concatenate(
        y_values
    )

    if event_weights is not None:
        weights = np.tile(
            event_weights,
            len(alpha_permutations),
        )
    else:
        weights = None

    return (
        x_values,
        y_values,
        weights,
        energies.shape[0],
    )


def plot_675_dalitz(root_file, save=None):
    """
    Plot the fully symmetrized Dalitz distribution for the 675-keV
    sequential three-alpha channel.

    The event selection is:

        branch_id == 1

    If resonance_id exists, the additional selection is:

        resonance_id == 675

    Parameters
    ----------
    root_file : str or pathlib.Path
        Input ROOT file.

    save : str or pathlib.Path or None
        Output filename. If None, the figure is displayed without saving.

    Returns
    -------
    matplotlib.figure.Figure
        Generated figure.
    """

    root_file = Path(root_file)

    if not root_file.is_file():
        raise FileNotFoundError(
            f"ROOT file not found: {root_file}"
        )

    energy_branches = [
        "e_3alpha_cm_alpha1",
        "e_3alpha_cm_alpha2",
        "e_3alpha_cm_alpha3",
    ]

    with uproot.open(root_file) as root:
        object_names = {
            str(name).split(";")[0]
            for name in root.keys()
        }

        if "reaction" not in object_names:
            raise RuntimeError(
                f"'reaction' tree was not found in {root_file}"
            )

        tree = root["reaction"]
        available = _get_branch_names(tree)

        required = (
            ["branch_id"]
            + energy_branches
        )

        missing = [
            name
            for name in required
            if name not in available
        ]

        if missing:
            raise RuntimeError(
                "Missing required branches: "
                + ", ".join(missing)
            )

        branches_to_read = required.copy()

        use_resonance_id = (
            "resonance_id" in available
        )

        if use_resonance_id:
            branches_to_read.append(
                "resonance_id"
            )

        use_event_weight = (
            "event_weight" in available
        )

        if use_event_weight:
            branches_to_read.append(
                "event_weight"
            )

        data = tree.arrays(
            branches_to_read,
            library="np",
        )

    branch_id = np.asarray(
        data["branch_id"],
        dtype=int,
    )

    selected = (
        branch_id == 1
    )

    if use_resonance_id:
        resonance_id = np.asarray(
            data["resonance_id"],
            dtype=int,
        )

        selected &= (
            resonance_id == 675
        )

    energies = np.column_stack([
        np.asarray(
            data[name],
            dtype=float,
        )
        for name in energy_branches
    ])

    energies = energies[selected]

    if use_event_weight:
        event_weights = np.asarray(
            data["event_weight"],
            dtype=float,
        )[selected]

        invalid_weight = (
            ~np.isfinite(event_weights)
            | (event_weights <= 0.0)
        )

        event_weights[invalid_weight] = 1.0

    else:
        event_weights = None

    if energies.shape[0] == 0:
        raise RuntimeError(
            "No selected 675-keV sequential events were found."
        )

    (
        dalitz_x,
        dalitz_y,
        dalitz_weights,
        selected_events,
    ) = _symmetrized_dalitz_coordinates(
        energies,
        event_weights,
    )

    edges = np.linspace(
        -DALITZ_LIMIT,
        DALITZ_LIMIT,
        DALITZ_BINS + 1,
    )

    histogram, x_edges, y_edges = np.histogram2d(
        dalitz_x,
        dalitz_y,
        bins=[edges, edges],
        weights=dalitz_weights,
    )

    # Normalize the complete Dalitz distribution to unit yield.
    total_yield = np.sum(
        histogram
    )

    if total_yield > 0.0:
        histogram = (
            histogram / total_yield
        )

    x_centers = 0.5 * (
        x_edges[:-1]
        + x_edges[1:]
    )

    y_centers = 0.5 * (
        y_edges[:-1]
        + y_edges[1:]
    )

    xx, yy = np.meshgrid(
        x_centers,
        y_centers,
        indexing="ij",
    )

    outside_boundary = (
        xx**2 + yy**2 > 1.0
    )

    empty_bins = (
        histogram <= 0.0
    )

    masked_histogram = np.ma.masked_where(
        outside_boundary | empty_bins,
        histogram,
    )

    cmap = plt.get_cmap(
        "viridis"
    ).copy()

    cmap.set_bad(
        "white"
    )

    plt.rcParams.update({
        "font.size": 14,
        "axes.labelsize": 16,
        "axes.titlesize": 16,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
        "axes.linewidth": 1.2,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    })

    fig, ax = plt.subplots(
        figsize=(16/3, 4)
    )

    ax.pcolormesh(
        x_edges,
        y_edges,
        masked_histogram.T,
        cmap=cmap,
        shading="flat",
        rasterized=True,
    )

    boundary_angle = np.linspace(
        0.0,
        2.0 * np.pi,
        600,
    )

    ax.plot(
        np.cos(boundary_angle),
        np.sin(boundary_angle),
        linestyle="--",
        linewidth=1.2,
        color="#777777",
    )

    ax.set_xlabel(
        r"$x=\sqrt{3}\,(E_2-E_3)/\sum_i E_i$"
    )

    ax.set_ylabel(
        r"$y=(2E_1-E_2-E_3)/\sum_i E_i$"
    )

    ax.set_xlim(
        -DALITZ_LIMIT,
        DALITZ_LIMIT,
    )

    ax.set_ylim(
        -DALITZ_LIMIT,
        DALITZ_LIMIT,
    )

    ax.set_aspect(
        "equal"
    )

    fig.subplots_adjust(**FIXED_MARGINS)

    if save is not None:
        save = Path(save)

        save.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        fig.savefig(
            save,
            dpi=250,
        )

        print(f"Saved: {save}")

    print(
        f"Selected events: {selected_events}"
    )

    print(
        f"Filled Dalitz entries: {dalitz_x.size}"
    )

    print(
        "Permutations per event: 6"
    )

    plt.show()

    return fig


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot the fully symmetrized Dalitz distribution "
            "for the 675-keV sequential three-alpha channel."
        )
    )

    parser.add_argument(
        "root_file",
        type=Path,
        help="Input ROOT file",
    )

    parser.add_argument(
        "--save",
        type=Path,
        default=None,
        help="Output image filename",
    )

    args = parser.parse_args()

    plot_675_dalitz(
        args.root_file,
        save=args.save,
    )


if __name__ == "__main__":
    main()
