#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot

FIXED_MARGINS = dict(left=0.18, right=0.86, bottom=0.15, top=0.90)


# Alpha-particle rest mass in MeV/c^2
ALPHA_MASS = 3727.379378


def _get_branch_names(tree):
    """Return branch names without ROOT cycle numbers."""
    return {
        str(name).split(";")[0]
        for name in tree.keys()
    }


def _find_branch_group(available, candidates):
    """
    Return the first complete branch-name group found in the ROOT tree.
    """
    for group in candidates:
        if all(name in available for name in group):
            return group

    return None


def _convert_angle_to_radians(values, name):
    """
    Automatically determine whether an angle is stored in radians or degrees.
    """
    values = np.asarray(values, dtype=float)

    finite = values[np.isfinite(values)]

    if finite.size == 0:
        raise RuntimeError(
            f"No finite {name} values were found."
        )

    maximum = np.max(np.abs(finite))

    # Radian range:
    # theta: 0 to pi
    # phi: -pi to pi or 0 to 2pi
    if maximum <= 2.0 * np.pi + 1.0e-6:
        return values

    # Degree range
    if maximum <= 360.0 + 1.0e-6:
        return np.deg2rad(values)

    raise RuntimeError(
        f"{name} values are outside the expected "
        "radian and degree ranges."
    )


def _boost_momenta_to_three_alpha_cm(
    energy_lab,
    momentum_lab,
):
    """
    Boost the three alpha momenta from LAB to the event-by-event
    three-alpha center-of-mass frame.

    Parameters
    ----------
    energy_lab : ndarray, shape (N, 3)
        Total relativistic energy of each alpha in LAB.

    momentum_lab : ndarray, shape (N, 3, 3)
        LAB momentum vectors. The final index represents px, py, pz.

    Returns
    -------
    ndarray, shape (N, 3, 3)
        Three alpha momentum vectors in the three-alpha CM frame.
    """

    event_total_energy = np.sum(
        energy_lab,
        axis=1,
    )

    event_total_momentum = np.sum(
        momentum_lab,
        axis=1,
    )

    beta = (
        event_total_momentum
        / event_total_energy[:, None]
    )

    beta_squared = np.sum(
        beta**2,
        axis=1,
    )

    if np.any(beta_squared >= 1.0):
        raise RuntimeError(
            "Unphysical CM velocity was obtained: beta^2 >= 1."
        )

    gamma = 1.0 / np.sqrt(
        1.0 - beta_squared
    )

    beta_dot_momentum = np.sum(
        momentum_lab * beta[:, None, :],
        axis=2,
    )

    coefficient = np.zeros_like(
        beta_dot_momentum
    )

    moving = beta_squared > 1.0e-30

    coefficient[moving] = (
        (
            gamma[moving, None] - 1.0
        )
        * beta_dot_momentum[moving]
        / beta_squared[moving, None]
        - gamma[moving, None]
        * energy_lab[moving]
    )

    momentum_cm = (
        momentum_lab
        + coefficient[:, :, None]
        * beta[:, None, :]
    )

    return momentum_cm


def _normalize_vectors(momentum):
    """
    Convert momentum vectors into unit vectors.
    """
    magnitude = np.linalg.norm(
        momentum,
        axis=2,
    )

    valid = np.all(
        magnitude > 0.0,
        axis=1,
    )

    if not np.any(valid):
        raise RuntimeError(
            "No events with nonzero CM alpha momenta were found."
        )

    unit_vectors = (
        momentum[valid]
        / magnitude[valid, :, None]
    )

    return unit_vectors, valid


def plot_675_alpha_pair_angular(
    root_file,
    save=None,
):
    """
    Plot tagged alpha-alpha opening-angle distributions for the
    675-keV sequential decay.

    Particle convention in the labelled ROOT file:

        alpha1 = primary, first-step alpha
        alpha2 = secondary alpha from 8Be decay
        alpha3 = secondary alpha from 8Be decay

    Two distributions are drawn on the same axes:

        1. secondary-secondary:
           alpha2 versus alpha3

        2. primary-secondary:
           alpha1 versus alpha2
           alpha1 versus alpha3

    All angles are calculated in the reconstructed three-alpha CM frame.

    Parameters
    ----------
    root_file : str or pathlib.Path
        Input labelled ROOT file.

    save : str or pathlib.Path or None
        Output image filename. If None, the figure is displayed
        without being saved.

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

    energy_names = [
        "e_alpha1",
        "e_alpha2",
        "e_alpha3",
    ]

    theta_candidates = [
        [
            "theta_lab_alpha1",
            "theta_lab_alpha2",
            "theta_lab_alpha3",
        ],
        [
            "alpha1_theta_lab",
            "alpha2_theta_lab",
            "alpha3_theta_lab",
        ],
    ]

    phi_candidates = [
        [
            "phi_lab_alpha1",
            "phi_lab_alpha2",
            "phi_lab_alpha3",
        ],
        [
            "alpha1_phi_lab",
            "alpha2_phi_lab",
            "alpha3_phi_lab",
        ],
        [
            "phi_alpha1_lab",
            "phi_alpha2_lab",
            "phi_alpha3_lab",
        ],
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

        missing_energy = [
            name
            for name in energy_names
            if name not in available
        ]

        if missing_energy:
            raise RuntimeError(
                "Missing alpha energy branches: "
                + ", ".join(missing_energy)
            )

        theta_names = _find_branch_group(
            available,
            theta_candidates,
        )

        if theta_names is None:
            raise RuntimeError(
                "Could not find the three LAB theta branches.\n"
                "Tried:\n"
                "  theta_lab_alpha1/2/3\n"
                "  alpha1/2/3_theta_lab"
            )

        phi_names = _find_branch_group(
            available,
            phi_candidates,
        )

        if phi_names is None:
            related_phi = sorted(
                name
                for name in available
                if "phi" in name.lower()
            )

            related_text = (
                "\n  ".join(related_phi)
                if related_phi
                else "(none)"
            )

            raise RuntimeError(
                "Could not find the three LAB phi branches.\n"
                "Tried:\n"
                "  phi_lab_alpha1/2/3\n"
                "  alpha1/2/3_phi_lab\n"
                "  phi_alpha1/2/3_lab\n\n"
                "Phi-related branches found:\n"
                f"  {related_text}"
            )

        branches = (
            energy_names
            + theta_names
            + phi_names
        )

        use_branch_id = (
            "branch_id" in available
        )

        if use_branch_id:
            branches.append(
                "branch_id"
            )

        use_resonance_id = (
            "resonance_id" in available
        )

        if use_resonance_id:
            branches.append(
                "resonance_id"
            )

        use_event_weight = (
            "event_weight" in available
        )

        if use_event_weight:
            branches.append(
                "event_weight"
            )

        data = tree.arrays(
            branches,
            library="np",
        )

    kinetic_energy = np.column_stack([
        np.asarray(
            data[name],
            dtype=float,
        )
        for name in energy_names
    ])

    theta_lab = np.column_stack([
        np.asarray(
            data[name],
            dtype=float,
        )
        for name in theta_names
    ])

    phi_lab = np.column_stack([
        np.asarray(
            data[name],
            dtype=float,
        )
        for name in phi_names
    ])

    theta_lab = _convert_angle_to_radians(
        theta_lab,
        "theta",
    )

    phi_lab = _convert_angle_to_radians(
        phi_lab,
        "phi",
    )

    selected = np.ones(
        kinetic_energy.shape[0],
        dtype=bool,
    )

    # Sequential channel
    if use_branch_id:
        branch_id = np.asarray(
            data["branch_id"],
            dtype=int,
        )

        selected &= (
            branch_id == 1
        )

    # 675-keV resonance, when the branch exists
    if use_resonance_id:
        resonance_id = np.asarray(
            data["resonance_id"],
            dtype=int,
        )

        selected &= (
            resonance_id == 675
        )

    valid = (
        selected
        & np.all(
            np.isfinite(kinetic_energy),
            axis=1,
        )
        & np.all(
            np.isfinite(theta_lab),
            axis=1,
        )
        & np.all(
            np.isfinite(phi_lab),
            axis=1,
        )
        & np.all(
            kinetic_energy > 0.0,
            axis=1,
        )
        & np.all(
            theta_lab >= 0.0,
            axis=1,
        )
        & np.all(
            theta_lab <= np.pi,
            axis=1,
        )
    )

    if not np.any(valid):
        raise RuntimeError(
            "No valid labelled 675-keV sequential events were found."
        )

    kinetic_energy = kinetic_energy[valid]
    theta_lab = theta_lab[valid]
    phi_lab = phi_lab[valid]

    if use_event_weight:
        event_weight = np.asarray(
            data["event_weight"],
            dtype=float,
        )[valid]

        bad_weight = (
            ~np.isfinite(event_weight)
            | (event_weight <= 0.0)
        )

        event_weight[bad_weight] = 1.0

    else:
        event_weight = None

    # Total relativistic energy
    energy_lab = (
        kinetic_energy
        + ALPHA_MASS
    )

    # Relativistic momentum magnitude
    momentum_magnitude = np.sqrt(
        kinetic_energy
        * (
            kinetic_energy
            + 2.0 * ALPHA_MASS
        )
    )

    sin_theta = np.sin(
        theta_lab
    )

    px_lab = (
        momentum_magnitude
        * sin_theta
        * np.cos(phi_lab)
    )

    py_lab = (
        momentum_magnitude
        * sin_theta
        * np.sin(phi_lab)
    )

    pz_lab = (
        momentum_magnitude
        * np.cos(theta_lab)
    )

    momentum_lab = np.stack(
        [
            px_lab,
            py_lab,
            pz_lab,
        ],
        axis=2,
    )

    momentum_cm = _boost_momenta_to_three_alpha_cm(
        energy_lab,
        momentum_lab,
    )

    unit_momentum, nonzero_momentum = _normalize_vectors(
        momentum_cm
    )

    if event_weight is not None:
        event_weight = event_weight[
            nonzero_momentum
        ]

    # alpha1 = primary
    primary = unit_momentum[:, 0, :]

    # alpha2 and alpha3 = secondary
    secondary_1 = unit_momentum[:, 1, :]
    secondary_2 = unit_momentum[:, 2, :]

    # Secondary-secondary: one entry per event
    cos_secondary_secondary = np.sum(
        secondary_1 * secondary_2,
        axis=1,
    )

    # Primary-secondary: two entries per event
    cos_primary_secondary_1 = np.sum(
        primary * secondary_1,
        axis=1,
    )

    cos_primary_secondary_2 = np.sum(
        primary * secondary_2,
        axis=1,
    )

    cos_primary_secondary = np.concatenate(
        [
            cos_primary_secondary_1,
            cos_primary_secondary_2,
        ]
    )

    cos_secondary_secondary = np.clip(
        cos_secondary_secondary,
        -1.0,
        1.0,
    )

    cos_primary_secondary = np.clip(
        cos_primary_secondary,
        -1.0,
        1.0,
    )

    if event_weight is not None:
        secondary_secondary_weight = event_weight

        primary_secondary_weight = np.concatenate(
            [
                event_weight,
                event_weight,
            ]
        )

    else:
        secondary_secondary_weight = None
        primary_secondary_weight = None

    bins = np.linspace(
        -1.0,
        1.0,
        41,
    )

    secondary_secondary_hist, _ = np.histogram(
        cos_secondary_secondary,
        bins=bins,
        weights=secondary_secondary_weight,
        density=True,
    )

    primary_secondary_hist, _ = np.histogram(
        cos_primary_secondary,
        bins=bins,
        weights=primary_secondary_weight,
        density=True,
    )

    plt.rcParams.update({
        "font.size": 14,
        "axes.labelsize": 16,
        "axes.titlesize": 16,
        "legend.fontsize": 13,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
        "axes.linewidth": 1.2,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    })

    fig, ax = plt.subplots(
        figsize=(5, 4)
    )

    ax.stairs(
        secondary_secondary_hist,
        bins,
        linewidth=2.2,
        color="#990000",
        label=r"secondary-secondary",
    )

    ax.stairs(
        primary_secondary_hist,
        bins,
        linewidth=2.2,
        color="#003366",
        label=r"primary-secondary",
    )

    ax.set_ylabel(
        "normalised yield"
    )

    ax.set_xlim(
        -1.0,
        1.0,
    )

    ax.set_ylim(
        0.0,
        2.2,
    )

    ax.legend(
        frameon=False,
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
            dpi=200,
        )

        print(f"Saved: {save}")

    number_of_events = (
        cos_secondary_secondary.size
    )

    print(
        f"Selected labelled events: {number_of_events}"
    )

    print(
        "Secondary-secondary entries: "
        f"{cos_secondary_secondary.size}"
    )

    print(
        "Primary-secondary entries: "
        f"{cos_primary_secondary.size}"
    )

    print(
        "Particle convention: "
        "alpha1 = primary, alpha2/alpha3 = secondary"
    )

    plt.show()

    return fig


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot tagged secondary-secondary and primary-secondary "
            "alpha opening-angle distributions in the three-alpha "
            "CM frame for the 675-keV resonance."
        )
    )

    parser.add_argument(
        "root_file",
        type=Path,
        help="Input labelled ROOT file",
    )

    parser.add_argument(
        "--save",
        type=Path,
        default=None,
        help="Output image filename",
    )

    args = parser.parse_args()

    plot_675_alpha_pair_angular(
        args.root_file,
        save=args.save,
    )


if __name__ == "__main__":
    main()
