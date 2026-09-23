#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot

FIXED_MARGINS = dict(left=0.18, right=0.86, bottom=0.15, top=0.90)


# alpha mass in MeV/c^2
ALPHA_MASS = 3727.379378


def _convert_theta_to_radians(theta):
    theta = np.asarray(theta, dtype=float)

    finite = theta[np.isfinite(theta)]
    if finite.size == 0:
        raise RuntimeError("No finite alpha angles were found.")

    maximum = np.max(np.abs(finite))

    if maximum <= np.pi + 1.0e-6:
        return theta

    if maximum <= 180.0 + 1.0e-6:
        return np.deg2rad(theta)

    raise RuntimeError(
        "The alpha angle values are outside both the radian and degree ranges."
    )


def plot_675_alpha_angular(root_file, save=None):
    """
    Plot LAB and three-alpha CM angular distributions on the same axis.

    All three alpha particles are filled into each distribution.
    The selection is fixed to branch_id == 1.

    Parameters
    ----------
    root_file : str or pathlib.Path
        Input ROOT file.

    save : str or pathlib.Path or None
        Output image filename. If None, the figure is only displayed.

    Returns
    -------
    matplotlib.figure.Figure
        Generated figure.
    """

    root_file = Path(root_file)

    if not root_file.exists():
        raise FileNotFoundError(f"ROOT file not found: {root_file}")

    branches = [
        "branch_id",
        "e_alpha1",
        "e_alpha2",
        "e_alpha3",
        "theta_lab_alpha1",
        "theta_lab_alpha2",
        "theta_lab_alpha3",
    ]

    with uproot.open(root_file) as root:
        object_names = {str(name).split(";")[0] for name in root.keys()}

        if "reaction" not in object_names:
            raise RuntimeError(f"'reaction' tree was not found in {root_file}")

        tree = root["reaction"]
        available = {str(name).split(";")[0] for name in tree.keys()}

        missing = [name for name in branches if name not in available]
        if missing:
            raise RuntimeError("Missing required branches: " + ", ".join(missing))

        data = tree.arrays(branches, library="np")

    branch_id = np.asarray(data["branch_id"], dtype=int)

    kinetic_energy = np.column_stack([
        np.asarray(data["e_alpha1"], dtype=float),
        np.asarray(data["e_alpha2"], dtype=float),
        np.asarray(data["e_alpha3"], dtype=float),
    ])

    theta_lab = np.column_stack([
        np.asarray(data["theta_lab_alpha1"], dtype=float),
        np.asarray(data["theta_lab_alpha2"], dtype=float),
        np.asarray(data["theta_lab_alpha3"], dtype=float),
    ])

    theta_lab = _convert_theta_to_radians(theta_lab)

    selected = (branch_id == 1)

    valid = (
        selected
        & np.all(np.isfinite(kinetic_energy), axis=1)
        & np.all(np.isfinite(theta_lab), axis=1)
        & np.all(kinetic_energy >= 0.0, axis=1)
        & np.all(theta_lab >= 0.0, axis=1)
        & np.all(theta_lab <= np.pi, axis=1)
    )

    if not np.any(valid):
        raise RuntimeError("No valid branch_id == 1 three-alpha events were found.")

    kinetic_energy = kinetic_energy[valid]
    theta_lab = theta_lab[valid]

    total_energy = kinetic_energy + ALPHA_MASS
    momentum = np.sqrt(kinetic_energy * (kinetic_energy + 2.0 * ALPHA_MASS))

    pt_lab = momentum * np.sin(theta_lab)
    pz_lab = momentum * np.cos(theta_lab)

    # LAB distribution
    cos_theta_lab = np.cos(theta_lab).reshape(-1)

    # three-alpha CM boost
    event_total_energy = np.sum(total_energy, axis=1)
    event_total_pz = np.sum(pz_lab, axis=1)

    beta_cm = event_total_pz / event_total_energy

    if np.any(np.abs(beta_cm) >= 1.0):
        raise RuntimeError("Unphysical CM velocity was obtained: |beta_CM| >= 1.")

    gamma_cm = 1.0 / np.sqrt(1.0 - beta_cm**2)

    pz_cm = gamma_cm[:, None] * (pz_lab - beta_cm[:, None] * total_energy)
    p_cm = np.sqrt(pt_lab**2 + pz_cm**2)

    if np.any(p_cm <= 0.0):
        raise RuntimeError("Zero CM momentum was found for one or more alpha particles.")

    cos_theta_cm = (pz_cm / p_cm).reshape(-1)

    cos_theta_lab = np.clip(cos_theta_lab, -1.0, 1.0)
    cos_theta_cm = np.clip(cos_theta_cm, -1.0, 1.0)

    bins = np.linspace(-1.0, 1.0, 41)
    centers = 0.5 * (bins[:-1] + bins[1:])

    lab_hist, _ = np.histogram(cos_theta_lab, bins=bins, density=True)
    cm_hist, _ = np.histogram(cos_theta_cm, bins=bins, density=True)

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

    fig, ax = plt.subplots(figsize=(5, 4))

    ax.step(
        centers,
        lab_hist,
        where="mid",
        linewidth=2.2,
        color="#990000",
        label="LAB",
    )

    ax.step(
        centers,
        cm_hist,
        where="mid",
        linewidth=2.2,
        color="#003366",
        label="three-alpha c.m.",
    )

    ax.axhline(
        0.5,
        linestyle="--",
        linewidth=1.2,
        color="#777777",
    )

    ax.set_xlabel(r"$\cos\theta_{\alpha}$")
    ax.set_ylabel("normalised yield")
    ax.set_xlim(-1.0, 1.0)

    ymax = max(np.max(lab_hist), np.max(cm_hist))
    ymax = max(0.65, 1.08 * ymax)
    ax.set_ylim(0.0, ymax)

    ax.legend(frameon=False)

    fig.subplots_adjust(**FIXED_MARGINS)

    if save is not None:
        save = Path(save)
        save.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save, dpi=200)
        print(f"Saved: {save}")

    print(f"Selected events: {np.count_nonzero(valid)}")
    print(f"Filled entries per frame: {cos_theta_lab.size}")
    print(f"Mean CM beta: {np.mean(beta_cm):.6e}")

    plt.show()

    return fig


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot LAB and three-alpha CM angular distributions on the same axis "
            "for all three alpha particles in the 675-keV sequential channel."
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

    plot_675_alpha_angular(
        args.root_file,
        save=args.save,
    )


if __name__ == "__main__":
    main()
