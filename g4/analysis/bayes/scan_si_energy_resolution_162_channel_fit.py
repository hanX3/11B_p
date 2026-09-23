#!/usr/bin/env python3
"""Si-detector energy-resolution scan for the 162-keV Dalitz fit.

The script keeps the three pure benchmark templates ideal and unsmeared.
A Gaussian energy response is applied only to the three alpha energies in
the mixed closure sample:

    sigma_E = FWHM_E / 2.354820045.

The same standard-normal random numbers are reused for every FWHM point,
which suppresses artificial point-to-point fluctuations. The analysis uses
full laboratory angular coverage and allows fixed theta/phi granularities.

Output:
    one tightly stacked three-row fitted-fraction figure.

Dependencies:
    python3 -m pip install --user numpy matplotlib uproot
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

import uproot


# ---------------------------------------------------------------------------
# Physics and ROOT conventions
# ---------------------------------------------------------------------------

ALPHA_PDG = 1000020040
ALPHA_MASS_MEV = 3727.379378
SOURCE_H11B_REACTION_PRODUCT = 1

CHANNEL_IDS = np.array([1, 2, 3], dtype=int)
CHANNEL_KEYS = ("162seq_Be_gs", "162seq_Be2plus", "162direct")
CHANNEL_LABELS = (
    r"$\alpha_0$: $^{8}\mathrm{Be}(\mathrm{g.s.})$",
    r"$\alpha_1$: $^{8}\mathrm{Be}(2^{+})$",
    r"Direct $3\alpha$",
)

DEFAULT_COMPARISON_YLIMS = (
    (0.192, 0.201),
    (0.578, 0.603),
    (0.195, 0.230),
)

INPANEL_LABEL_POSITIONS = (
    (0.06, 0.12),  # alpha0: lower-left
    (0.06, 0.12),  # alpha1: lower-left
    (0.06, 0.86),  # direct: upper-left
)

SPHERE_BRANCHES = (
    "event_id",
    "pdg",
    "kinetic_energy_MeV",
    "px_MeV_c",
    "py_MeV_c",
    "pz_MeV_c",
    "h11b_particle_source",
    "generator_particle_index",
)

SCRIPT_VERSION = "2026-07-20-si-energy-resolution-v1"

SPHERE_POSITION_BRANCHES = (
    "x_mm",
    "y_mm",
    "z_mm",
)


# ---------------------------------------------------------------------------
# ROOT helpers
# ---------------------------------------------------------------------------

def object_names(root_file: Any) -> set[str]:
    return {str(name).split(";")[0] for name in root_file.keys()}


def branch_names(tree: Any) -> set[str]:
    return {str(name).split(";")[0] for name in tree.keys()}


def require_tree(root_file: Any, name: str) -> Any:
    if name not in object_names(root_file):
        raise RuntimeError(f"missing tree '{name}'")
    return root_file[name]


def require_branches(
    tree: Any,
    needed: tuple[str, ...],
    tree_name: str,
) -> None:
    available = branch_names(tree)
    missing = sorted(set(needed) - available)
    if missing:
        raise RuntimeError(
            f"tree '{tree_name}' is missing branches: "
            f"{', '.join(missing)}"
        )


def load_reaction_event_ids(root_path: Path) -> np.ndarray:
    with uproot.open(root_path) as root_file:
        tree = require_tree(root_file, "reaction")
        available = branch_names(tree)
        if "event" in available:
            branch = "event"
        elif "event_id" in available:
            branch = "event_id"
        else:
            raise RuntimeError(
                f"{root_path}: reaction tree has no event/event_id branch"
            )
        return np.asarray(
            tree[branch].array(library="np"),
            dtype=np.int64,
        )


def load_reaction_channels(root_path: Path) -> np.ndarray:
    with uproot.open(root_path) as root_file:
        tree = require_tree(root_file, "reaction")
        if "h11b_reaction_channel" not in branch_names(tree):
            raise RuntimeError(
                f"{root_path}: reaction tree is missing "
                "h11b_reaction_channel"
            )
        return np.asarray(
            tree["h11b_reaction_channel"].array(library="np"),
            dtype=int,
        )


def load_complete_three_alpha(
    root_path: Path,
    load_hit_position: bool = False,
) -> dict[str, np.ndarray]:
    """Load one ordered row triplet per complete H11B three-alpha event."""
    requested_branches = list(SPHERE_BRANCHES)
    if load_hit_position:
        requested_branches.extend(SPHERE_POSITION_BRANCHES)

    with uproot.open(root_path) as root_file:
        tree = require_tree(root_file, "virtual_sphere")
        require_branches(
            tree,
            tuple(requested_branches),
            "virtual_sphere",
        )
        data = tree.arrays(
            tuple(requested_branches),
            library="np",
        )

    selected = (
        (np.asarray(data["pdg"]) == ALPHA_PDG)
        & (
            np.asarray(data["h11b_particle_source"])
            == SOURCE_H11B_REACTION_PRODUCT
        )
        & (
            np.asarray(data["generator_particle_index"])
            >= 0
        )
        & (
            np.asarray(data["generator_particle_index"])
            <= 2
        )
    )
    if not np.any(selected):
        raise RuntimeError(
            f"{root_path}: no labeled H11B reaction-product alphas"
        )

    event = np.asarray(
        data["event_id"][selected],
        dtype=np.int64,
    )
    generator_index = np.asarray(
        data["generator_particle_index"][selected],
        dtype=np.int32,
    )
    kinetic = np.asarray(
        data["kinetic_energy_MeV"][selected],
        dtype=float,
    )
    momentum = np.column_stack(
        [
            np.asarray(
                data["px_MeV_c"][selected],
                dtype=float,
            ),
            np.asarray(
                data["py_MeV_c"][selected],
                dtype=float,
            ),
            np.asarray(
                data["pz_MeV_c"][selected],
                dtype=float,
            ),
        ]
    )

    hit_position = None
    if load_hit_position:
        hit_position = np.column_stack(
            [
                np.asarray(data["x_mm"][selected], dtype=float),
                np.asarray(data["y_mm"][selected], dtype=float),
                np.asarray(data["z_mm"][selected], dtype=float),
            ]
        )

    order = np.lexsort((generator_index, event))
    event = event[order]
    generator_index = generator_index[order]
    kinetic = kinetic[order]
    momentum = momentum[order]
    if hit_position is not None:
        hit_position = hit_position[order]

    unique_event, first, multiplicity = np.unique(
        event,
        return_index=True,
        return_counts=True,
    )
    starts = first[multiplicity == 3]
    candidate_events = unique_event[multiplicity == 3]
    if starts.size == 0:
        raise RuntimeError(
            f"{root_path}: no complete three-alpha events"
        )

    rows = starts[:, None] + np.arange(3, dtype=np.int64)[None, :]
    correct_indices = np.all(
        generator_index[rows]
        == np.array([0, 1, 2], dtype=generator_index.dtype),
        axis=1,
    )
    rows = rows[correct_indices]
    complete_event = candidate_events[correct_indices]
    if rows.size == 0:
        raise RuntimeError(
            f"{root_path}: no events with generator indices 0, 1 and 2"
        )

    result = {
        "event_id": np.asarray(complete_event, dtype=np.int64),
        "kinetic": np.asarray(kinetic[rows], dtype=float),
        "momentum": np.asarray(momentum[rows], dtype=float),
    }
    if hit_position is not None:
        result["hit_position"] = np.asarray(
            hit_position[rows],
            dtype=float,
        )
    return result


# ---------------------------------------------------------------------------
# Event quality and fixed event selection
# ---------------------------------------------------------------------------

def quality_mask(
    kinetic: np.ndarray,
    momentum: np.ndarray,
    mass_shell_tolerance: float,
) -> np.ndarray:
    finite = (
        np.all(np.isfinite(kinetic), axis=1)
        & np.all(np.isfinite(momentum), axis=(1, 2))
    )
    nonnegative = np.all(kinetic >= 0.0, axis=1)

    momentum2 = np.sum(momentum**2, axis=2)
    expected2 = kinetic * (
        kinetic + 2.0 * ALPHA_MASS_MEV
    )
    residual = np.abs(
        momentum2 - expected2
    ) / np.maximum(expected2, 1.0)
    mass_shell = (
        np.max(residual, axis=1)
        <= mass_shell_tolerance
    )
    direction_defined = np.all(
        np.sqrt(momentum2) > 0.0,
        axis=1,
    )
    return (
        finite
        & nonnegative
        & mass_shell
        & direction_defined
    )


def apply_mask(
    events: dict[str, np.ndarray],
    mask: np.ndarray,
) -> dict[str, np.ndarray]:
    result = {
        "event_id": events["event_id"][mask],
        "kinetic": events["kinetic"][mask],
        "momentum": events["momentum"][mask],
    }
    if "hit_position" in events:
        result["hit_position"] = events["hit_position"][mask]
    return result


def select_fixed_mixed_sample(
    events: dict[str, np.ndarray],
    reaction_event_ids: np.ndarray,
    n_fit: int,
    seed: int,
) -> dict[str, np.ndarray]:
    """Choose the same accepted closure events for every scan point."""
    accepted_ids = np.asarray(events["event_id"], dtype=np.int64)
    rng = np.random.default_rng(seed)
    shuffled_reaction_ids = rng.permutation(
        np.asarray(reaction_event_ids, dtype=np.int64)
    )
    accepted_in_order = np.isin(
        shuffled_reaction_ids,
        accepted_ids,
        assume_unique=True,
    )
    accepted_positions = np.flatnonzero(accepted_in_order)
    if accepted_positions.size < n_fit:
        raise RuntimeError(
            f"Requested N_fit={n_fit:,}, but only "
            f"{accepted_positions.size:,} accepted complete events exist"
        )

    generated_prefix = shuffled_reaction_ids[
        : accepted_positions[n_fit - 1] + 1
    ]
    selected_mask = np.isin(
        accepted_ids,
        generated_prefix,
        assume_unique=True,
    )
    selected_events = apply_mask(events, selected_mask)
    if selected_events["event_id"].size != n_fit:
        raise RuntimeError(
            "Internal fixed-sample selection mismatch: "
            f"expected {n_fit}, found "
            f"{selected_events['event_id'].size}"
        )

    print(
        f"Fixed mixed sample: N_generated={generated_prefix.size:,}, "
        f"N_fit={n_fit:,}"
    )
    return selected_events


# ---------------------------------------------------------------------------
# Angular quantization and momentum reconstruction
# ---------------------------------------------------------------------------

def unit_directions(momentum: np.ndarray) -> np.ndarray:
    magnitude = np.linalg.norm(momentum, axis=2)
    direction = np.full_like(momentum, np.nan, dtype=float)
    np.divide(
        momentum,
        magnitude[:, :, None],
        out=direction,
        where=magnitude[:, :, None] > 0.0,
    )
    return direction


def quantize_theta(
    theta_deg: np.ndarray,
    width_deg: float,
    offset_deg: float,
) -> np.ndarray:
    if width_deg <= 0.0:
        return np.asarray(theta_deg, dtype=float).copy()

    shifted = np.asarray(theta_deg, dtype=float) - offset_deg
    cell = np.floor(shifted / width_deg)
    centre = offset_deg + (cell + 0.5) * width_deg
    return np.clip(centre, 0.0, 180.0)


def quantize_phi(
    phi_deg: np.ndarray,
    width_deg: float,
    offset_deg: float,
) -> np.ndarray:
    phi = np.mod(
        np.asarray(phi_deg, dtype=float),
        360.0,
    )
    if width_deg <= 0.0:
        return phi

    shifted = np.mod(phi - offset_deg, 360.0)
    cell = np.floor(shifted / width_deg)
    centre = offset_deg + (cell + 0.5) * width_deg
    return np.mod(centre, 360.0)


def reconstruct_with_granularity(
    events: dict[str, np.ndarray],
    theta_width_deg: float,
    phi_width_deg: float,
    theta_offset_deg: float,
    phi_offset_deg: float,
    angle_source: str,
    sphere_center_mm: np.ndarray,
) -> tuple[dict[str, np.ndarray], dict[str, float]]:
    """Rebuild on-shell momenta using quantized reconstructed directions."""
    kinetic = np.asarray(events["kinetic"], dtype=float)

    if angle_source == "momentum":
        direction = unit_directions(
            np.asarray(events["momentum"], dtype=float)
        )
    elif angle_source == "sphere-center":
        if "hit_position" not in events:
            raise RuntimeError(
                "sphere-center reconstruction requires x_mm, y_mm and z_mm"
            )
        displacement = (
            np.asarray(events["hit_position"], dtype=float)
            - sphere_center_mm[None, None, :]
        )
        direction = unit_directions(displacement)
    else:
        raise RuntimeError(
            f"unsupported angle source: {angle_source}"
        )

    theta = np.rad2deg(
        np.arccos(
            np.clip(direction[:, :, 2], -1.0, 1.0)
        )
    )
    phi = np.mod(
        np.rad2deg(
            np.arctan2(
                direction[:, :, 1],
                direction[:, :, 0],
            )
        ),
        360.0,
    )

    theta_rec = quantize_theta(
        theta,
        theta_width_deg,
        theta_offset_deg,
    )
    phi_rec = quantize_phi(
        phi,
        phi_width_deg,
        phi_offset_deg,
    )

    theta_rad = np.deg2rad(theta_rec)
    phi_rad = np.deg2rad(phi_rec)
    reconstructed_direction = np.empty_like(direction)
    reconstructed_direction[:, :, 0] = (
        np.sin(theta_rad) * np.cos(phi_rad)
    )
    reconstructed_direction[:, :, 1] = (
        np.sin(theta_rad) * np.sin(phi_rad)
    )
    reconstructed_direction[:, :, 2] = np.cos(theta_rad)

    momentum_magnitude = np.sqrt(
        kinetic * (
            kinetic + 2.0 * ALPHA_MASS_MEV
        )
    )
    reconstructed_momentum = (
        reconstructed_direction
        * momentum_magnitude[:, :, None]
    )

    delta_theta = theta_rec - theta
    delta_phi = (
        phi_rec - phi + 180.0
    ) % 360.0 - 180.0

    diagnostics = {
        "theta_rms_deg": float(
            np.sqrt(np.mean(delta_theta**2))
        ),
        "phi_rms_deg": float(
            np.sqrt(np.mean(delta_phi**2))
        ),
        "theta_max_abs_deg": float(
            np.max(np.abs(delta_theta))
        ),
        "phi_max_abs_deg": float(
            np.max(np.abs(delta_phi))
        ),
    }

    return (
        {
            "event_id": events["event_id"],
            "kinetic": kinetic,
            "momentum": reconstructed_momentum,
        },
        diagnostics,
    )


# ---------------------------------------------------------------------------
# Four-momentum CM boost and ordered Dalitz observable
# ---------------------------------------------------------------------------

def three_alpha_cm_kinetic(
    kinetic: np.ndarray,
    momentum: np.ndarray,
) -> np.ndarray:
    energy = kinetic + ALPHA_MASS_MEV
    total_energy = np.sum(energy, axis=1)
    total_momentum = np.sum(momentum, axis=1)

    beta = total_momentum / total_energy[:, None]
    beta2 = np.sum(beta**2, axis=1)
    valid = (
        np.all(np.isfinite(energy), axis=1)
        & np.all(np.isfinite(momentum), axis=(1, 2))
        & np.isfinite(total_energy)
        & (total_energy > 0.0)
        & (beta2 < 1.0 - 1.0e-14)
    )
    if not np.all(valid):
        raise RuntimeError(
            f"{np.count_nonzero(~valid)} reconstructed events have "
            "invalid total four-momentum"
        )

    gamma = 1.0 / np.sqrt(1.0 - beta2)
    beta_dot_p = np.einsum(
        "ni,nji->nj",
        beta,
        momentum,
    )
    energy_cm = gamma[:, None] * (
        energy - beta_dot_p
    )
    kinetic_cm = energy_cm - ALPHA_MASS_MEV

    if np.any(kinetic_cm < -1.0e-7):
        minimum = float(np.min(kinetic_cm))
        raise RuntimeError(
            "Negative reconstructed CM kinetic energy beyond numerical "
            f"roundoff: minimum={minimum:.6e} MeV"
        )
    return np.maximum(kinetic_cm, 0.0)


def ordered_dalitz_xy(
    events: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    kinetic_cm = three_alpha_cm_kinetic(
        events["kinetic"],
        events["momentum"],
    )
    ordered = np.sort(kinetic_cm, axis=1)[:, ::-1]
    high, middle, low = ordered.T
    energy_sum = high + middle + low

    x = (
        np.sqrt(3.0)
        * (middle - low)
        / energy_sum
    )
    y = (
        2.0 * high - middle - low
    ) / energy_sum

    finite = np.isfinite(x) & np.isfinite(y)
    if not np.all(finite):
        raise RuntimeError(
            f"{np.count_nonzero(~finite)} non-finite Dalitz points"
        )
    return x, y


def dalitz_histogram(
    events: dict[str, np.ndarray],
    bins: int,
    limits: tuple[float, float],
) -> np.ndarray:
    x, y = ordered_dalitz_xy(events)
    histogram, _, _ = np.histogram2d(
        x,
        y,
        bins=bins,
        range=[
            [limits[0], limits[1]],
            [limits[0], limits[1]],
        ],
    )
    if int(np.sum(histogram)) != events["event_id"].size:
        raise RuntimeError(
            "Some Dalitz points fell outside the histogram range. "
            "Increase --dalitz-limit."
        )
    return histogram.astype(float)


# ---------------------------------------------------------------------------
# Bayesian template fit with adaptive local posterior grid
# ---------------------------------------------------------------------------

def simplex_grid(grid_size: int) -> np.ndarray:
    values = np.linspace(0.0, 1.0, grid_size)
    f1, f2 = np.meshgrid(
        values,
        values,
        indexing="ij",
    )
    valid = f1 + f2 <= 1.0 + 1.0e-12
    return np.column_stack(
        [
            f1[valid],
            f2[valid],
            np.maximum(
                1.0 - f1[valid] - f2[valid],
                0.0,
            ),
        ]
    )


def evaluate_posterior(
    fractions: np.ndarray,
    template_probability: np.ndarray,
    efficiencies: np.ndarray,
    data_counts: np.ndarray,
    batch_size: int = 512,
) -> tuple[np.ndarray, np.ndarray]:
    log_likelihood = np.empty(
        fractions.shape[0],
        dtype=float,
    )

    for start in range(
        0,
        fractions.shape[0],
        batch_size,
    ):
        stop = min(
            start + batch_size,
            fractions.shape[0],
        )
        generator_fraction = fractions[start:stop]
        accepted_weight = (
            generator_fraction
            * efficiencies[None, :]
        )
        accepted_weight /= np.sum(
            accepted_weight,
            axis=1,
            keepdims=True,
        )
        bin_probability = (
            accepted_weight @ template_probability
        )
        bin_probability = np.clip(
            bin_probability,
            1.0e-300,
            None,
        )
        log_likelihood[start:stop] = (
            np.log(bin_probability) @ data_counts
        )

    shifted = log_likelihood - np.max(log_likelihood)
    posterior = np.exp(shifted)
    posterior /= np.sum(posterior)
    return posterior, log_likelihood


def refine_map(
    initial: np.ndarray,
    global_step: float,
    template_probability: np.ndarray,
    efficiencies: np.ndarray,
    data_counts: np.ndarray,
    refinements: int = 3,
    grid_size: int = 81,
) -> np.ndarray:
    centre = np.asarray(initial, dtype=float).copy()
    half_width = 2.0 * global_step

    for _ in range(refinements):
        f1_axis = np.linspace(
            max(0.0, centre[0] - half_width),
            min(1.0, centre[0] + half_width),
            grid_size,
        )
        f2_axis = np.linspace(
            max(0.0, centre[1] - half_width),
            min(1.0, centre[1] + half_width),
            grid_size,
        )
        f1, f2 = np.meshgrid(
            f1_axis,
            f2_axis,
            indexing="ij",
        )
        valid = f1 + f2 <= 1.0 + 1.0e-14
        trial = np.column_stack(
            [
                f1[valid],
                f2[valid],
                1.0 - f1[valid] - f2[valid],
            ]
        )
        _, log_likelihood = evaluate_posterior(
            trial,
            template_probability,
            efficiencies,
            data_counts,
        )
        centre = trial[int(np.argmax(log_likelihood))]
        half_width /= 8.0

    return centre


def local_covariance(
    mode: np.ndarray,
    global_step: float,
    template_probability: np.ndarray,
    efficiencies: np.ndarray,
    data_counts: np.ndarray,
) -> np.ndarray | None:
    h = max(global_step / 100.0, 1.0e-6)
    f1, f2 = float(mode[0]), float(mode[1])
    if min(f1, f2, 1.0 - f1 - f2) <= 2.5 * h:
        return None

    points2 = np.array(
        [
            [f1, f2],
            [f1 + h, f2],
            [f1 - h, f2],
            [f1, f2 + h],
            [f1, f2 - h],
            [f1 + h, f2 + h],
            [f1 + h, f2 - h],
            [f1 - h, f2 + h],
            [f1 - h, f2 - h],
        ],
        dtype=float,
    )
    points = np.column_stack(
        [
            points2,
            1.0 - points2[:, 0] - points2[:, 1],
        ]
    )
    _, log_likelihood = evaluate_posterior(
        points,
        template_probability,
        efficiencies,
        data_counts,
    )
    nll = -log_likelihood

    h11 = (
        nll[1] - 2.0 * nll[0] + nll[2]
    ) / h**2
    h22 = (
        nll[3] - 2.0 * nll[0] + nll[4]
    ) / h**2
    h12 = (
        nll[5] - nll[6] - nll[7] + nll[8]
    ) / (4.0 * h**2)

    hessian = np.array(
        [[h11, h12], [h12, h22]],
        dtype=float,
    )
    if (
        not np.all(np.isfinite(hessian))
        or np.any(np.linalg.eigvalsh(hessian) <= 0.0)
    ):
        return None
    return np.linalg.inv(hessian)


def adaptive_local_grid(
    mode: np.ndarray,
    covariance: np.ndarray | None,
    global_step: float,
    sigma_span: float,
    grid_size: int,
) -> np.ndarray:
    if covariance is None:
        half1 = half2 = max(
            10.0 * global_step,
            0.15,
        )
    else:
        sigma = np.sqrt(
            np.maximum(
                np.diag(covariance),
                0.0,
            )
        )
        half1 = max(
            sigma_span * sigma[0],
            global_step / 5.0,
        )
        half2 = max(
            sigma_span * sigma[1],
            global_step / 5.0,
        )

    f1_axis = np.linspace(
        max(0.0, mode[0] - half1),
        min(1.0, mode[0] + half1),
        grid_size,
    )
    f2_axis = np.linspace(
        max(0.0, mode[1] - half2),
        min(1.0, mode[1] + half2),
        grid_size,
    )
    f1, f2 = np.meshgrid(
        f1_axis,
        f2_axis,
        indexing="ij",
    )
    valid = f1 + f2 <= 1.0 + 1.0e-14
    return np.column_stack(
        [
            f1[valid],
            f2[valid],
            1.0 - f1[valid] - f2[valid],
        ]
    )


def weighted_quantiles(
    values: np.ndarray,
    weights: np.ndarray,
    probabilities: tuple[float, ...],
) -> np.ndarray:
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cumulative = np.cumsum(weights)
    cumulative /= cumulative[-1]
    return np.interp(
        np.asarray(probabilities),
        cumulative,
        values,
    )


def fit_histogram(
    data_histogram: np.ndarray,
    template_histograms: np.ndarray,
    efficiencies: np.ndarray,
    global_fractions: np.ndarray,
    grid_size: int,
    local_grid_size: int,
    local_sigma_span: float,
    template_pseudocount: float,
) -> dict[str, np.ndarray]:
    active = (
        np.any(template_histograms > 0.0, axis=0)
        | (data_histogram > 0.0)
    )
    data_counts = data_histogram[active].astype(np.int64)
    template_counts = template_histograms[:, active]

    template_probability = (
        template_counts + template_pseudocount
    )
    template_probability /= np.sum(
        template_probability,
        axis=1,
        keepdims=True,
    )

    _, global_log_likelihood = evaluate_posterior(
        global_fractions,
        template_probability,
        efficiencies,
        data_counts,
    )
    global_map = global_fractions[
        int(np.argmax(global_log_likelihood))
    ]
    global_step = 1.0 / (grid_size - 1)

    refined_map = refine_map(
        global_map,
        global_step,
        template_probability,
        efficiencies,
        data_counts,
    )
    covariance = local_covariance(
        refined_map,
        global_step,
        template_probability,
        efficiencies,
        data_counts,
    )
    local_fractions = adaptive_local_grid(
        refined_map,
        covariance,
        global_step,
        local_sigma_span,
        local_grid_size,
    )
    posterior, log_likelihood = evaluate_posterior(
        local_fractions,
        template_probability,
        efficiencies,
        data_counts,
    )

    median = np.empty(3)
    lower95 = np.empty(3)
    upper95 = np.empty(3)
    for channel in range(3):
        (
            lower95[channel],
            median[channel],
            upper95[channel],
        ) = weighted_quantiles(
            local_fractions[:, channel],
            posterior,
            (0.025, 0.50, 0.975),
        )

    return {
        "median": median,
        "lower95": lower95,
        "upper95": upper95,
        "map": local_fractions[
            int(np.argmax(log_likelihood))
        ],
    }




# ---------------------------------------------------------------------------
# Generated-prefix truth
# ---------------------------------------------------------------------------

def prepare_reaction_channel_lookup(
    reaction_event_ids: np.ndarray,
    reaction_channels: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Prepare an exact event-id to reaction-channel lookup."""
    event_ids = np.asarray(
        reaction_event_ids,
        dtype=np.int64,
    )
    channels = np.asarray(
        reaction_channels,
        dtype=int,
    )
    if event_ids.shape != channels.shape:
        raise RuntimeError(
            "reaction event-id and channel arrays have different shapes"
        )
    if event_ids.size == 0:
        raise RuntimeError(
            "reaction tree contains no events"
        )

    order = np.argsort(event_ids)
    sorted_ids = event_ids[order]
    sorted_channels = channels[order]
    if np.any(np.diff(sorted_ids) == 0):
        raise RuntimeError(
            "reaction event IDs are not unique"
        )
    return sorted_ids, sorted_channels


def generated_prefix_truth(
    shuffled_reaction_ids: np.ndarray,
    generated_count: int,
    sorted_reaction_ids: np.ndarray,
    sorted_reaction_channels: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return channel counts and fractions for the actual generated prefix."""
    if generated_count <= 0:
        raise RuntimeError(
            "generated_count must be positive"
        )

    prefix_ids = np.asarray(
        shuffled_reaction_ids[:generated_count],
        dtype=np.int64,
    )
    positions = np.searchsorted(
        sorted_reaction_ids,
        prefix_ids,
    )
    valid = (
        (positions >= 0)
        & (positions < sorted_reaction_ids.size)
    )
    if not np.all(valid):
        raise RuntimeError(
            "generated prefix contains event IDs absent from reaction lookup"
        )
    if not np.array_equal(
        sorted_reaction_ids[positions],
        prefix_ids,
    ):
        raise RuntimeError(
            "reaction-channel lookup failed for generated prefix"
        )

    prefix_channels = sorted_reaction_channels[positions]
    counts = np.asarray(
        [
            np.count_nonzero(
                prefix_channels == channel_id
            )
            for channel_id in CHANNEL_IDS
        ],
        dtype=np.int64,
    )
    recognized = int(np.sum(counts))
    if recognized != generated_count:
        unknown = generated_count - recognized
        raise RuntimeError(
            f"generated prefix contains {unknown} events outside "
            "the three fitted 162-keV channels"
        )
    fractions = counts.astype(float) / float(recognized)
    return counts, fractions


# ---------------------------------------------------------------------------
# One-dimensional response handling
# ---------------------------------------------------------------------------

def apply_closure_energy_response(
    reconstructed: dict[str, np.ndarray],
    energy_fwhm_kev: float,
    threshold_kev: float,
    standard_normal: np.ndarray,
) -> tuple[dict[str, np.ndarray], np.ndarray, dict[str, float]]:
    """Smear closure alpha energies, apply threshold, and rebuild momenta."""
    true_kinetic = np.asarray(
        reconstructed["kinetic"],
        dtype=float,
    )
    direction = unit_directions(
        np.asarray(reconstructed["momentum"], dtype=float)
    )

    sigma_mev = (
        float(energy_fwhm_kev)
        / 1000.0
        / 2.3548200450309493
    )
    measured_kinetic = (
        true_kinetic
        + sigma_mev * np.asarray(standard_normal, dtype=float)
    )

    threshold_mev = float(threshold_kev) / 1000.0
    finite = np.all(
        np.isfinite(measured_kinetic),
        axis=1,
    )
    accepted = (
        finite
        & np.all(
            measured_kinetic >= threshold_mev,
            axis=1,
        )
    )

    nonnegative_kinetic = np.maximum(
        measured_kinetic,
        0.0,
    )
    magnitude = np.sqrt(
        nonnegative_kinetic
        * (
            nonnegative_kinetic
            + 2.0 * ALPHA_MASS_MEV
        )
    )
    measured_momentum = (
        direction * magnitude[:, :, None]
    )

    result = {
        "event_id": np.asarray(
            reconstructed["event_id"],
            dtype=np.int64,
        )[accepted],
        "kinetic": measured_kinetic[accepted],
        "momentum": measured_momentum[accepted],
    }
    diagnostics = {
        "energy_sigma_kev": float(
            sigma_mev * 1000.0
        ),
        "threshold_survival": float(
            np.count_nonzero(accepted)
            / max(accepted.size, 1)
        ),
        "available_after_threshold": float(
            np.count_nonzero(accepted)
        ),
    }
    return result, accepted, diagnostics


def select_exact_nfit(
    accepted_events: dict[str, np.ndarray],
    shuffled_reaction_ids: np.ndarray,
    n_fit: int,
) -> tuple[dict[str, np.ndarray], int]:
    """Select the generated-event prefix yielding exactly n_fit events."""
    accepted_ids = np.asarray(
        accepted_events["event_id"],
        dtype=np.int64,
    )
    accepted_in_order = np.isin(
        shuffled_reaction_ids,
        accepted_ids,
        assume_unique=True,
    )
    positions = np.flatnonzero(accepted_in_order)
    if positions.size < n_fit:
        raise RuntimeError(
            f"Requested N_fit={n_fit:,}, but only "
            f"{positions.size:,} events survive this response point"
        )

    generated_count = int(positions[n_fit - 1] + 1)
    generated_prefix = shuffled_reaction_ids[
        :generated_count
    ]
    selected_mask = np.isin(
        accepted_ids,
        generated_prefix,
        assume_unique=True,
    )
    selected = apply_mask(
        accepted_events,
        selected_mask,
    )
    actual = int(selected["event_id"].size)
    if actual != n_fit:
        raise RuntimeError(
            "Internal exact-statistics mismatch: "
            f"expected {n_fit}, found {actual}"
        )
    return selected, generated_count


def response_parameters(
    scan_name: str,
    value: float,
    fixed_theta_deg: float,
    fixed_phi_deg: float,
    fixed_energy_fwhm_kev: float,
    fixed_threshold_kev: float,
) -> tuple[float, float, float, float]:
    theta = float(fixed_theta_deg)
    phi = float(fixed_phi_deg)
    fwhm = float(fixed_energy_fwhm_kev)
    threshold = float(fixed_threshold_kev)

    if scan_name == "theta":
        theta = float(value)
    elif scan_name == "phi":
        phi = float(value)
    elif scan_name == "energy-fwhm":
        fwhm = float(value)
    elif scan_name == "threshold":
        threshold = float(value)
    else:
        raise RuntimeError(
            f"Unsupported scan parameter: {scan_name}"
        )
    return theta, phi, fwhm, threshold


def default_scan_values(scan_name: str) -> np.ndarray:
    if scan_name == "theta":
        return np.asarray(
            [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0],
            dtype=float,
        )
    if scan_name == "phi":
        return np.asarray(
            [0.0, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0],
            dtype=float,
        )
    if scan_name == "energy-fwhm":
        return np.asarray(
            [0.0, 20.0, 50.0, 100.0, 150.0, 200.0],
            dtype=float,
        )
    if scan_name == "threshold":
        return np.asarray(
            [0.0, 50.0, 100.0, 150.0, 200.0, 300.0, 500.0],
            dtype=float,
        )
    raise RuntimeError(
        f"Unsupported scan parameter: {scan_name}"
    )


def scan_axis_metadata(
    scan_name: str,
) -> tuple[str, str, str]:
    if scan_name == "theta":
        return (
            r"$\Delta\theta$ granularity (deg)",
            "dtheta",
            "Angular granularity",
        )
    if scan_name == "phi":
        return (
            r"$\Delta\phi$ granularity (deg)",
            "dphi",
            "Angular granularity",
        )
    if scan_name == "energy-fwhm":
        return (
            r"Alpha-energy FWHM (keV)",
            "energy_fwhm",
            "Closure energy resolution",
        )
    if scan_name == "threshold":
        return (
            r"Per-alpha energy threshold (keV)",
            "threshold",
            "Closure energy threshold",
        )
    raise RuntimeError(
        f"Unsupported scan parameter: {scan_name}"
    )


# ---------------------------------------------------------------------------
# One-dimensional plotting and CSV output
# ---------------------------------------------------------------------------

def configure_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.labelsize": 14,
            "axes.titlesize": 14,
            "legend.fontsize": 10,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "axes.linewidth": 1.1,
            "lines.linewidth": 1.7,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )


def plot_fraction_scan(
    values: np.ndarray,
    median: np.ndarray,
    lower95: np.ndarray,
    upper95: np.ndarray,
    truth: np.ndarray,
    xlabel: str,
    fixed_ylims: tuple[tuple[float, float], ...] | None = None,
) -> None:
    colors = ("#0072B2", "#D55E00", "#009E73")

    truth_array = np.asarray(truth, dtype=float)
    if truth_array.ndim == 1:
        if truth_array.shape != (3,):
            raise RuntimeError(
                "one-dimensional truth must contain exactly three fractions"
            )
        truth_array = np.repeat(
            truth_array[None, :],
            values.size,
            axis=0,
        )
    if truth_array.shape != median.shape:
        raise RuntimeError(
            "truth array must have shape (number of scan points, 3)"
        )

    figure, axes = plt.subplots(
        3,
        1,
        figsize=(7.2, 8.6),
        sharex=True,
    )

    for channel, axis in enumerate(axes):
        y = median[:, channel]
        lower = lower95[:, channel]
        upper = upper95[:, channel]
        truth_channel = truth_array[:, channel]

        axis.errorbar(
            values,
            y,
            yerr=np.vstack(
                [
                    y - lower,
                    upper - y,
                ]
            ),
            marker="o",
            linestyle="-",
            color=colors[channel],
            capsize=3,
            markersize=5,
            label="Posterior median and 95% CI",
        )
        axis.plot(
            values,
            truth_channel,
            color="black",
            linestyle=":",
            linewidth=1.4,
            label="Generated-prefix truth",
        )

        label_x, label_y = INPANEL_LABEL_POSITIONS[channel]
        label_va = "top" if channel == 2 else "bottom"
        axis.text(
            label_x,
            label_y,
            CHANNEL_LABELS[channel],
            transform=axis.transAxes,
            ha="left",
            va=label_va,
            fontsize=16,
            bbox=dict(
                facecolor="white",
                edgecolor="none",
                alpha=0.75,
                pad=2.0,
            ),
        )
        axis.grid(alpha=0.25)
        if values.size >= 2:
            min_step = float(np.min(np.diff(np.sort(values))))
            xpad = max(
                0.35 * min_step,
                0.02 * (float(np.max(values)) - float(np.min(values))),
            )
        else:
            xpad = 0.5
        axis.set_xlim(
            float(np.min(values)) - xpad,
            float(np.max(values)) + xpad,
        )

        if fixed_ylims is not None:
            ymin, ymax = fixed_ylims[channel]
            axis.set_ylim(
                float(ymin),
                float(ymax),
            )
        else:
            low_limit = min(
                float(np.min(lower)),
                float(np.min(truth_channel)),
            )
            high_limit = max(
                float(np.max(upper)),
                float(np.max(truth_channel)),
            )
            span = max(
                high_limit - low_limit,
                1.0e-3,
            )
            axis.set_ylim(
                max(0.0, low_limit - 0.15 * span),
                min(1.0, high_limit + 0.15 * span),
            )

    figure.supxlabel(
        xlabel,
        y=0.035,
        fontsize=15,
    )
    figure.supylabel(
        "Fitted generator fraction",
        x=0.02,
        fontsize=15,
    )

    figure.subplots_adjust(
        left=0.14,
        right=0.97,
        bottom=0.08,
        top=0.985,
        hspace=0.12,
    )
    return figure




# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Scan the effect of Si-detector alpha-energy resolution on the "
            "extracted 162-keV channel fractions. Pure benchmark templates "
            "remain ideal and unsmeared; Gaussian energy broadening is applied "
            "only to the mixed closure sample."
        )
    )
    parser.add_argument(
        "seq_be_gs",
        type=Path,
        help="pure 162seq_Be_gs ROOT file",
    )
    parser.add_argument(
        "seq_be2plus",
        type=Path,
        help="pure 162seq_Be2plus ROOT file",
    )
    parser.add_argument(
        "direct",
        type=Path,
        help="pure 162direct ROOT file",
    )
    parser.add_argument(
        "mixed",
        type=Path,
        help="mixed closure ROOT file",
    )

    parser.add_argument(
        "--fwhm-values-keV",
        type=float,
        nargs="+",
        default=[
            0.0,
            20.0,
            50.0,
            100.0,
            150.0,
            200.0,
            300.0,
        ],
        metavar="FWHM",
        help=(
            "Si alpha-energy FWHM values in keV "
            "(default: 0 20 50 100 150 200 300)"
        ),
    )
    parser.add_argument(
        "--theta-granularity-deg",
        type=float,
        default=0.0,
        help="fixed closure Delta-theta granularity (default: 0 deg)",
    )
    parser.add_argument(
        "--phi-granularity-deg",
        type=float,
        default=0.0,
        help="fixed closure Delta-phi granularity (default: 0 deg)",
    )
    parser.add_argument(
        "--energy-threshold-keV",
        type=float,
        default=0.0,
        help="fixed per-alpha measured-energy threshold (default: 0 keV)",
    )

    parser.add_argument(
        "--closure-angle-source",
        choices=(
            "auto",
            "momentum",
            "sphere-center",
        ),
        default="auto",
        help=(
            "closure angle source. 'auto' uses momentum only for a point-beam "
            "file with Delta-theta=Delta-phi=0; otherwise sphere-center "
            "(default: auto)"
        ),
    )
    parser.add_argument(
        "--theta-offset-deg",
        type=float,
        default=0.0,
        help="theta segmentation origin (default: 0 deg)",
    )
    parser.add_argument(
        "--phi-offset-deg",
        type=float,
        default=0.0,
        help="phi segmentation origin (default: 0 deg)",
    )
    parser.add_argument(
        "--sphere-center-x-mm",
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "--sphere-center-y-mm",
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "--sphere-center-z-mm",
        type=float,
        default=170.0,
    )

    parser.add_argument(
        "--n-fit",
        type=int,
        default=100000,
        help=(
            "exact accepted complete three-alpha events at every FWHM point "
            "(default: 100000)"
        ),
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=162,
        help="generated-event ordering seed (default: 162)",
    )
    parser.add_argument(
        "--energy-smear-seed",
        type=int,
        default=1620,
        help=(
            "common standard-normal energy-smearing seed "
            "(default: 1620)"
        ),
    )

    parser.add_argument(
        "--mass-shell-tolerance",
        type=float,
        default=1.0e-4,
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=60,
    )
    parser.add_argument(
        "--dalitz-limit",
        type=float,
        default=1.05,
    )
    parser.add_argument(
        "--grid-size",
        type=int,
        default=101,
    )
    parser.add_argument(
        "--local-grid-size",
        type=int,
        default=201,
    )
    parser.add_argument(
        "--local-sigma-span",
        type=float,
        default=6.0,
    )
    parser.add_argument(
        "--template-pseudocount",
        type=float,
        default=0.5,
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "scan_162_si_energy_resolution_fractions.png"
        ),
        help="output PNG",
    )
    parser.add_argument(
        "--auto-y-range",
        action="store_true",
        help="use independent automatic y-axis limits",
    )
    parser.add_argument(
        "--alpha0-ylim",
        type=float,
        nargs=2,
        metavar=("YMIN", "YMAX"),
        default=(0.185, 0.205),
        help="fixed alpha0 y-axis range (default: 0.185 0.205)",
    )
    parser.add_argument(
        "--alpha1-ylim",
        type=float,
        nargs=2,
        metavar=("YMIN", "YMAX"),
        default=(0.540, 0.605),
        help="fixed alpha1 y-axis range (default: 0.540 0.605)",
    )
    parser.add_argument(
        "--direct-ylim",
        type=float,
        nargs=2,
        metavar=("YMIN", "YMAX"),
        default=(0.195, 0.270),
        help="fixed direct-channel y-axis range (default: 0.195 0.270)",
    )
    return parser.parse_args(argv)


def infer_beam_fwhm_mm(
    path: Path,
) -> float:
    import re

    name = path.stem.lower()
    match = re.search(
        r"beam[_-]?fwhm([0-9]+(?:\.[0-9]+)?)mm",
        name,
    )
    if match:
        return float(match.group(1))
    return 0.0


def validate_args(
    args: argparse.Namespace,
) -> np.ndarray:
    for path in (
        args.seq_be_gs,
        args.seq_be2plus,
        args.direct,
        args.mixed,
    ):
        if not path.is_file():
            raise FileNotFoundError(f"Input file not found: {path}")

    values = np.unique(
        np.asarray(
            args.fwhm_values_keV,
            dtype=float,
        )
    )
    if (
        values.size == 0
        or np.any(~np.isfinite(values))
        or np.any(values < 0.0)
    ):
        raise ValueError(
            "--fwhm-values-keV must contain finite non-negative values"
        )

    fixed_values = (
        args.theta_granularity_deg,
        args.phi_granularity_deg,
        args.energy_threshold_keV,
    )
    if any(
        (not np.isfinite(value)) or value < 0.0
        for value in fixed_values
    ):
        raise ValueError(
            "fixed response parameters must be finite and non-negative"
        )
    if args.theta_granularity_deg > 180.0:
        raise ValueError(
            "--theta-granularity-deg must not exceed 180"
        )
    if args.phi_granularity_deg > 360.0:
        raise ValueError(
            "--phi-granularity-deg must not exceed 360"
        )
    if args.n_fit <= 0:
        raise ValueError("--n-fit must be positive")
    if args.bins < 10:
        raise ValueError("--bins must be at least 10")
    if args.grid_size < 21:
        raise ValueError("--grid-size must be at least 21")
    if args.local_grid_size < 51:
        raise ValueError("--local-grid-size must be at least 51")
    if args.template_pseudocount <= 0.0:
        raise ValueError(
            "--template-pseudocount must be positive"
        )

    for name, limits in (
        ("--alpha0-ylim", args.alpha0_ylim),
        ("--alpha1-ylim", args.alpha1_ylim),
        ("--direct-ylim", args.direct_ylim),
    ):
        ymin, ymax = map(float, limits)
        if (
            not np.isfinite(ymin)
            or not np.isfinite(ymax)
            or ymax <= ymin
        ):
            raise ValueError(
                f"{name} requires finite YMIN YMAX with YMAX > YMIN"
            )

    return values


def resolve_closure_angle_source(
    args: argparse.Namespace,
) -> str:
    if args.closure_angle_source != "auto":
        return args.closure_angle_source

    beam_fwhm = infer_beam_fwhm_mm(
        args.mixed
    )
    if (
        abs(beam_fwhm) < 1.0e-12
        and abs(args.theta_granularity_deg) < 1.0e-12
        and abs(args.phi_granularity_deg) < 1.0e-12
    ):
        return "momentum"
    return "sphere-center"


def _run(args, save=False, outdir="."):
    print(f"Script version: {SCRIPT_VERSION}")
    fwhm_values = validate_args(args)
    configure_plot_style()

    sphere_center_mm = np.asarray(
        [
            args.sphere_center_x_mm,
            args.sphere_center_y_mm,
            args.sphere_center_z_mm,
        ],
        dtype=float,
    )
    limits = (
        -float(args.dalitz_limit),
        float(args.dalitz_limit),
    )

    pure_paths = (
        args.seq_be_gs,
        args.seq_be2plus,
        args.direct,
    )

    print("Loading pure benchmark samples...")
    pure_events: list[dict[str, np.ndarray]] = []
    efficiencies = np.empty(3, dtype=float)
    for channel, path in enumerate(pure_paths):
        raw = load_complete_three_alpha(
            path,
            load_hit_position=False,
        )
        mask = quality_mask(
            raw["kinetic"],
            raw["momentum"],
            args.mass_shell_tolerance,
        )
        accepted = apply_mask(
            raw,
            mask,
        )
        n_reaction = load_reaction_event_ids(
            path
        ).size
        efficiencies[channel] = (
            accepted["event_id"].size
            / float(n_reaction)
        )
        pure_events.append(accepted)
        print(
            f"  {CHANNEL_KEYS[channel]:<22} "
            f"accepted={accepted['event_id'].size:,}, "
            f"efficiency={efficiencies[channel]:.6%}"
        )

    print("Loading mixed closure sample...")
    mixed_raw = load_complete_three_alpha(
        args.mixed,
        load_hit_position=True,
    )
    mixed_mask = quality_mask(
        mixed_raw["kinetic"],
        mixed_raw["momentum"],
        args.mass_shell_tolerance,
    )
    mixed_events = apply_mask(
        mixed_raw,
        mixed_mask,
    )

    reaction_event_ids = load_reaction_event_ids(
        args.mixed
    )
    reaction_channels = load_reaction_channels(
        args.mixed
    )
    (
        sorted_reaction_ids,
        sorted_reaction_channels,
    ) = prepare_reaction_channel_lookup(
        reaction_event_ids,
        reaction_channels,
    )

    rng_order = np.random.default_rng(
        args.seed
    )
    shuffled_reaction_ids = rng_order.permutation(
        reaction_event_ids
    )

    rng_energy = np.random.default_rng(
        args.energy_smear_seed
    )
    common_standard_normal = rng_energy.normal(
        size=mixed_events["kinetic"].shape
    )

    full_file_counts = np.asarray(
        [
            np.count_nonzero(
                reaction_channels == channel_id
            )
            for channel_id in CHANNEL_IDS
        ],
        dtype=float,
    )
    if np.sum(full_file_counts) <= 0.0:
        raise RuntimeError(
            "No recognized 162 channel labels in mixed reaction tree"
        )
    flat_truth = (
        full_file_counts / np.sum(full_file_counts)
    )
    plot_truth = np.repeat(
        flat_truth[None, :],
        fwhm_values.size,
        axis=0,
    )

    closure_angle_source = (
        resolve_closure_angle_source(args)
    )
    print(
        "Analysis settings:"
    )
    print(
        f"  full angular coverage"
    )
    print(
        f"  Delta-theta={args.theta_granularity_deg:g} deg"
    )
    print(
        f"  Delta-phi={args.phi_granularity_deg:g} deg"
    )
    print(
        f"  threshold={args.energy_threshold_keV:g} keV"
    )
    print(
        f"  closure angle source={closure_angle_source}"
    )
    print(
        f"  N_fit={args.n_fit:,}"
    )
    print(
        "  flat truth reference="
        + np.array2string(
            flat_truth,
            precision=8,
            separator=", ",
        )
    )

    def build_ideal_template_histograms() -> np.ndarray:
        histograms = np.empty(
            (3, args.bins, args.bins),
            dtype=float,
        )
        for channel, events in enumerate(pure_events):
            reconstructed, _ = reconstruct_with_granularity(
                events,
                0.0,
                0.0,
                args.theta_offset_deg,
                args.phi_offset_deg,
                angle_source="momentum",
                sphere_center_mm=sphere_center_mm,
            )
            histograms[channel] = dalitz_histogram(
                reconstructed,
                args.bins,
                limits,
            )
        return histograms

    ideal_template_histograms = (
        build_ideal_template_histograms()
    )
    global_fractions = simplex_grid(
        args.grid_size
    )

    reconstructed_closure, _ = reconstruct_with_granularity(
        mixed_events,
        args.theta_granularity_deg,
        args.phi_granularity_deg,
        args.theta_offset_deg,
        args.phi_offset_deg,
        angle_source=closure_angle_source,
        sphere_center_mm=sphere_center_mm,
    )

    n_points = fwhm_values.size
    medians = np.empty((n_points, 3), dtype=float)
    lower95 = np.empty((n_points, 3), dtype=float)
    upper95 = np.empty((n_points, 3), dtype=float)

    for index, fwhm_kev in enumerate(
        fwhm_values
    ):
        print(
            "\n"
            f"[{index + 1}/{n_points}] "
            f"Si energy FWHM={fwhm_kev:g} keV"
        )

        response_events, _, diagnostics = (
            apply_closure_energy_response(
                reconstructed_closure,
                float(fwhm_kev),
                args.energy_threshold_keV,
                common_standard_normal,
            )
        )
        selected, generated_count = select_exact_nfit(
            response_events,
            shuffled_reaction_ids,
            args.n_fit,
        )
        mixed_histogram = dalitz_histogram(
            selected,
            args.bins,
            limits,
        )

        prefix_counts, prefix_truth = generated_prefix_truth(
            shuffled_reaction_ids,
            generated_count,
            sorted_reaction_ids,
            sorted_reaction_channels,
        )

        fit = fit_histogram(
            mixed_histogram,
            ideal_template_histograms,
            efficiencies,
            global_fractions,
            args.grid_size,
            args.local_grid_size,
            args.local_sigma_span,
            args.template_pseudocount,
        )

        median = np.asarray(
            fit["median"]
        )
        low = np.asarray(
            fit["lower95"]
        )
        high = np.asarray(
            fit["upper95"]
        )
        medians[index] = median
        lower95[index] = low
        upper95[index] = high

        print(
            "  generated-prefix truth="
            + np.array2string(
                prefix_truth,
                precision=7,
                separator=", ",
            )
            + "  counts="
            + np.array2string(
                prefix_counts,
                separator=", ",
            )
        )
        print(
            "  median="
            + np.array2string(
                median,
                precision=7,
                separator=", ",
            )
        )
        print(
            "  95% widths="
            + np.array2string(
                high - low,
                precision=7,
                separator=", ",
            )
        )
        print(
            f"  N_generated={generated_count:,}, "
            f"threshold survival="
            f"{diagnostics['threshold_survival']:.6%}"
        )
        print(
            f"  fitted-sum={np.sum(median):.12f}"
        )

    fixed_ylims = None
    if not args.auto_y_range:
        fixed_ylims = (
            tuple(map(float, args.alpha0_ylim)),
            tuple(map(float, args.alpha1_ylim)),
            tuple(map(float, args.direct_ylim)),
        )

    figure = plot_fraction_scan(
        fwhm_values,
        medians,
        lower95,
        upper95,
        plot_truth,
        "Si energy resolution FWHM (keV)",
        fixed_ylims=fixed_ylims,
    )

    output_path = Path(outdir) / Path(args.output).name
    if save:
        Path(outdir).mkdir(parents=True, exist_ok=True)
        figure.savefig(output_path, dpi=300, bbox_inches="tight")
        print("\nWrote:")
        print(f"  {output_path}")

    return {
        "figure": figure,
        "fwhm_values": fwhm_values,
        "medians": medians,
        "lower95": lower95,
        "upper95": upper95,
        "truth": plot_truth,
    }


def _default_args() -> argparse.Namespace:
    """Argparse defaults with placeholder positionals (for code calls)."""
    return parse_args(["_", "_", "_", "_"])


def scan_si_energy_resolution_162_channel_fit(
    seq_be_gs, seq_be2plus, direct, mixed,
    save=False, outdir=".", **overrides,
):
    """Si energy-resolution response scan of the 162-keV channel-fraction fit.

    Parameters
    ----------
    seq_be_gs, seq_be2plus, direct : str or pathlib.Path
        The three pure single-channel benchmark ROOT files (templates).
    mixed : str or pathlib.Path
        A single mixed pseudodata ROOT file. Different Si energy-resolution
        FWHM values are applied internally (see fwhm_values_keV), so only one
        mixed file is needed.
    save : bool, optional
        False (default) -> do not write a file.
        True            -> write the scan figure into outdir.
    outdir : str or pathlib.Path, optional
        Directory for the saved PNG (created if needed; default ".").
    **overrides :
        Any parser option, e.g. fwhm_values_keV=[0, 50, 100, 200],
        n_fit=100000, bins=60, output="my.png", auto_y_range=True.

    Returns
    -------
    dict with keys: figure, fwhm_values, medians, lower95, upper95, truth.
    """
    args = _default_args()
    args.seq_be_gs = Path(seq_be_gs)
    args.seq_be2plus = Path(seq_be2plus)
    args.direct = Path(direct)
    args.mixed = Path(mixed)
    for key, value in overrides.items():
        if not hasattr(args, key):
            raise TypeError(f"unknown option: {key!r}")
        setattr(args, key, value)
    return _run(args, save=save, outdir=outdir)


if __name__ == "__main__":
    args = parse_args()
    _run(args, save=True, outdir=".")
