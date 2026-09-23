#!/usr/bin/env python3
"""Bayesian closure fit of the three 162-keV three-alpha channels.

The four positional inputs are, in this exact order:

1. pure 162seq_Be_gs benchmark ROOT file;
2. pure 162seq_Be2plus benchmark ROOT file;
3. pure 162direct benchmark ROOT file;
4. mixed ROOT file treated as pseudodata.

Only complete three-alpha kinematics from ``virtual_sphere`` enter the fit.
The only event-quality requirements are three unique alpha tracks, finite
four-vectors, positive kinetic energies and per-alpha mass-shell consistency.
Energy closure, missing momentum, invariant-energy regions and channel labels
are not used as selection cuts.
The mixed-file reaction-channel labels are read only after the fit to report
the closure truth.  A uniform Dirichlet(1,1,1) prior is evaluated directly on
the two-dimensional fraction simplex; MCMC is neither needed nor used.

Outputs
-------
bayes_162_pull.png
bayes_162_posterior_contour.png

No JSON file or Python bytecode cache is written.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np

import uproot


ALPHA_PDG = 1000020040
ALPHA_MASS_MEV = 3727.379378
SOURCE_H11B_REACTION_PRODUCT = 1

CHANNEL_IDS = np.array([1, 2, 3], dtype=np.int32)
CHANNEL_KEYS = ("162seq_Be_gs", "162seq_Be2plus", "162direct")
CHANNEL_LABELS = (
    r"162 seq. $^{8}$Be(g.s.)",
    r"162 seq. $^{8}$Be($2^{+}$)",
    "162 direct",
)
CHANNEL_COLORS = ("#0072B2", "#D55E00", "#009E73")

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
REACTION_BRANCHES = ("event", "h11b_reaction_channel")


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.labelsize": 14,
            "axes.titlesize": 14,
            "legend.fontsize": 11,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "axes.linewidth": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )


def tree_names(root_file) -> set[str]:
    return {str(key).split(";")[0] for key in root_file.keys()}


def require_tree(root_file, name: str):
    if name not in tree_names(root_file):
        raise RuntimeError(
            f"Missing ROOT tree '{name}'. Available objects: "
            + ", ".join(sorted(tree_names(root_file)))
        )
    return root_file[name]


def require_branches(tree, required: tuple[str, ...], tree_name: str) -> None:
    available = {str(key).split(";")[0] for key in tree.keys()}
    missing = [name for name in required if name not in available]
    if missing:
        raise RuntimeError(
            f"Tree '{tree_name}' is missing branches: {', '.join(missing)}"
        )


def load_reaction_channels(root_path: Path) -> np.ndarray:
    with uproot.open(root_path) as root_file:
        tree = require_tree(root_file, "reaction")
        require_branches(tree, REACTION_BRANCHES, "reaction")
        arrays = tree.arrays(REACTION_BRANCHES, library="np")

    event = np.asarray(arrays["event"], dtype=np.int64)
    channel = np.asarray(arrays["h11b_reaction_channel"], dtype=np.int32)
    if event.size == 0:
        raise RuntimeError(f"Reaction tree is empty in {root_path}")
    if np.unique(event).size != event.size:
        raise RuntimeError(f"Reaction event IDs are not unique in {root_path}")
    return channel


def load_complete_three_alpha(root_path: Path) -> dict[str, np.ndarray | int]:
    """Load complete sphere events without reading their channel labels."""
    with uproot.open(root_path) as root_file:
        tree = require_tree(root_file, "virtual_sphere")
        require_branches(tree, SPHERE_BRANCHES, "virtual_sphere")
        data = tree.arrays(SPHERE_BRANCHES, library="np")

    selected = (
        (np.asarray(data["pdg"]) == ALPHA_PDG)
        & (
            np.asarray(data["h11b_particle_source"])
            == SOURCE_H11B_REACTION_PRODUCT
        )
        & (np.asarray(data["generator_particle_index"]) >= 0)
        & (np.asarray(data["generator_particle_index"]) <= 2)
    )
    if not np.any(selected):
        raise RuntimeError(
            f"No labeled H11B reaction-product alpha crossings in {root_path}"
        )

    event = np.asarray(data["event_id"][selected], dtype=np.int64)
    index = np.asarray(data["generator_particle_index"][selected], dtype=np.int32)
    kinetic = np.asarray(data["kinetic_energy_MeV"][selected], dtype=float)
    momentum = np.column_stack(
        [
            np.asarray(data["px_MeV_c"][selected], dtype=float),
            np.asarray(data["py_MeV_c"][selected], dtype=float),
            np.asarray(data["pz_MeV_c"][selected], dtype=float),
        ]
    )

    order = np.lexsort((index, event))
    event = event[order]
    index = index[order]
    kinetic = kinetic[order]
    momentum = momentum[order]

    unique_event, first, multiplicity = np.unique(
        event, return_index=True, return_counts=True
    )
    candidate_first = first[multiplicity == 3]
    candidate_event = unique_event[multiplicity == 3]
    if candidate_event.size == 0:
        raise RuntimeError(f"No complete three-alpha events in {root_path}")

    rows = candidate_first[:, None] + np.arange(3, dtype=np.int64)[None, :]
    correct_indices = np.all(index[rows] == np.array([0, 1, 2]), axis=1)
    rows = rows[correct_indices]
    complete_event = candidate_event[correct_indices]
    if complete_event.size == 0:
        raise RuntimeError(
            f"No event has unique generator indices 0, 1 and 2 in {root_path}"
        )

    return {
        "event_id": complete_event,
        "kinetic": kinetic[rows],
        "momentum": momentum[rows],
        "n_complete": int(complete_event.size),
        "multiplicity_counts": {
            int(value): int(number)
            for value, number in zip(*np.unique(multiplicity, return_counts=True))
        },
    }


def event_quality(
    kinetic: np.ndarray,
    momentum: np.ndarray,
) -> dict[str, np.ndarray]:
    """Calculate only channel-independent per-alpha validity quantities."""
    finite = np.all(np.isfinite(kinetic), axis=1) & np.all(
        np.isfinite(momentum), axis=(1, 2)
    )
    positive = np.all(kinetic >= 0.0, axis=1)

    momentum2 = np.sum(momentum**2, axis=2)
    expected_momentum2 = kinetic * (kinetic + 2.0 * ALPHA_MASS_MEV)
    mass_shell_residual = np.max(
        np.abs(momentum2 - expected_momentum2)
        / np.maximum(expected_momentum2, 1.0),
        axis=1,
    )

    return {
        "finite": finite,
        "positive": positive,
        "mass_shell": mass_shell_residual,
    }


def apply_quality_cuts(
    events: dict[str, np.ndarray | int | dict[int, int]],
    args: argparse.Namespace,
) -> tuple[dict[str, np.ndarray | int | dict[int, int]], list[tuple[str, int]]]:
    """Apply exactly the same channel-blind cuts to templates and pseudodata."""
    kinetic = np.asarray(events["kinetic"], dtype=float)
    momentum = np.asarray(events["momentum"], dtype=float)
    quality = event_quality(kinetic, momentum)
    stages = (
        ("finite four-vectors", quality["finite"]),
        ("non-negative alpha kinetic energies", quality["positive"]),
        ("alpha mass shell", quality["mass_shell"] <= args.mass_shell_tolerance),
    )
    accepted = np.ones(kinetic.shape[0], dtype=bool)
    flow: list[tuple[str, int]] = [("complete three-alpha input", int(accepted.size))]
    for name, stage_mask in stages:
        accepted &= stage_mask
        flow.append((name, int(np.count_nonzero(accepted))))

    if not np.any(accepted):
        raise RuntimeError("All complete three-alpha events were rejected by quality cuts")
    filtered = {
        "event_id": np.asarray(events["event_id"])[accepted],
        "kinetic": kinetic[accepted],
        "momentum": momentum[accepted],
        "n_complete": int(events["n_complete"]),
        "n_accepted": int(np.count_nonzero(accepted)),
        "multiplicity_counts": events["multiplicity_counts"],
        "quality": quality,
    }
    return filtered, flow


def three_alpha_cm_kinetic(
    kinetic: np.ndarray, momentum: np.ndarray
) -> np.ndarray:
    energy = kinetic + ALPHA_MASS_MEV
    total_energy = np.sum(energy, axis=1)
    total_momentum = np.sum(momentum, axis=1)
    beta = total_momentum / total_energy[:, None]
    beta2 = np.sum(beta**2, axis=1)
    gamma = 1.0 / np.sqrt(np.clip(1.0 - beta2, 1.0e-15, None))

    beta_dot_p = np.sum(momentum * beta[:, None, :], axis=2)
    energy_cm = gamma[:, None] * (energy - beta_dot_p)
    return energy_cm - ALPHA_MASS_MEV


def ordered_dalitz(energy_cm: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return one statistically independent, energy-ordered point per event."""
    ordered = np.sort(energy_cm, axis=1)[:, ::-1]
    high, middle, low = ordered.T
    energy_sum = high + middle + low
    valid = np.isfinite(energy_sum) & (energy_sum > 0.0)
    x = np.sqrt(3.0) * (middle[valid] - low[valid]) / energy_sum[valid]
    y = (2.0 * high[valid] - middle[valid] - low[valid]) / energy_sum[valid]
    return x, y


def dalitz_histogram(
    events: dict[str, np.ndarray | int], bins: int, limits: tuple[float, float]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    energy_cm = three_alpha_cm_kinetic(events["kinetic"], events["momentum"])
    x, y = ordered_dalitz(energy_cm)
    histogram, x_edges, y_edges = np.histogram2d(
        x, y, bins=bins, range=[limits, limits]
    )
    return histogram.astype(float), x_edges, y_edges


def pure_template(
    root_path: Path,
    expected_channel: int,
    bins: int,
    limits: tuple[float, float],
    args: argparse.Namespace,
) -> dict[str, object]:
    reaction_channel = load_reaction_channels(root_path)
    unique, nums = np.unique(reaction_channel, return_counts=True)
    channel_counts = {int(k): int(v) for k, v in zip(unique, nums)}
    if not np.all(reaction_channel == expected_channel):
        raise RuntimeError(
            f"{root_path} is not a pure channel-{expected_channel} sample: "
            f"counts={channel_counts}"
        )

    complete_events = load_complete_three_alpha(root_path)
    events, cut_flow = apply_quality_cuts(complete_events, args)
    histogram, x_edges, y_edges = dalitz_histogram(events, bins, limits)
    efficiency = events["n_accepted"] / reaction_channel.size
    return {
        "histogram": histogram,
        "x_edges": x_edges,
        "y_edges": y_edges,
        "efficiency": float(efficiency),
        "n_reaction": int(reaction_channel.size),
        "n_complete": int(complete_events["n_complete"]),
        "n_accepted": int(events["n_accepted"]),
        "cut_flow": cut_flow,
        "multiplicity_counts": complete_events["multiplicity_counts"],
    }


def simplex_grid(grid_size: int) -> np.ndarray:
    values = np.linspace(0.0, 1.0, grid_size)
    f1, f2 = np.meshgrid(values, values, indexing="ij")
    valid = f1 + f2 <= 1.0 + 1.0e-12
    f1 = f1[valid]
    f2 = f2[valid]
    f3 = 1.0 - f1 - f2
    f3 = np.maximum(f3, 0.0)
    return np.column_stack([f1, f2, f3])


def evaluate_posterior(
    fractions: np.ndarray,
    template_probability: np.ndarray,
    efficiencies: np.ndarray,
    data_counts: np.ndarray,
    batch_size: int = 512,
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate uniform-prior posterior on the generator-fraction simplex."""
    log_likelihood = np.empty(fractions.shape[0], dtype=float)
    for start in range(0, fractions.shape[0], batch_size):
        stop = min(start + batch_size, fractions.shape[0])
        generator_fraction = fractions[start:stop]
        selected_weight = generator_fraction * efficiencies[None, :]
        selected_weight /= np.sum(selected_weight, axis=1, keepdims=True)
        bin_probability = selected_weight @ template_probability
        log_likelihood[start:stop] = np.log(bin_probability) @ data_counts

    shifted = log_likelihood - np.max(log_likelihood)
    posterior = np.exp(shifted)
    posterior /= np.sum(posterior)
    return posterior, log_likelihood


def refine_map_fraction(
    initial: np.ndarray,
    global_step: float,
    template_probability: np.ndarray,
    efficiencies: np.ndarray,
    data_counts: np.ndarray,
    refinements: int = 3,
    grid_size: int = 81,
) -> np.ndarray:
    """Iteratively refine the maximum inside the two-parameter simplex."""
    center = np.asarray(initial, dtype=float).copy()
    half_width = 2.0 * global_step
    for _ in range(refinements):
        f1_axis = np.linspace(
            max(0.0, center[0] - half_width),
            min(1.0, center[0] + half_width),
            grid_size,
        )
        f2_axis = np.linspace(
            max(0.0, center[1] - half_width),
            min(1.0, center[1] + half_width),
            grid_size,
        )
        f1_mesh, f2_mesh = np.meshgrid(f1_axis, f2_axis, indexing="ij")
        valid = f1_mesh + f2_mesh <= 1.0 + 1.0e-14
        trial = np.column_stack(
            [
                f1_mesh[valid],
                f2_mesh[valid],
                1.0 - f1_mesh[valid] - f2_mesh[valid],
            ]
        )
        _, log_likelihood = evaluate_posterior(
            trial, template_probability, efficiencies, data_counts
        )
        center = trial[int(np.argmax(log_likelihood))]
        half_width /= 8.0
    return center


def estimate_local_covariance(
    mode: np.ndarray,
    global_step: float,
    template_probability: np.ndarray,
    efficiencies: np.ndarray,
    data_counts: np.ndarray,
) -> tuple[np.ndarray | None, float]:
    """Estimate the two-parameter inverse Hessian around an interior mode."""
    h = max(global_step / 100.0, 1.0e-6)
    f1, f2 = float(mode[0]), float(mode[1])
    if min(f1, f2, 1.0 - f1 - f2) <= 2.5 * h:
        return None, h

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
        ]
    )
    points = np.column_stack(
        [points2, 1.0 - points2[:, 0] - points2[:, 1]]
    )
    _, log_likelihood = evaluate_posterior(
        points, template_probability, efficiencies, data_counts
    )
    nll = -log_likelihood
    h11 = (nll[1] - 2.0 * nll[0] + nll[2]) / h**2
    h22 = (nll[3] - 2.0 * nll[0] + nll[4]) / h**2
    h12 = (nll[5] - nll[6] - nll[7] + nll[8]) / (4.0 * h**2)
    hessian = np.array([[h11, h12], [h12, h22]], dtype=float)
    eigenvalues = np.linalg.eigvalsh(hessian)
    if not np.all(np.isfinite(hessian)) or np.any(eigenvalues <= 0.0):
        return None, h
    return np.linalg.inv(hessian), h


def adaptive_local_grid(
    mode: np.ndarray,
    covariance2: np.ndarray | None,
    global_step: float,
    sigma_span: float,
    grid_size: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Create a fine rectangular grid around the posterior mode."""
    if covariance2 is None:
        half1 = half2 = 2.0 * global_step
    else:
        sigma = np.sqrt(np.maximum(np.diag(covariance2), 0.0))
        half1 = max(sigma_span * sigma[0], global_step / 5.0)
        half2 = max(sigma_span * sigma[1], global_step / 5.0)

    f1_axis = np.linspace(max(0.0, mode[0] - half1), min(1.0, mode[0] + half1), grid_size)
    f2_axis = np.linspace(max(0.0, mode[1] - half2), min(1.0, mode[1] + half2), grid_size)
    f1_mesh, f2_mesh = np.meshgrid(f1_axis, f2_axis, indexing="ij")
    valid = f1_mesh + f2_mesh <= 1.0 + 1.0e-14
    fractions = np.column_stack(
        [
            f1_mesh[valid],
            f2_mesh[valid],
            1.0 - f1_mesh[valid] - f2_mesh[valid],
        ]
    )
    return fractions, f1_axis, f2_axis, f1_mesh, valid


def covariance_and_correlation(
    covariance: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    standard_deviation = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    denominator = np.outer(standard_deviation, standard_deviation)
    correlation = np.divide(
        covariance,
        denominator,
        out=np.zeros_like(covariance),
        where=denominator > 0.0,
    )
    return standard_deviation, correlation


def hpd_density_threshold(posterior: np.ndarray, probability: float) -> float:
    """Posterior-height threshold enclosing the requested HPD probability."""
    order = np.argsort(posterior)[::-1]
    cumulative = np.cumsum(posterior[order])
    index = min(int(np.searchsorted(cumulative, probability)), posterior.size - 1)
    return float(posterior[order[index]])


def weighted_quantile(
    values: np.ndarray, weights: np.ndarray, probabilities: tuple[float, ...]
) -> np.ndarray:
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cumulative = np.cumsum(weights)
    cumulative /= cumulative[-1]
    return np.interp(np.asarray(probabilities), cumulative, values)


def posterior_summary(
    fractions: np.ndarray, posterior: np.ndarray, log_likelihood: np.ndarray
) -> dict[str, np.ndarray]:
    mean = np.sum(fractions * posterior[:, None], axis=0)
    map_fraction = fractions[int(np.argmax(log_likelihood))]
    lower = np.empty(3)
    median = np.empty(3)
    upper = np.empty(3)
    lower95 = np.empty(3)
    upper95 = np.empty(3)
    for channel in range(3):
        (
            lower95[channel],
            lower[channel],
            median[channel],
            upper[channel],
            upper95[channel],
        ) = weighted_quantile(
            fractions[:, channel], posterior, (0.025, 0.16, 0.50, 0.84, 0.975)
        )
    centered = fractions - mean[None, :]
    covariance = (centered * posterior[:, None]).T @ centered
    return {
        "mean": mean,
        "map": map_fraction,
        "median": median,
        "lower": lower,
        "upper": upper,
        "lower95": lower95,
        "upper95": upper95,
        "covariance": covariance,
    }


def selected_bin_probability(
    generator_fraction: np.ndarray,
    template_probability: np.ndarray,
    efficiencies: np.ndarray,
) -> np.ndarray:
    weight = generator_fraction * efficiencies
    weight /= np.sum(weight)
    return weight @ template_probability


def posterior_predictive_p_value(
    fractions: np.ndarray,
    posterior: np.ndarray,
    template_probability: np.ndarray,
    efficiencies: np.ndarray,
    data_counts: np.ndarray,
    draws: int,
    seed: int,
) -> float:
    if draws <= 0:
        return float("nan")
    rng = np.random.default_rng(seed)
    chosen = rng.choice(fractions.shape[0], size=draws, p=posterior)
    n_total = int(np.sum(data_counts))
    more_extreme = 0
    for index in chosen:
        probability = selected_bin_probability(
            fractions[index], template_probability, efficiencies
        )
        expected = n_total * probability
        denominator = np.maximum(expected, 1.0)
        observed_statistic = np.sum((data_counts - expected) ** 2 / denominator)
        replicated = rng.multinomial(n_total, probability)
        replicated_statistic = np.sum((replicated - expected) ** 2 / denominator)
        more_extreme += replicated_statistic >= observed_statistic
    return more_extreme / draws


def truth_from_mixed_reaction(root_path: Path) -> tuple[np.ndarray, dict[int, int]]:
    """Read truth only for the final closure comparison, never for fitting."""
    channel = load_reaction_channels(root_path)
    valid = np.isin(channel, CHANNEL_IDS)
    if not np.all(valid):
        bad = channel[~valid]
        unique, nums = np.unique(bad, return_counts=True)
        raise RuntimeError(
            "Mixed reaction tree contains channels outside 1,2,3: "
            + str({int(k): int(v) for k, v in zip(unique, nums)})
        )
    channel_counts = np.array(
        [np.count_nonzero(channel == channel_id) for channel_id in CHANNEL_IDS],
        dtype=float,
    )
    truth = channel_counts / np.sum(channel_counts)
    return truth, {
        int(channel_id): int(count)
        for channel_id, count in zip(CHANNEL_IDS, channel_counts)
    }


def plot_pull(
    data_histogram: np.ndarray,
    model_histogram: np.ndarray,
    active_mask: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
):
    pull = np.full(data_histogram.shape, np.nan)
    valid = active_mask & (model_histogram >= 1.0)
    pull[valid] = (
        data_histogram[valid] - model_histogram[valid]
    ) / np.sqrt(model_histogram[valid])

    fig, axis = plt.subplots(figsize=(6.5, 5.5))
    pull_image = axis.pcolormesh(
        x_edges,
        y_edges,
        np.ma.masked_invalid(pull.T),
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-5.0, vcenter=0.0, vmax=5.0),
        shading="auto",
    )
    axis.set_title("Pull")
    axis.set_xlabel(r"$\sqrt{3}(E_M-E_L)/\sum E_i$")
    axis.set_ylabel(r"$(2E_H-E_M-E_L)/\sum E_i$")
    axis.set_aspect("equal", adjustable="box")
    pull_bar = fig.colorbar(pull_image, ax=axis, pad=0.02)
    pull_bar.set_label(r"$(n-\mu)/\sqrt{\mu}$")
    fig.tight_layout()
    return fig


def plot_posterior_contour(
    f1_axis: np.ndarray,
    f2_axis: np.ndarray,
    valid: np.ndarray,
    posterior: np.ndarray,
    mode: np.ndarray,
):
    f1_mesh, f2_mesh = np.meshgrid(f1_axis, f2_axis, indexing="ij")
    posterior_matrix = np.full(valid.shape, np.nan)
    posterior_matrix[valid] = posterior
    relative = posterior_matrix / np.nanmax(posterior_matrix)
    threshold68 = hpd_density_threshold(posterior, 0.68) / np.max(posterior)
    threshold95 = hpd_density_threshold(posterior, 0.95) / np.max(posterior)

    fig, axis = plt.subplots(figsize=(6.6, 5.5))
    image = axis.pcolormesh(
        f1_mesh,
        f2_mesh,
        np.ma.masked_invalid(relative),
        cmap="viridis",
        shading="auto",
        vmin=0.0,
        vmax=1.0,
    )
    contour95 = axis.contour(
        f1_mesh,
        f2_mesh,
        relative,
        levels=[threshold95],
        colors=["white"],
        linewidths=1.6,
    )
    contour68 = axis.contour(
        f1_mesh,
        f2_mesh,
        relative,
        levels=[threshold68],
        colors=["#D55E00"],
        linewidths=2.0,
    )
    axis.clabel(contour95, fmt={threshold95: "95%"}, inline=True, fontsize=10)
    axis.clabel(contour68, fmt={threshold68: "68%"}, inline=True, fontsize=10)
    axis.scatter(
        mode[0], mode[1], marker="*", s=115, color="#D55E00", zorder=4
    )
    axis.annotate(
        rf"$({mode[0]:.5f},\ {mode[1]:.5f})$",
        xy=(mode[0], mode[1]),
        xytext=(10, -20),
        textcoords="offset points",
        color="#D55E00",
        fontsize=11,
        ha="left",
        va="bottom",
        bbox={"boxstyle": "round,pad=0.22", "fc": "white", "ec": "none", "alpha": 0.82},
    )
    axis.set_xlabel(r"$f_{\mathrm{g.s.}}$")
    axis.set_ylabel(r"$f_{2^+}$")
    axis.set_title("Joint posterior for channel fractions")
    fig.colorbar(image, ax=axis, pad=0.02, label="Posterior / posterior maximum")
    fig.tight_layout()
    return fig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Bayesian Dalitz-template closure fit for three 162-keV channels."
    )
    parser.add_argument("seq_be_gs", type=Path, help="pure 162seq_Be_gs ROOT file")
    parser.add_argument(
        "seq_be2plus", type=Path, help="pure 162seq_Be2plus ROOT file"
    )
    parser.add_argument("direct", type=Path, help="pure 162direct ROOT file")
    parser.add_argument("mixed", type=Path, help="mixed pseudodata ROOT file")
    parser.add_argument("--mass-shell-tolerance", type=float, default=1.0e-4)
    parser.add_argument("--bins", type=int, default=60, help="bins per Dalitz axis")
    parser.add_argument(
        "--grid-size",
        type=int,
        default=201,
        help="global fraction-grid points along each simplex edge",
    )
    parser.add_argument(
        "--local-grid-size",
        type=int,
        default=401,
        help="adaptive local posterior-grid points per axis",
    )
    parser.add_argument(
        "--local-sigma-span",
        type=float,
        default=6.0,
        help="local-grid half-width in estimated marginal sigmas",
    )
    parser.add_argument(
        "--template-pseudocount",
        type=float,
        default=0.5,
        help="finite-template pseudocount per active bin",
    )
    parser.add_argument(
        "--ppc-draws",
        type=int,
        default=300,
        help="posterior-predictive replicated datasets",
    )
    parser.add_argument("--seed", type=int, default=162, help="PPC random seed")
    parser.add_argument(
        "--prefix", default="bayes_162", help="output image filename prefix"
    )
    return parser.parse_args()


def _run(args, save=False, outdir="."):
    if args.bins < 10:
        raise ValueError("bins must be at least 10")
    if args.grid_size < 21:
        raise ValueError("grid_size must be at least 21")
    if args.local_grid_size < 101:
        raise ValueError("local_grid_size must be at least 101")
    if args.local_sigma_span <= 2.0:
        raise ValueError("local_sigma_span must be greater than 2")
    if args.template_pseudocount <= 0.0:
        raise ValueError("template_pseudocount must be positive")
    for path in (args.seq_be_gs, args.seq_be2plus, args.direct, args.mixed):
        if not Path(path).is_file():
            raise FileNotFoundError(f"Input file not found: {path}")

    configure_style()
    limits = (-1.05, 1.05)
    pure_paths = (args.seq_be_gs, args.seq_be2plus, args.direct)

    print("Loading pure-channel templates...")
    templates = [
        pure_template(path, int(channel_id), args.bins, limits, args)
        for path, channel_id in zip(pure_paths, CHANNEL_IDS)
    ]
    raw_template_histograms = np.stack(
        [item["histogram"] for item in templates], axis=0
    )
    x_edges = templates[0]["x_edges"]
    y_edges = templates[0]["y_edges"]
    efficiencies = np.array([item["efficiency"] for item in templates])

    print("Loading mixed pseudodata without channel labels...")
    mixed_complete_events = load_complete_three_alpha(args.mixed)
    mixed_events, mixed_cut_flow = apply_quality_cuts(mixed_complete_events, args)
    data_histogram, data_x_edges, data_y_edges = dalitz_histogram(
        mixed_events, args.bins, limits
    )
    if not np.array_equal(x_edges, data_x_edges) or not np.array_equal(
        y_edges, data_y_edges
    ):
        raise RuntimeError("Internal Dalitz bin-edge mismatch")

    active_mask = (
        np.any(raw_template_histograms > 0.0, axis=0) | (data_histogram > 0.0)
    )
    data_counts = data_histogram[active_mask].astype(np.int64)
    template_counts = raw_template_histograms[:, active_mask]
    template_probability = template_counts + args.template_pseudocount
    template_probability /= np.sum(template_probability, axis=1, keepdims=True)

    global_fractions = simplex_grid(args.grid_size)
    print(
        f"Evaluating global posterior at {global_fractions.shape[0]:,} points "
        "with Dirichlet(1,1,1) prior..."
    )
    global_posterior, global_log_likelihood = evaluate_posterior(
        global_fractions,
        template_probability,
        efficiencies,
        data_counts,
    )
    global_map = global_fractions[int(np.argmax(global_log_likelihood))]
    global_step = 1.0 / (args.grid_size - 1)
    refined_map = refine_map_fraction(
        global_map,
        global_step,
        template_probability,
        efficiencies,
        data_counts,
    )
    covariance2, hessian_step = estimate_local_covariance(
        refined_map,
        global_step,
        template_probability,
        efficiencies,
        data_counts,
    )
    local_fractions, f1_axis, f2_axis, _, local_valid = adaptive_local_grid(
        refined_map,
        covariance2,
        global_step,
        args.local_sigma_span,
        args.local_grid_size,
    )
    print(
        f"Evaluating adaptive local posterior at {local_fractions.shape[0]:,} points..."
    )
    local_posterior, local_log_likelihood = evaluate_posterior(
        local_fractions,
        template_probability,
        efficiencies,
        data_counts,
    )
    summary = posterior_summary(
        local_fractions, local_posterior, local_log_likelihood
    )
    posterior_std, posterior_correlation = covariance_and_correlation(
        summary["covariance"]
    )

    posterior_bin_probability = selected_bin_probability(
        summary["mean"], template_probability, efficiencies
    )
    model_histogram = np.zeros_like(data_histogram, dtype=float)
    model_histogram[active_mask] = np.sum(data_counts) * posterior_bin_probability

    ppc_p_value = posterior_predictive_p_value(
        local_fractions,
        local_posterior,
        template_probability,
        efficiencies,
        data_counts,
        args.ppc_draws,
        args.seed,
    )

    # Reveal the truth only after all fit quantities have been calculated.
    truth, truth_counts = truth_from_mixed_reaction(args.mixed)
    _, truth_log_likelihood = evaluate_posterior(
        truth[None, :], template_probability, efficiencies, data_counts
    )
    truth_delta = max(
        0.0,
        float(-2.0 * (truth_log_likelihood[0] - np.max(local_log_likelihood))),
    )

    fig_pull = plot_pull(
        data_histogram, model_histogram, active_mask, x_edges, y_edges,
    )
    fig_contour = plot_posterior_contour(
        f1_axis, f2_axis, local_valid, local_posterior, summary["map"],
    )
    pull_output = Path(outdir) / f"{args.prefix}_pull.png"
    contour_output = Path(outdir) / f"{args.prefix}_posterior_contour.png"
    if save:
        Path(outdir).mkdir(parents=True, exist_ok=True)
        fig_pull.savefig(pull_output, dpi=300, bbox_inches="tight")
        fig_contour.savefig(contour_output, dpi=300, bbox_inches="tight")

    print("\nBayesian 162-keV channel-fraction closure")
    print("=" * 124)
    print(
        f"{'channel':<24} {'truth':>9} {'posterior median':>18} "
        f"{'68% credible interval':>25} {'95% credible interval':>25} {'MAP':>10}"
    )
    print("-" * 124)
    for i, key in enumerate(CHANNEL_KEYS):
        print(
            f"{key:<24} {truth[i]:9.5f} {summary['median'][i]:18.5f} "
            f"[{summary['lower'][i]:.5f}, {summary['upper'][i]:.5f}]"
            f" [{summary['lower95'][i]:.5f}, {summary['upper95'][i]:.5f}]"
            f" {summary['map'][i]:10.6f}"
        )
    print("-" * 124)
    print(
        "Quality cuts: exactly three unique alphas, finite four-vectors, "
        f"positive kinetic energies, mass-shell<={args.mass_shell_tolerance:g}"
    )
    print("Pure-template accepted-event efficiencies:")
    for key, efficiency, item in zip(CHANNEL_KEYS, efficiencies, templates):
        print(
            f"  {key:<22} {efficiency:.6%} "
            f"({item['n_accepted']}/{item['n_reaction']}); "
            f"complete before cuts={item['n_complete']}"
        )
    print("Mixed pseudodata quality-cut flow:")
    previous = mixed_cut_flow[0][1]
    for name, number in mixed_cut_flow:
        step_efficiency = number / previous if previous else 0.0
        print(f"  {name:<38} {number:>9}  step={step_efficiency:9.5%}")
        previous = number
    print(f"Mixed accepted events used in fit: {int(np.sum(data_counts)):,}")
    print(f"Mixed reaction truth counts: {truth_counts}")
    print(f"Truth -2DeltaLogPosterior: {truth_delta:.6g}")
    print(
        f"Posterior-predictive p-value ({args.ppc_draws} draws): "
        f"{ppc_p_value:.4f}"
    )
    print("Posterior standard deviations (templates treated as fixed):")
    for key, value in zip(CHANNEL_KEYS, posterior_std):
        print(f"  {key:<22} {value:.8f}")
    print("Posterior covariance matrix:")
    for row in summary["covariance"]:
        print("  " + " ".join(f"{value: .6e}" for value in row))
    print("Posterior correlation matrix:")
    for row in posterior_correlation:
        print("  " + " ".join(f"{value: .6f}" for value in row))
    local_step1 = float(f1_axis[1] - f1_axis[0])
    local_step2 = float(f2_axis[1] - f2_axis[0])
    bins_per_sigma = np.array(
        [
            posterior_std[0] / local_step1,
            posterior_std[1] / local_step2,
        ]
    )
    print(
        f"Grid diagnostics: global step={global_step:.6g}, "
        f"local steps=({local_step1:.6g}, {local_step2:.6g}), "
        f"bins/sigma=({bins_per_sigma[0]:.2f}, {bins_per_sigma[1]:.2f}), "
        f"Hessian step={hessian_step:.3g}"
    )
    if np.any(bins_per_sigma < 5.0):
        print("WARNING: local posterior resolution is below five bins per sigma; increase --local-grid-size.")
    print("Uniform prior: Dirichlet(1,1,1)")
    print("Inference: global simplex scan plus adaptive local posterior grid; no MCMC")
    print("Reported intervals condition on the finite pure templates as fixed shapes.")
    if save:
        print("\nWrote:")
        print(f"  {pull_output}")
        print(f"  {contour_output}")

    return {
        "fig_pull": fig_pull,
        "fig_contour": fig_contour,
        "summary": summary,
        "truth": truth,
        "ppc_p_value": ppc_p_value,
        "posterior_std": posterior_std,
    }


def bayes_fit_162_channels(
    seq_be_gs, seq_be2plus, direct, mixed,
    save=False, outdir=".",
    bins=60, grid_size=201, local_grid_size=401, local_sigma_span=6.0,
    template_pseudocount=0.5, ppc_draws=300, seed=162,
    mass_shell_tolerance=1.0e-4, prefix="bayes_162",
):
    """Bayesian Dalitz-template closure fit of the three 162-keV channels.

    Parameters
    ----------
    seq_be_gs, seq_be2plus, direct : str or pathlib.Path
        The three pure single-channel benchmark ROOT files (templates).
    mixed : str or pathlib.Path
        Mixed ROOT file treated as pseudodata.
    save : bool, optional
        False (default) -> do not write files.
        True            -> write <prefix>_pull.png and
                           <prefix>_posterior_contour.png into outdir.
    outdir : str or pathlib.Path, optional
        Directory for the two output PNGs (created if needed; default ".").
    bins, grid_size, local_grid_size, local_sigma_span, template_pseudocount,
    ppc_draws, seed, mass_shell_tolerance, prefix :
        Fit / grid / output controls (same meaning as the CLI flags).

    Returns
    -------
    dict with keys: fig_pull, fig_contour, summary, truth, ppc_p_value,
    posterior_std.
    """
    args = argparse.Namespace(
        seq_be_gs=Path(seq_be_gs), seq_be2plus=Path(seq_be2plus),
        direct=Path(direct), mixed=Path(mixed),
        mass_shell_tolerance=mass_shell_tolerance, bins=bins,
        grid_size=grid_size, local_grid_size=local_grid_size,
        local_sigma_span=local_sigma_span,
        template_pseudocount=template_pseudocount, ppc_draws=ppc_draws,
        seed=seed, prefix=prefix,
    )
    return _run(args, save=save, outdir=outdir)


if __name__ == "__main__":
    args = parse_args()
    _run(args, save=True, outdir=".")
