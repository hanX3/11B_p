#!/usr/bin/env python3
"""Plot comparisons for 675-keV alpha strict-model validation."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


FIELDS = [
    "reaction_channel",
    "e_3alpha_cm_alpha1",
    "e_3alpha_cm_alpha2",
    "e_3alpha_cm_alpha3",
    "eaa_8Be",
    "cos_chi_secondary_8be",
    "h11b675_strict_weight",
    "h11b675_strict_weight_max",
    "h11b675_strict_sampling_attempts",
]


def load_uproot():
    try:
        import uproot  # type: ignore
    except ImportError as exc:
        raise SystemExit("uproot is not installed; cannot plot ROOT files.") from exc
    return uproot


def load_pyroot():
    try:
        import ROOT  # type: ignore
    except ImportError as exc:
        raise SystemExit("PyROOT is not installed; cannot plot ROOT files with --backend pyroot.") from exc
    return ROOT


def load_dataset_uproot(path: Path) -> dict[str, np.ndarray] | None:
    uproot = load_uproot()
    with uproot.open(path) as root_file:
        if "tr" not in root_file:
            print(f"{path}: missing tree 'tr'; skipped")
            return None
        tree = root_file["tr"]
        branches = set(tree.keys())
        missing = [field for field in FIELDS if field not in branches]
        if missing:
            print(f"{path}: missing fields {', '.join(missing)}; skipped")
            return None
        data = {field: tree[field].array(library="np") for field in FIELDS}

    mask = data["reaction_channel"] == 1
    if not np.any(mask):
        print(f"{path}: no 675-alpha events; skipped")
        return None
    return {field: values[mask] for field, values in data.items() if field != "reaction_channel"}


def load_dataset_pyroot(path: Path) -> dict[str, np.ndarray] | None:
    ROOT = load_pyroot()
    root_file = ROOT.TFile.Open(str(path))
    if not root_file or root_file.IsZombie():
        print(f"{path}: could not open ROOT file; skipped")
        return None
    tree = root_file.Get("tr")
    if not tree:
        print(f"{path}: missing tree 'tr'; skipped")
        root_file.Close()
        return None
    branches = {branch.GetName() for branch in tree.GetListOfBranches()}
    missing = [field for field in FIELDS if field not in branches]
    if missing:
        print(f"{path}: missing fields {', '.join(missing)}; skipped")
        root_file.Close()
        return None

    values = {field: [] for field in FIELDS}
    for index in range(int(tree.GetEntries())):
        tree.GetEntry(index)
        if int(getattr(tree, "reaction_channel")) != 1:
            continue
        for field in FIELDS:
            if field == "reaction_channel":
                continue
            values[field].append(getattr(tree, field))
    root_file.Close()

    if not values["e_3alpha_cm_alpha1"]:
        print(f"{path}: no 675-alpha events; skipped")
        return None
    return {field: np.asarray(field_values) for field, field_values in values.items() if field != "reaction_channel"}


def load_dataset(path: Path, backend: str) -> dict[str, np.ndarray] | None:
    if backend == "uproot":
        return load_dataset_uproot(path)
    return load_dataset_pyroot(path)


def finite(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return values[np.isfinite(values)]


def single_alpha_energy(data: dict[str, np.ndarray]) -> np.ndarray:
    return finite(
        np.concatenate(
            [
                data["e_3alpha_cm_alpha1"],
                data["e_3alpha_cm_alpha2"],
                data["e_3alpha_cm_alpha3"],
            ]
        )
    )


def hist_l1(a: np.ndarray, b: np.ndarray, bins: int = 80) -> float:
    a = finite(a)
    b = finite(b)
    if a.size == 0 or b.size == 0:
        return float("nan")
    lo = float(min(np.min(a), np.min(b)))
    hi = float(max(np.max(a), np.max(b)))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return float("nan")
    ha, edges = np.histogram(a, bins=bins, range=(lo, hi), density=True)
    hb, _ = np.histogram(b, bins=edges, density=True)
    width = np.diff(edges)
    return float(np.sum(np.abs(ha - hb) * width))


def save_dalitz(label: str, data: dict[str, np.ndarray], outdir: Path) -> None:
    e1 = np.asarray(data["e_3alpha_cm_alpha1"], dtype=float)
    e2 = np.asarray(data["e_3alpha_cm_alpha2"], dtype=float)
    e3 = np.asarray(data["e_3alpha_cm_alpha3"], dtype=float)
    esum = e1 + e2 + e3
    mask = np.isfinite(esum) & (esum > 0)
    x = (e2[mask] - e3[mask]) / (np.sqrt(3.0) * esum[mask])
    y = e1[mask] / esum[mask] - 1.0 / 3.0

    plt.figure(figsize=(6, 5))
    plt.hist2d(x, y, bins=120, cmap="viridis")
    plt.xlabel("(E2 - E3) / (sqrt(3) Esum)")
    plt.ylabel("E1 / Esum - 1/3")
    plt.title(f"Dalitz {label}")
    plt.colorbar(label="entries")
    plt.tight_layout()
    plt.savefig(outdir / f"dalitz_{label}.png", dpi=180)
    plt.close()


def compare_hist(
    datasets: dict[str, dict[str, np.ndarray]],
    value_getter,
    outpath: Path,
    xlabel: str,
    bins: int = 100,
    logy: bool = False,
) -> None:
    values_by_label = {label: finite(value_getter(data)) for label, data in datasets.items()}
    nonempty = [values for values in values_by_label.values() if values.size]
    if not nonempty:
        print(f"{outpath.name}: no data; skipped")
        return
    lo = float(min(np.min(values) for values in nonempty))
    hi = float(max(np.max(values) for values in nonempty))
    if hi <= lo:
        hi = lo + 1.0

    plt.figure(figsize=(8, 5))
    for label, values in values_by_label.items():
        if values.size == 0:
            continue
        plt.hist(values, bins=bins, range=(lo, hi), histtype="step", density=True, label=label)
    plt.xlabel(xlabel)
    plt.ylabel("normalized entries")
    if logy:
        plt.yscale("log")
    plt.legend(fontsize="small")
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def strict_ratio(data: dict[str, np.ndarray]) -> np.ndarray:
    weight = np.asarray(data["h11b675_strict_weight"], dtype=float)
    weight_max = np.asarray(data["h11b675_strict_weight_max"], dtype=float)
    ratio = np.divide(weight, weight_max, out=np.full_like(weight, np.nan), where=weight_max > 0)
    return ratio


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy", type=Path)
    parser.add_argument("--strict", type=Path)
    parser.add_argument("--incoherent", type=Path)
    parser.add_argument("--nosym", type=Path)
    parser.add_argument("--l1only", type=Path)
    parser.add_argument("--l3only", type=Path)
    parser.add_argument("--phase0", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--backend", choices=("uproot", "pyroot"), default="uproot")
    args = parser.parse_args()

    paths = {
        "legacy": args.legacy,
        "strict": args.strict,
        "incoherent": args.incoherent,
        "nosym": args.nosym,
        "l1only": args.l1only,
        "l3only": args.l3only,
        "phase0": args.phase0,
    }

    datasets: dict[str, dict[str, np.ndarray]] = {}
    for label, path in paths.items():
        if path is None:
            continue
        data = load_dataset(path, args.backend)
        if data is not None:
            datasets[label] = data

    args.outdir.mkdir(parents=True, exist_ok=True)
    for label, data in datasets.items():
        save_dalitz(label, data, args.outdir)

    compare_hist(
        datasets,
        single_alpha_energy,
        args.outdir / "single_alpha_cm_energy_compare.png",
        "single-alpha CM kinetic energy",
    )
    compare_hist(
        datasets,
        lambda data: data["eaa_8Be"],
        args.outdir / "eaa_8be_compare.png",
        "eaa_8Be",
    )
    print("eaa_8Be is the stored 8Be-pair relative energy; all alpha pairings are not reconstructed by this script.")
    compare_hist(
        datasets,
        lambda data: data["cos_chi_secondary_8be"],
        args.outdir / "cos_chi_compare.png",
        "cos_chi_secondary_8be",
    )

    strict_datasets = {label: data for label, data in datasets.items() if label != "legacy"}
    compare_hist(
        strict_datasets,
        lambda data: data["h11b675_strict_sampling_attempts"],
        args.outdir / "strict_sampling_attempts.png",
        "strict sampling attempts",
        bins=80,
        logy=True,
    )
    compare_hist(
        strict_datasets,
        strict_ratio,
        args.outdir / "strict_weight_ratio.png",
        "strict weight / weight_max",
        bins=80,
        logy=True,
    )

    comparisons = [
        ("legacy", "strict"),
        ("strict", "incoherent"),
        ("strict", "nosym"),
        ("l1only", "l3only"),
        ("strict", "phase0"),
    ]
    for left, right in comparisons:
        if left not in datasets or right not in datasets:
            continue
        print(f"{left} vs {right}:")
        print(f"  single-alpha L1 distance = {hist_l1(single_alpha_energy(datasets[left]), single_alpha_energy(datasets[right])):.12g}")
        print(f"  eaa_8Be L1 distance = {hist_l1(datasets[left]['eaa_8Be'], datasets[right]['eaa_8Be']):.12g}")
        print(f"  cos_chi L1 distance = {hist_l1(datasets[left]['cos_chi_secondary_8be'], datasets[right]['cos_chi_secondary_8be']):.12g}")

    print(f"plots written to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
