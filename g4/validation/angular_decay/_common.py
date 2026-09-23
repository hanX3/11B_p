#!/usr/bin/env python3
"""Shared helpers for angular-decay validation scripts."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np

TREE_NAME = "tr"


class RootTree:
    def __init__(self, path: str | Path):
        import ROOT  # type: ignore

        ROOT.gROOT.SetBatch(True)
        self.ROOT = ROOT
        self.path = Path(path)
        self.file = ROOT.TFile.Open(str(self.path))
        if not self.file or self.file.IsZombie():
            raise SystemExit(f"{self.path}: could not open ROOT file")
        self.tree = self.file.Get(TREE_NAME)
        if not self.tree:
            raise SystemExit(f"{self.path}: missing tree '{TREE_NAME}'")
        self.num_entries = int(self.tree.GetEntries())
        self._keys = [b.GetName() for b in self.tree.GetListOfBranches()]

    def keys(self) -> list[str]:
        return list(self._keys)

    def arrays(self, names: Iterable[str]) -> dict[str, np.ndarray]:
        names = [name for name in dict.fromkeys(names) if name in self._keys]
        if not names:
            return {}
        data = self.ROOT.RDataFrame(TREE_NAME, str(self.path)).AsNumpy(names)
        return {name: np.asarray(data[name]) for name in names}


def load_arrays(path: str | Path, names: Iterable[str]) -> tuple[dict[str, np.ndarray], list[str], int]:
    tree = RootTree(path)
    available = set(tree.keys())
    names = list(dict.fromkeys(names))
    missing = [name for name in names if name not in available]
    return tree.arrays([name for name in names if name in available]), missing, tree.num_entries


def apply_selection(data: dict[str, np.ndarray], selections: list[str]) -> np.ndarray:
    n = len(next(iter(data.values()))) if data else 0
    mask = np.ones(n, dtype=bool)
    for selection in selections:
        if "=" not in selection:
            raise SystemExit(f"selection must be branch=value: {selection}")
        branch, value = selection.split("=", 1)
        if branch not in data:
            raise SystemExit(f"selection branch missing: {branch}")
        target = float(value)
        mask &= np.isclose(data[branch].astype(float), target, atol=1e-9)
    return mask


def finite(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return values[np.isfinite(values)]


def p1(x):
    return x


def p2(x):
    return 0.5 * (3.0 * x * x - 1.0)


def p4(x):
    x2 = x * x
    return 0.125 * (35.0 * x2 * x2 - 30.0 * x2 + 3.0)


def weight_a1a2(x, a1, a2):
    return 1.0 + a1 * p1(x) + a2 * p2(x)


def weight_a2a4(x, a2, a4):
    return 1.0 + a2 * p2(x) + a4 * p4(x)


def expected_probabilities(edges: np.ndarray, kind: str, coeffs: tuple[float, ...]) -> tuple[np.ndarray, float, float]:
    grid = np.linspace(-1.0, 1.0, 20001)
    if kind == "a1a2":
        raw = weight_a1a2(grid, coeffs[0], coeffs[1])
    elif kind == "a2a4":
        raw = weight_a2a4(grid, coeffs[0], coeffs[1])
    elif kind == "uniform":
        raw = np.ones_like(grid)
    else:
        raise ValueError(kind)
    clipped = np.clip(raw, 0.0, None)
    norm = np.trapz(clipped, grid)
    probs = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        local = np.linspace(lo, hi, 201)
        if kind == "a1a2":
            w = np.clip(weight_a1a2(local, coeffs[0], coeffs[1]), 0.0, None)
        elif kind == "a2a4":
            w = np.clip(weight_a2a4(local, coeffs[0], coeffs[1]), 0.0, None)
        else:
            w = np.ones_like(local)
        probs.append(float(np.trapz(w, local) / norm))
    probs = np.asarray(probs)
    return probs / probs.sum(), float(np.min(raw)), float(np.max(clipped))


def theory_cdf(x: np.ndarray, kind: str, coeffs: tuple[float, ...]) -> np.ndarray:
    x = np.asarray(x)
    grid = np.linspace(-1.0, 1.0, 40001)
    if kind == "a1a2":
        w = np.clip(weight_a1a2(grid, coeffs[0], coeffs[1]), 0.0, None)
    elif kind == "a2a4":
        w = np.clip(weight_a2a4(grid, coeffs[0], coeffs[1]), 0.0, None)
    elif kind == "uniform":
        w = np.ones_like(grid)
    else:
        raise ValueError(kind)
    cumulative = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(grid))])
    cumulative /= cumulative[-1]
    return np.interp(np.clip(x, -1.0, 1.0), grid, cumulative)


def ks_one_sample(values: np.ndarray, kind: str, coeffs: tuple[float, ...]) -> float:
    values = np.sort(finite(values))
    n = values.size
    if n == 0:
        return float("nan")
    cdf = theory_cdf(values, kind, coeffs)
    empirical_hi = np.arange(1, n + 1) / n
    empirical_lo = np.arange(0, n) / n
    return float(np.max(np.maximum(np.abs(empirical_hi - cdf), np.abs(cdf - empirical_lo))))


def ks_two_sample(a: np.ndarray, b: np.ndarray) -> float:
    a = np.sort(finite(a))
    b = np.sort(finite(b))
    if a.size == 0 or b.size == 0:
        return float("nan")
    values = np.sort(np.concatenate([a, b]))
    ca = np.searchsorted(a, values, side="right") / a.size
    cb = np.searchsorted(b, values, side="right") / b.size
    return float(np.max(np.abs(ca - cb)))


def chi2_ndf(values: np.ndarray, edges: np.ndarray, probs: np.ndarray, n_parameters: int = 0) -> tuple[float, int]:
    counts, _ = np.histogram(finite(values), bins=edges)
    expected = probs * counts.sum()
    mask = expected > 0.0
    ndf = max(1, int(np.sum(mask)) - 1 - n_parameters)
    chi2 = float(np.sum((counts[mask] - expected[mask]) ** 2 / expected[mask]))
    return chi2 / ndf, ndf


def moment_summary(values: np.ndarray) -> dict:
    values = finite(values)
    return {
        "n": int(values.size),
        "mean": float(np.mean(values)) if values.size else float("nan"),
        "mean_p1": float(np.mean(p1(values))) if values.size else float("nan"),
        "mean_p2": float(np.mean(p2(values))) if values.size else float("nan"),
        "mean_p4": float(np.mean(p4(values))) if values.size else float("nan"),
    }


def write_json(path: str | Path, data) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=2, sort_keys=True)
        f.write("\n")


def status_rank(status: str) -> int:
    return {"PASS": 0, "WARN": 1, "FAIL": 2}[status]


def combine_status(statuses: Iterable[str]) -> str:
    result = "PASS"
    for status in statuses:
        if status_rank(status) > status_rank(result):
            result = status
    return result


def safe_median(data: dict[str, np.ndarray], branch: str, mask: np.ndarray, default=float("nan")) -> float:
    if branch not in data:
        return default
    values = finite(data[branch][mask])
    return float(np.median(values)) if values.size else default

