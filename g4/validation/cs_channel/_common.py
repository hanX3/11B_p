"""Shared helpers for validation/cs_channel/*.py.

Loads branches from the "tr" TTree in a reaction_merged_*.root file with
uproot, reports MISSING BRANCH instead of crashing, and prints results in the
[PASS]/[WARN]/[FAIL] format required by gpt/codex_cs_channel_validation.md.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Iterable

import numpy as np

TREE_NAME = "tr"
DEFAULT_TOL = 1.0e-9


def load_uproot():
    try:
        import uproot  # type: ignore
    except ImportError as exc:
        raise SystemExit("uproot is not installed; cannot validate ROOT files.") from exc
    return uproot


def open_tree(path: Path):
    uproot = load_uproot()
    root_file = uproot.open(path)
    if TREE_NAME not in root_file:
        raise SystemExit(f"{path}: missing tree '{TREE_NAME}'")
    return root_file[TREE_NAME]


def load_branches(tree, names: Iterable[str]) -> tuple[dict[str, np.ndarray], list[str]]:
    """Return (data, missing) where data holds only the branches that exist."""
    available = set(tree.keys())
    missing = [name for name in names if name not in available]
    for name in missing:
        print(f"MISSING BRANCH: {name}")
    data = {name: tree[name].array(library="np") for name in names if name in available}
    return data, missing


def decode_char_branch(values: np.ndarray) -> np.ndarray:
    decoded = []
    for value in values:
        if isinstance(value, (bytes, bytearray)):
            decoded.append(value.decode("utf-8", errors="replace").rstrip("\x00"))
        else:
            decoded.append(str(value).rstrip("\x00"))
    return np.asarray(decoded, dtype=object)


def max_abs_diff(a: np.ndarray, b: np.ndarray) -> float:
    if a.size == 0:
        return float("nan")
    return float(np.nanmax(np.abs(a - b)))


def rel_diff(observed: float, expected: float) -> float:
    denom = abs(expected)
    if denom < 1e-300:
        return float("nan") if abs(observed) > 1e-300 else 0.0
    return abs(observed - expected) / denom


class Report:
    """Collects [PASS]/[WARN]/[FAIL] lines and tracks overall exit status."""

    def __init__(self, path: Path):
        self.path = path
        self.status = "PASS"
        self.lines: list[str] = []

    def _bump(self, level: str) -> None:
        order = {"PASS": 0, "WARN": 1, "FAIL": 2}
        if order[level] > order[self.status]:
            self.status = level

    def emit(self, level: str, label: str) -> None:
        self._bump(level)
        line = f"[{level}] {label}"
        print(line)
        self.lines.append(line)

    def fail_detail(
        self,
        branch: str,
        expected,
        observed,
        entry_index=None,
        related_setting: str | None = None,
    ) -> None:
        abs_diff = abs(observed - expected) if _is_number(observed) and _is_number(expected) else float("nan")
        relative = rel_diff(observed, expected) if _is_number(observed) and _is_number(expected) else float("nan")
        print(f"    branch name       : {branch}")
        print(f"    expected value    : {expected}")
        print(f"    observed value    : {observed}")
        print(f"    absolute difference: {abs_diff}")
        print(f"    relative difference: {relative}")
        if entry_index is not None:
            print(f"    entry index       : {entry_index}")
        if related_setting is not None:
            print(f"    related macro setting: {related_setting}")

    def check_close(
        self,
        label: str,
        observed: np.ndarray,
        expected: np.ndarray,
        branch: str,
        tol: float = DEFAULT_TOL,
        rtol: float = 1.0e-9,
        related_setting: str | None = None,
    ) -> bool:
        # Cross sections here span many orders of magnitude (internal G4 units),
        # so a pure absolute tolerance falsely flags float-roundoff on large
        # values; use atol + rtol*|expected|, matching numpy.isclose semantics.
        diff = max_abs_diff(observed, expected)
        scale = float(np.nanmax(np.abs(expected))) if expected.size else 0.0
        threshold = tol + rtol * scale
        ok = np.isfinite(diff) and diff < threshold
        if ok:
            self.emit("PASS", f"{label} (max abs diff = {diff:.3g}, threshold = {threshold:.3g})")
        else:
            self.emit("FAIL", f"{label} (max abs diff = {diff:.3g}, threshold = {threshold:.3g})")
            bad_idx = int(np.nanargmax(np.abs(observed - expected))) if observed.size else None
            self.fail_detail(
                branch,
                float(expected[bad_idx]) if bad_idx is not None else expected,
                float(observed[bad_idx]) if bad_idx is not None else observed,
                entry_index=bad_idx,
                related_setting=related_setting,
            )
        return ok

    def summary_line(self) -> str:
        return f"[{self.status}] {self.path}"


def _is_number(value) -> bool:
    try:
        float(value)
        return True
    except (TypeError, ValueError):
        return False


def first_value(data: dict[str, np.ndarray], name: str):
    arr = data.get(name)
    if arr is None or arr.size == 0:
        return None
    return arr[0]
