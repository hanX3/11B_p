#!/usr/bin/env python3
"""Thin plotting wrapper around validate_legendre_distributions.py."""

from __future__ import annotations

import subprocess
import sys


if __name__ == "__main__":
    raise SystemExit(subprocess.call([sys.executable, "validation/angular_decay/validate_legendre_distributions.py", *sys.argv[1:]]))

