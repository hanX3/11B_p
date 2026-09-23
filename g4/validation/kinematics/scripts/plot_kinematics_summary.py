#!/usr/bin/env python3
"""Small helper that runs the plotting-capable kinematics validators.

Most plots are produced by the channel-specific validation scripts.  This
wrapper exists so the validation directory has the task-doc entry point while
keeping plotting logic close to the checks that own the residuals.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True)
    parser.add_argument("--outdir", default="validation/kinematics/plots")
    parser.add_argument("--channel", default="all")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    commands = [
        [sys.executable, str(script_dir / "validate_four_momentum_closure.py"), "--root", args.root, "--channel", args.channel, "--outdir", args.outdir],
    ]
    if args.channel in ("all", "2"):
        commands.append([sys.executable, str(script_dir / "validate_direct_3alpha_phase_space.py"), "--root", args.root, "--channel", "2", "--outdir", args.outdir])

    status = 0
    for command in commands:
        rc = subprocess.call(command)
        if rc != 0:
            status = rc
    return status


if __name__ == "__main__":
    raise SystemExit(main())
