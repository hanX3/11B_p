#!/usr/bin/env python3
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


def main():
  script = Path(__file__).resolve().parents[1] / "validation_angle_distribution" / "scripts" / "plot_compare_angular_cos.py"
  spec = spec_from_file_location("validation_angle_plot_compare", script)
  if spec is None or spec.loader is None:
    print(f"ERROR: unable to load {script}", file=sys.stderr)
    return 1
  module = module_from_spec(spec)
  spec.loader.exec_module(module)
  return module.main()


if __name__ == "__main__":
  raise SystemExit(main())
