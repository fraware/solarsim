"""Backward-compatible launcher. Prefer: ``python -m solarsim`` or the ``solarsim`` console script."""

from __future__ import annotations

import runpy
import sys
from pathlib import Path

if __name__ == "__main__":
    src = Path(__file__).resolve().parent / "src"
    sys.path.insert(0, str(src))
    runpy.run_module("solarsim", run_name="__main__")
