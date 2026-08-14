#!/usr/bin/env python3
"""Run NS5A-specific profile-alignment QC with fixed gene thresholds."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    repo_root = Path(__file__).resolve().parents[3]
    shared_validator = repo_root / "hcv-ns3-comet-build-workflow" / "scripts" / "validate_ns3_profile_alignment.py"
    command = [
        sys.executable,
        str(shared_validator),
        *sys.argv[1:],
        "--reference-gene",
        "NS5A_NTD",
        "--ras-positions",
        "24,26,28,29,30,31,32,38,58,62,92,93",
        "--high-divergence-percent",
        "30",
        "--min-divergence-coverage",
        "150",
    ]
    return subprocess.run(command, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
