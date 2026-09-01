#!/usr/bin/env python3
"""Run the shared subtype RAS consensus-difference report for NS5B."""

from __future__ import annotations

import runpy
import sys
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "hcv-ns3-one-ras-comet-build-workflow"
    / "scripts"
    / "build_ns3_subtype_ras_consensus_difference_summary.py"
)
sys.argv = [str(SCRIPT), "--gene", "NS5B", *sys.argv[1:]]
runpy.run_path(str(SCRIPT), run_name="__main__")
