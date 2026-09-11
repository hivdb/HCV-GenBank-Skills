#!/usr/bin/env python3
"""Collect publication workbooks from active HCV Comet workflow outputs."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]

FILES = (
    (
        "outputs/Ref_same/12_pubmed_metadata_update/SubmissionSets.xlsx",
        "Table S1 - References.xlsx",
    ),
    (
        "outputs/comet-NS3-one-ras/17_merge-subtype-complete-profiles/NS3_Subtype_CompleteProfiles_Pos36_175.xlsx",
        "Table S2 - NS3.xlsx",
    ),
    (
        "outputs/comet-NS5A-one-ras/17_merge-subtype-complete-profiles/NS5A_Subtype_CompleteProfiles_Pos24_93.xlsx",
        "Table S3 - NS5A.xlsx",
    ),
    (
        "outputs/comet-NS5B-position-282-four-ras/17_merge-subtype-complete-profiles/NS5B_Subtype_CompleteProfiles_Pos150_321.xlsx",
        "Table S4 - NS5B.xlsx",
    ),
    (
        "outputs/comet-NS3-one-ras/23_build-subtype-ras-profile/NS3_Subtype_RAS_Profiles.xlsx",
        "Figure S1 - NS3 Profile - all subtypes.xlsx",
    ),
    (
        "outputs/comet-NS5A-one-ras/23_build-subtype-ras-profile/NS5A_Subtype_RAS_Profiles.xlsx",
        "Figure S3 - NS5A Profile - all subtypes.xlsx",
    ),
    (
        "outputs/comet-NS5B-position-282-four-ras/23_build-subtype-ras-profile/NS5B_Subtype_RAS_Profiles.xlsx",
        "Figure S5 - NS5B Profile - all subtypes.xlsx",
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--output-dir", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = args.repo_root.expanduser().resolve()
    output_dir = (args.output_dir or repo_root / "for_paper").expanduser()
    if not output_dir.is_absolute():
        output_dir = repo_root / output_dir
    output_dir = output_dir.resolve()

    missing = [str(repo_root / source) for source, _ in FILES if not (repo_root / source).is_file()]
    if missing:
        raise SystemExit("Missing source files:\n" + "\n".join(missing))

    copied = []
    output_dir.mkdir(parents=True, exist_ok=True)
    for source, destination_name in FILES:
        destination = output_dir / destination_name
        shutil.copy2(repo_root / source, destination)
        copied.append(str(destination))

    print(json.dumps({"output_dir": str(output_dir), "copied_files": copied}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
