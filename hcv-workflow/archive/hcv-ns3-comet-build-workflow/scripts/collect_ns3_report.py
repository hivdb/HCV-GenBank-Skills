#!/usr/bin/env python3
"""Copy the final NS3 report artifacts into one report directory."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--readme-docx", type=Path, required=True)
    parser.add_argument("--annotated-combined-profile", type=Path, required=True)
    parser.add_argument("--combined-profile", type=Path, required=True)
    parser.add_argument("--subtype-ras-profile", type=Path, required=True)
    parser.add_argument("--subtype-complete-profiles", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sources = (
        args.readme_docx,
        args.annotated_combined_profile,
        args.combined_profile,
        args.subtype_ras_profile,
        args.subtype_complete_profiles,
    )
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        raise SystemExit("Missing report artifact(s): " + "; ".join(missing))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for source in sources:
        shutil.copy2(source, args.output_dir / source.name)
        print(f"copied={source.name}")


if __name__ == "__main__":
    main()
