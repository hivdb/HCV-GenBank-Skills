#!/usr/bin/env python3
"""Copy final workflow artifacts into its output report directory."""

from __future__ import annotations
import argparse
import shutil
from pathlib import Path

README_GENE = {"NS3": "NS3", "NS5A": "NS5A_NTD", "NS5B": "NS5B"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workflow-output-dir", type=Path, required=True)
    parser.add_argument("--gene", choices=tuple(README_GENE), required=True)
    args = parser.parse_args()
    root = args.workflow_output_dir

    def artifact(steps: tuple[str, ...], filename: str) -> Path:
        for step in steps:
            matches = sorted(root.glob(f"*_{step}/{filename}"))
            if len(matches) == 1:
                return matches[0]
            if len(matches) > 1:
                raise SystemExit(
                    f"Expected one {step}/{filename} under {root}, found {len(matches)}"
                )
        raise SystemExit(f"Missing {' or '.join(steps)}/{filename} under {root}")

    gene = args.gene
    sources = (
        artifact(
            ("compare-reference-consensus", "publish-ictv-report"),
            f"README_Subtype_Consensus_Mutations_{README_GENE[gene]}.docx",
        ),
        artifact(
            ("add-nonconsensus-row",), f"{gene}_Combined_RAS_Profiles_Annotated.xlsx"
        ),
        artifact(("build-combined-ras-reports",), f"{gene}_Combined_RAS_Profiles.xlsx"),
        artifact(("build-subtype-ras-profile",), f"{gene}_Subtype_RAS_Profiles.xlsx"),
        artifact(
            ("merge-subtype-complete-profiles",),
            f"{gene}_Subtype_CompleteProfiles_Merged.xlsx",
        ),
    )
    destination = root / "report"
    destination.mkdir(parents=True, exist_ok=True)
    for source in sources:
        shutil.copy2(source, destination / source.name)
        print(f"copied={source.name}")


if __name__ == "__main__":
    main()
