#!/usr/bin/env python3
"""Prepare the run-specific directories used by the NS3 COMET pipeline."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--clean-dir",
        action="append",
        required=True,
        help="Directory to recreate empty; repeat for each directory.",
    )
    parser.add_argument(
        "--remove-file",
        action="append",
        default=[],
        help="File to remove if it exists; repeat for each file.",
    )
    return parser.parse_args()


def validate_directory(path: Path) -> None:
    resolved = path.resolve()
    if resolved == Path(resolved.anchor):
        raise RuntimeError(f"Refusing to recreate filesystem root: {resolved}")
    if resolved == Path.cwd().resolve():
        raise RuntimeError(f"Refusing to recreate the current directory: {resolved}")


def recreate_directory(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    validate_directory(resolved)
    if resolved.exists():
        if not resolved.is_dir():
            raise RuntimeError(f"Expected a directory to recreate, found a file: {resolved}")
        shutil.rmtree(resolved)
    resolved.mkdir(parents=True)
    return resolved


def remove_file(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if resolved.exists():
        if not resolved.is_file():
            raise RuntimeError(f"Expected a file to remove, found a directory: {resolved}")
        resolved.unlink()
    return resolved


def main() -> int:
    args = parse_args()
    clean_dirs = [recreate_directory(Path(value)) for value in args.clean_dir]
    removed_files = [remove_file(Path(value)) for value in args.remove_file]

    for path in clean_dirs:
        print(f"recreated_directory={display_path(path)}")
    for path in removed_files:
        print(f"removed_file={display_path(path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
