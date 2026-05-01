"""Placeholder script for the IJmuiden Ver real-data application.

Run this script to reproduce the IJmuiden Ver results from the paper:

    python scripts/ijmuiden_ver.py

Requires the IJmuiden Ver dataset to be available in ``data/processed/``.
See ``data/README.md`` for instructions on obtaining the data.
"""

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="GeoMix — IJmuiden Ver application."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("data/processed"),
        help="Directory containing pre-processed data (default: data/processed/).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("figures"),
        help="Directory for saving figures (default: figures/).",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("IJmuiden Ver application — not yet implemented.")
    print(f"Data directory : {args.data_dir}")
    print(f"Output directory: {args.output_dir}")


if __name__ == "__main__":
    main()
