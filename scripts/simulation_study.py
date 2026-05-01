"""Placeholder script for the simulation study experiments and figures.

Run this script to reproduce all simulation study results from the paper:

    python scripts/simulation_study.py

Figures are saved to the ``figures/`` directory.
"""

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="GeoMix simulation study.")
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

    print("Simulation study — not yet implemented.")
    print(f"Figures will be saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
