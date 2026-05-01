"""Convenience script to regenerate all paper figures.

    python scripts/make_figures.py
"""

import subprocess
import sys
from pathlib import Path

SCRIPTS = [
    "scripts/simulation_study.py",
    "scripts/ijmuiden_ver.py",
]


def main() -> None:
    for script in SCRIPTS:
        print(f"Running {script} ...")
        result = subprocess.run(
            [sys.executable, script], check=False, capture_output=False
        )
        if result.returncode != 0:
            print(f"  WARNING: {script} exited with code {result.returncode}.")
    print("Done. Figures written to figures/.")


if __name__ == "__main__":
    main()
