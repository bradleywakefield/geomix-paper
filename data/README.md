# Data

This directory contains datasets used by GeoMix.

```
data/
├── raw/         # Original, unmodified input data
└── processed/   # Pre-processed data ready for model ingestion
```

## IJmuiden Ver Dataset

The IJmuiden Ver offshore wind farm dataset used in the paper is subject to
third-party licensing restrictions and cannot be redistributed in this
repository. To reproduce the real-data results, please follow the steps below:

1. Contact the data provider (details in the paper) to obtain access.
2. Place the raw data files in `data/raw/`.
3. Run the pre-processing script:
   ```bash
   # TODO: add the command/script you use to preprocess the IJmuiden data
   ```
   This will write the prepared files to `data/processed/`.

## Simulation Data

Synthetic datasets used in the simulation studies are generated programmatically
and are **not** stored in this directory. Run:

```bash
# TODO: add the command/script you use to generate simulation data
```

to reproduce the simulation results from scratch.
