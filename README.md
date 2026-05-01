# geomix-paper

Reproducible analysis code for the paper:

> **GeoMix: A Bayesian Hierarchical Model for Joint Inversion of Geotechnical and Geophysical Data**  
> Bradley Wakefield

The `geomix` R package implementing the core model is at <https://github.com/bradleywakefield/geomix>.

---

## Repository structure

```
geomix-paper/
├── paper/                          # LaTeX source and compiled figures
├── scripts/
│   ├── run_analysis.R              # Main entry point — runs the full pipeline
│   ├── 01_process_data.R           # Data processing: CPT + geophysics → data/processed/
│   ├── 02_fit_models.R             # MCMC inference + predictions (simulation & IJV application)
│   ├── 03_simulation_results.R     # Figures for the simulation study
│   ├── 04_application_results.R    # Figures for the IJV wind farm application
│   ├── 05_map_figures.R            # Study area map figures + saves data/line_df.rds
│   ├── utils/
│   │   ├── prediction_scoring.R    # Scoring rules and evaluation metrics
│   │   ├── geowarp_builder.R       # GeoWarp competing model factory
│   │   ├── fit_competing_models.R  # Fits GeoWarp / GP / LM competing models
│   │   └── hex_grouping.R          # Hexagonal grouping of lattice points
│   └── exploratory/
│       ├── eda.R                   # Study area maps (standalone version of 05_map_figures.R)
│       ├── cross_section_plots.R   # Cross-section diagnostic plots
│       ├── misspec.R               # Misspecification illustration
│       └── generate_simulation_data.R  # Synthetic data generator (demo)
├── data/
│   └── README.md                   # Required input files and generated outputs
└── results/
    ├── figures/                    # PDF/PNG figures saved by the scripts
    ├── simulation/                 # MCMC samples and predictions (simulation study)
    └── application/                # MCMC samples and predictions (IJV application)
```

---

## Reproducing the results

### 1. Install the geomix package

```r
devtools::install_github("bradleywakefield/geomix")
```

### 2. Place data files

See [data/README.md](data/README.md) for the list of required input files and where to put them.

### 3. Run the analysis

Open an R session with the project root as the working directory and run:

```r
source("scripts/run_analysis.R")
```

This executes the full pipeline in order:
1. `01_process_data.R` — processes raw CPT and geophysics data
2. `02_fit_models.R` — runs MCMC chains and produces all figures; internally runs:
   - `03_simulation_results.R` — simulation study figures
   - `05_map_figures.R` — study area maps; saves `data/line_df.rds` (required by step below)
   - `04_application_results.R` — IJV wind farm application figures

> **Runtime note:** MCMC sampling is computationally intensive and was originally run on a 50-core server. `confirm_run()` prompts in `02_fit_models.R` let you skip sampling and load pre-saved results instead.

### 4. Optional: additional diagnostic figures

These scripts can be run independently after the main pipeline:

```r
source("scripts/exploratory/cross_section_plots.R")  # Cross-section diagnostic plots
source("scripts/exploratory/misspec.R")              # Misspecification figure
```

---

## Output figures

All figures are saved to `results/figures/`. Key outputs:

| File | Description |
|------|-------------|
| `syn-cross.pdf` | Simulation study cross-section |
| `syn-post.pdf` | Posterior sample histograms (simulation) |
| `syn-params.pdf` | Parameter recovery (simulation) |
| `syn_trace*.png` | MCMC trace plots (simulation) |
| `ijv-test.pdf` | Test CPT predictions (IJV application) |
| `ijv-z1.pdf`, `ijv-z2.pdf`, `ijv-y1.pdf` | Grid prediction maps |
| `ijv-cross*.pdf` | Cross-section plots along seismic lines |
| `ijv_trace*.png` | MCMC trace plots (IJV application) |
| `misspec.pdf` | Misspecification illustration |
| `map*.pdf` | Study area location maps |