# geomix-paper

Reproducible analysis code for the paper:

> **GeoMix: A Bayesian Hierarchical Model for Joint Inversion of Geotechnical and Geophysical Data**  
> Bradley Wakefield

The `geomix` R package implementing the core model is at <https://github.com/bradleywakefield/geomix>.

---

## Repository structure

```
geomix-paper/
├── run_analysis.R              # Full pipeline entry point (simulation + application)
├── run_simulation.R            # Simulation-only entry point (no proprietary data required)
├── run_application.R           # IJV application entry point
├── renv.lock                   # Package dependency lockfile
├── scripts/
│   ├── application/
│   │   ├── 01_process_data.R   # CPT + geophysics → data/processed/
│   │   ├── 02_map_figures.R    # Study area map figures
│   │   ├── 03_fit_models.R     # MCMC inference + predictions (IJV application)
│   │   └── 04_results.R        # Figures and tables for the IJV application
│   ├── simulation/
│   │   ├── 01_generate_data.R  # Synthetic dataset generator
│   │   ├── 02_fit_models.R     # MCMC inference + predictions (simulation study)
│   │   └── 03_results.R        # Figures and tables for the simulation study
│   ├── utils/
│   │   ├── ask_prompts.R           # Interactive run-mode prompts
│   │   ├── prediction_scoring.R    # Scoring rules and evaluation metrics
│   │   ├── geowarp_builder.R       # GeoWarp competing model factory
│   │   ├── fit_competing_models.R  # Fits GeoWarp / GP / LM competing models
│   │   └── hex_grouping.R          # Hexagonal grouping of lattice points
│   └── exploratory/
│       ├── eda.R                   # Study area maps (standalone EDA figures)
│       ├── cross_section_plots.R   # Cross-section diagnostic plots
│       └── misspec.R               # Misspecification illustration
├── data/
│   ├── application/                # IJV wind farm raw input data (place before running)
│   ├── simulation/                 # Simulation study data
│   ├── processed/                  # All files generated automatically by the scripts
│   └── README.md                   # Required input files and generated outputs
└── results/
    ├── figures/                    # PDF/PNG figures saved by the scripts
    ├── tables/                     # CSV summary tables
    ├── simulation/                 # MCMC samples and predictions (simulation study)
    └── application/                # MCMC samples and predictions (IJV application)
```

---

## Dependencies

### R packages

All CRAN packages and their transitive dependencies are captured in [`renv.lock`](renv.lock). The two GitHub-only packages must be installed from source:

| Package | Source | Purpose |
|---------|--------|---------|
| `geomix` | `bradleywakefield/geomix` (GitHub) | Core GeoMix model |
| `geowarp` | `mbertolacci/geowarp` (GitHub) | GeoWarp competing model |
| `nimble` | CRAN | MCMC engine |
| `posterior` | CRAN | MCMC diagnostics |
| `scoringRules` | CRAN | Proper scoring rules |
| `mvtnorm` | CRAN | Multivariate normal utilities |
| `sf` | CRAN | Spatial data handling |
| `rnaturalearth` + `rnaturalearthdata` | CRAN | Basemap data |
| `concaveman` | CRAN | Concave hull computation |
| `tidyverse` | CRAN | Data wrangling and plotting |
| `ggplot2` / `ggtext` / `ggh4x` / `patchwork` / `cowplot` | CRAN | Figure layout |
| `matrixStats` | CRAN | Vectorised matrix statistics |
| `reshape2` | CRAN | Data reshaping |
| `slider` | CRAN | Sliding-window operations |
| `tictoc` | CRAN | Timing utilities |
| `yardstick` | CRAN | Model performance metrics |

---

## Reproducing the results

### 1. Restore the R environment

The easiest way to install all dependencies at the exact versions used for this paper is via [`renv`](https://rstudio.github.io/renv/):

```r
install.packages("renv")
renv::restore()
```

This reads [`renv.lock`](renv.lock) and installs every package (CRAN and GitHub) into a project-local library.

Alternatively, install the two GitHub packages manually and let `run_simulation.R` / `run_application.R` install the rest on first run:

```r
install.packages("devtools")
devtools::install_github("bradleywakefield/geomix")
devtools::install_github("mbertolacci/geowarp")
```

### 2. Place data files

See [data/README.md](data/README.md) for the list of required input files and where to put them.

### 3. Run the analysis

Open an R session with the project root as the working directory and run:

```r
source("run_analysis.R")
```

To run only the simulation study (no proprietary data required):

```r
source("run_simulation.R")
```

`run_analysis.R` executes the full pipeline in order:

**Simulation study** (`run_simulation.R`):
1. `scripts/simulation/02_fit_models.R` — MCMC inference on synthetic data; saves samples and predictions to `results/simulation/`
2. `scripts/simulation/03_results.R` — produces all simulation figures and tables

**IJV application** (`run_application.R`):
1. `scripts/application/01_process_data.R` — processes raw CPT and geophysics data; saves `data/processed/`
2. `scripts/application/02_map_figures.R` — study area map figures
3. `scripts/application/03_fit_models.R` — MCMC inference on IJV data; saves samples and predictions to `results/application/`
4. `scripts/application/04_results.R` — produces all application figures and tables

> **Runtime note:** MCMC sampling is computationally intensive and was originally run on a 50-core server. `ask_run_mode()` / `ask_pred_mode()` prompts in the fit-models scripts let you skip sampling and load pre-saved results instead.

### 4. Optional: additional diagnostic figures

These scripts can be run independently after the main pipeline:

```r
source("scripts/exploratory/cross_section_plots.R")  # Cross-section diagnostic plots
source("scripts/exploratory/misspec.R")              # Misspecification figure
source("scripts/exploratory/eda.R")                  # Study area EDA figures
```

---

## Output figures and tables

Figures are saved to `results/figures/`; summary tables to `results/tables/`.

### Simulation study

| File | Description |
|------|-------------|
| `syn-cross.pdf` | Cross-section of simulated stratigraphy and predictions |
| `syn-post.pdf` | Posterior sample histograms |
| `syn-params.pdf` | Parameter recovery plots |
| `syn_trace1–4.png` | MCMC parameter trace plots (4 chains) |
| `syn_traceY1.png` | MCMC stratigraphy trace plot |
| `results/tables/sim_stats.csv` | Y2 prediction metrics (all models) |
| `results/tables/sim_prob_stats.csv` | Y1 classification metrics |

### IJV application

| File | Description |
|------|-------------|
| `ijv-test.pdf` | Test CPT predictions |
| `ijv-test-small.pdf` | Compact version of test CPT predictions |
| `ijv-z1.pdf` | Grid prediction maps (Z1 stratigraphy) |
| `ijv-z2y1.pdf` | Grid prediction maps (Z2 geotechnical + Y1) |
| `ijv-cross.pdf` | Full cross-section plot |
| `ijv-crossA.pdf` | Cross-section (subset A) |
| `ijv-crossB.pdf` | Cross-section (subset B) |
| `ijv_trace1–4.png` | MCMC parameter trace plots (4 chains) |
| `ijv_traceY1.png` | MCMC stratigraphy trace plot |

### Study area maps and exploratory

| File | Description |
|------|-------------|
| `map1.pdf`, `map2.pdf`, `map3.pdf` | Study area location maps |
| `map2sub.pdf`, `map3sub.pdf` | Inset maps |
| `misspec.pdf` | Misspecification illustration |
| `misspec2.pdf` | Misspecification illustration (alternative view) |
