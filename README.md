# GeoMix: A Bayesian Framework for Stratigraphic and Geotechnical Inference

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the code, data, and paper source files for reproducing all results and figures in the paper:

> **GeoMix: A Bayesian Framework for Stratigraphic and Geotechnical Inference**
> Bradley Wakefield et al.

---

## Abstract

Accurate characterisation of subsea sediments is essential for engineering projects such as offshore wind farms, but the challenging environment and large spatial domains involved make dense geotechnical measurements impractical. Complementing geotechnical measurements are seismic traces, which are relatively cheap and abundant. However, geotechnical interpretation of seismic traces relies on subjective judgement, which can lead to biased predictions of the underlying sediment structure with limited scope for validation.

To address this challenge, we introduce **GeoMix**, a fully Bayesian hierarchical framework for jointly inferring geological stratigraphy and geotechnical properties from seismic and geotechnical data. GeoMix uses a Potts Markov random field to represent spatially coherent geology organised into distinct stratigraphic layers. Labels derived from seismic interpretation are linked to the latent stratigraphy through an observation model that captures elevated uncertainty near transitions between strata. Conditional on the latent stratigraphy, geotechnical properties are modelled using class-specific Gaussian processes, with scalable likelihood evaluation enabled by a grouped Vecchia approximation.

This construction yields a unified latent-variable model capable of recovering stratigraphy, quantifying uncertainty, and predicting geotechnical properties across large three-dimensional domains. Simulation studies and an application to the IJmuiden Ver offshore wind farm demonstrate that GeoMix produces spatially coherent stratigraphic estimates and well-calibrated predictions of sediment characteristics that are substantially better than existing approaches.

---

## Repository Structure

```
geomix-paper/
├── data/                  # Datasets (see Data section below)
│   ├── raw/               # Raw input data
│   └── processed/         # Pre-processed data
├── figures/               # Generated figures (output)
├── paper/                 # LaTeX source files for the paper
├── R/                     # Core GeoMix source code
│   ├── model.R            # Bayesian hierarchical model
│   ├── potts.R            # Potts Markov random field
│   ├── gp.R               # Class-specific Gaussian processes
│   ├── vecchia.R          # Grouped Vecchia approximation
│   ├── inference.R        # MCMC / variational inference routines
│   └── utils.R            # Utility functions (plotting, I/O, etc.)
├── scripts/               # End-to-end experiment and figure scripts
│   ├── simulation_study.R
│   ├── ijmuiden_ver.R
│   └── make_figures.R
├── tests/                 # Unit and integration tests (testthat)
│   └── testthat/
├── .gitignore
├── .Rprofile              # Activates renv on project open
├── CITATION.cff
├── CONTRIBUTING.md
├── DESCRIPTION            # R package metadata
├── LICENSE
├── README.md
└── renv.lock              # Reproducible dependency snapshot (renv)
```

---

## Installation

### Prerequisites

- R ≥ 4.3
- [renv](https://rstudio.github.io/renv/) for dependency management

### Setup

```r
# Clone the repository, then in R:
install.packages("renv")    # if not already installed
renv::restore()             # install all dependencies from renv.lock
```

> **Note:** `renv.lock` is a placeholder until the full dependency list is
> finalised. Run `renv::init()` followed by `renv::snapshot()` to populate it
> once the package dependencies are known.

---

## Data

The IJmuiden Ver dataset used in the paper is subject to third-party licensing restrictions and cannot be redistributed directly. Please see `data/README.md` for instructions on how to obtain and prepare the data.

Synthetic data for the simulation studies is generated programmatically by the scripts in `scripts/`.

---

## Reproducing Results

### Simulation Studies

```bash
Rscript scripts/simulation_study.R
```

### IJmuiden Ver Application

```bash
Rscript scripts/ijmuiden_ver.R
```

### Generating All Figures

```bash
Rscript scripts/make_figures.R
```

Generated figures will be written to `figures/`.

---

## Running Tests

```r
# In R:
devtools::test()
```

or from the command line:

```bash
Rscript -e "devtools::test()"
```

---

## Citation

If you use GeoMix or this code in your work, please cite:

```bibtex
@article{wakefield2025geomix,
  title   = {{GeoMix}: A {Bayesian} Framework for Stratigraphic and Geotechnical Inference},
  author  = {Wakefield, Bradley and others},
  journal = {arXiv preprint},
  year    = {2025},
  url     = {https://github.com/bradleywakefield/geomix-paper}
}
```

A `CITATION.cff` file is also provided for software citation.

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Contributing

Contributions, bug reports, and feature requests are welcome. Please read [CONTRIBUTING.md](CONTRIBUTING.md) before opening an issue or pull request.

---

## Contact

For questions about the paper or code, please open a GitHub issue or contact the corresponding author via the details in the paper.
