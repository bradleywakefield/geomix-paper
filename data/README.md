# Data

This directory holds the raw and processed data files required to reproduce the paper results. Raw data files are not included in the repository due to confidentiality/size constraints.

## Directory structure

```
data/
├── application/       # IJV wind farm raw input data (place before running)
├── simulation/        # Simulation study raw input data (place before running)
├── processed/         # All files generated automatically by the scripts
└── README.md
```

---

## Required files (place before running)

### IJV Wind Farm application — `data/application/`

| File | Description |
|------|-------------|
| `cpt_profiles.rds` | CPT measurements at 10 cm depth resolution. Columns: `CPT_name`, `Easting_m`, `Northing_m`, `Bathymetry_UHR_5m_LAT`, `Bathymetry_MBES_0p5m_LAT`, `Depth_m_bsf`, `qn_MPa` |
| `cpt_stratigraphy.rds` | CPT sample data with synthetic ground model stratigraphy annotations. Columns include `name`, `depth`, `soilUnitID`, `soilUnit`, `eastingS`, `northingS` |
| `seismic_cdp.rds` | Seismic geophysics data (Common Depth Points). Columns include `id`, `easting`, `northing`, `depth`, `soilUnitID`, `soilUnit` |
| `misspec_cpt.csv` | Single CPT profile with Y1 (true stratigraphy) and Z1 (ground model stratigraphy) columns, used for the misspecification illustration figure |

### Simulation study — `data/simulation/`

| File | Description |
|------|-------------|
| `synthetic_data.rds` | Synthetic dataset generated for the simulation study. Contains `$data`, `$K`, `$dims`, `$beta`, `$penalty`, `$alpha`, `$sigma2`, `$tau`, `$lD` |

---

## Generated files (created automatically by scripts)

All generated files go into `data/processed/`. Running `scripts/01_process_data.R` creates:

| File | Description |
|------|-------------|
| `processed/locs.rds` | Combined spatial locations (CPT + CDP points) |
| `processed/cell_grid.rds` | 35 km spatial grid used for model lattice (local CRS, no projection) |
| `processed/data3D.RData` | Processed 3D dataset for the IJV application: `data`, `full_df`, `test_locs`, `cpt_locs`, `dims`, `K` |
| `processed/points_grid.rds` | Hexagonal grouping lattice (created by `scripts/utils/hex_grouping.R` after model setup) |

Running `scripts/05_map_figures.R` (called automatically from `scripts/02_fit_models.R`) additionally creates:

| File | Description |
|------|-------------|
| `processed/line_df.rds` | Cross-section line definitions (grid cell indices for each seismic line) — required by `scripts/04_application_results.R` |
| `processed/cell_grid_map.rds` | WGS84 version of the spatial grid for map visualisation |
