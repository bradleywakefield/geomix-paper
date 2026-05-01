# Data

This directory holds the raw and processed data files required to reproduce the paper results. The raw data files are not included in the repository due to confidentiality/size constraints.

## Required files

Place the following files here before running the analysis scripts.

### IJV Wind Farm (application study)

| File | Description |
|------|-------------|
| `data/ijv_10cm.rds` | CPT measurements at 10 cm depth resolution. Columns: `CPT_name`, `Easting_m`, `Northing_m`, `Bathymetry_UHR_5m_LAT`, `Bathymetry_MBES_0p5m_LAT`, `Depth_m_bsf`, `qn_MPa` |
| `data/sample_wsynth.rds` | CPT sample data with synthetic ground model stratigraphy annotations. Columns include `name`, `depth`, `soilUnitID`, `soilUnit`, `eastingS`, `northingS` |
| `data/geophys/synth_df.rds` | Seismic geophysics data (Common Depth Points). Columns include `id`, `easting`, `northing`, `depth`, `soilUnitID`, `soilUnit` |
| `data/misspec.csv` | Single CPT profile with Y1 (true stratigraphy) and Z1 (ground model stratigraphy) columns, used for the misspecification illustration figure |

### Simulation study

| File | Description |
|------|-------------|
| `data/synthetic_list.rds` | Synthetic dataset generated for the simulation study. Contains `$data`, `$K`, `$dims`, `$beta`, `$penalty`, `$alpha`, `$sigma2`, `$tau`, `$lD` |

## Generated files (created automatically)

Running `scripts/01_process_data.R` creates:

| File | Description |
|------|-------------|
| `data/locs.rds` | Combined spatial locations (CPT + CDP points) |
| `data/cell_grid.rds` | 35 km spatial grid used for model lattice (local CRS, no projection) |
| `data/processed/data3D.RData` | Processed 3D dataset for the IJV application: `data`, `full_df`, `test_locs`, `cpt_locs`, `dims`, `K` |
| `data/processed/points_grid.rds` | Hexagonal grouping lattice (created by `scripts/utils/group_lattice.R` after model setup) |

Running `scripts/exploratory/eda.R` additionally creates:

| File | Description |
|------|-------------|
| `data/line_df.rds` | Cross-section line definitions (grid cell indices for each seismic line) — required by `scripts/04_application_results.R` |
| `data/cell_grid_map.rds` | WGS84 version of the spatial grid for map visualisation |
