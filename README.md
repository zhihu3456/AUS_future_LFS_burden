# Fire-related air pollution health impact analysis (Australia)

Research code for quantifying the health burden attributable to fire-related
PM2.5 and O3 in Australia under historical and future SSP/GCM scenarios.

The repository contains:

- `association/` — the health-burden and inequality analysis code that is
  currently shared. Scripts are numbered and contain only analytical steps
  (no figure/plotting code).
- `deep_learning/` — model training, SSP inference, and quantile-delta-mapping
  (QDM) calibration notebooks.

## Repository layout

```text
.
├── association/
│   ├── 00_install_packages.R
│   ├── 01_calculate_burden_spatial_cdr.R
│   │     Per-SA2 burden calculation under the spatial_cdr assumption;
│   │     writes raw SA2 burden Parquet files to
│   │     `burden/sa2_burden_by_demography/spatial_cdr/`.
│   ├── 02_summarise_burden_by_state.R
│   │     Aggregates SA2 burden to Australian states (and national totals).
│   ├── 03_summarise_burden_by_period.R
│   │     Aggregates SA2 burden into multi-year periods
│   │     (2000-2020, 2021-2040, ... 2081-2100).
│   ├── 04_gbd_factor_decomposition_baseline.R
│   │     Baseline-to-future three-factor GBD decomposition with uncertainty.
│   └── 05_gdp_amr_inequality_analysis.R
│         GDP-AMR SII/RII inequality analysis; writes analysis tables
│         (no figure/plot code).
└── deep_learning/
    ├── 01_training_firepm25_rawflux_2000_2020_unet_random_10fold_external_2021_2023.ipynb
    ├── 02_predict_firepm25_ssp_17models_4ssp_to_daily_tif.ipynb
    └── 03_qdm_grid_calibration_from_yearly_tif.ipynb
```

## Association workflow

1. `01_calculate_burden_spatial_cdr.R`
   calculates per-SA2 annual attributable mortality burden (PM2.5 and O3,
   long- and short-term) for each coefficient simulation and GCM under the
   `spatial_cdr` demographic assumption.
2. `02_summarise_burden_by_state.R`
   aggregates the SA2 burden files into Australian-state summaries.
3. `03_summarise_burden_by_period.R`
   aggregates the SA2 burden files into multi-year period summaries.
4. `04_gbd_factor_decomposition_baseline.R`
   decomposes burden changes from the baseline period into population,
   exposure and baseline-mortality components, with Monte Carlo/GCM
   uncertainty. It reads a factor-input file (see Data below).
5. `05_gdp_amr_inequality_analysis.R`
   computes the GDP-AMR SII/RII inequality metrics and writes the annual
   results and manuscript tables. Plotting code is intentionally not shared.

The adjacent-period GBD decomposition script is not included; only the
baseline decomposition is kept.

## Data

Raw and intermediate data are **not included** in this repository because they
are large and may have separate licensing terms. The scripts expect the
following sibling data layout (customise with the `AU_FIRE_ROOT` environment
variable):

```text
data/
├── population/
│   ├── exposure_worldpop/
│   │   ├── hist_yearly/                 # historical exposure parquet files
│   │   └── ssp_yearly/                  # future exposure parquet files
│   ├── sa2_pop_cdr_cali_2000_2100.feather
│   └── sa2_hist_MorDoy_ratio.feather
└── burden/
    ├── betasim_short_long.feather       # coefficient simulations
    ├── sa2_used_statistics.rds          # SA2 grouping/lookup data
    ├── sa2_GDPpc_2000_2100_ssp123.feather
    ├── sa2_burden_by_demography/
    │   └── spatial_cdr/                 # produced by script 01
    ├── gbd_3factor_decomposition_total_uncertainty_adjacent/
    │   └── gbd_factor_inputs_by_sim_model_period.parquet
    └── gbd_3factor_decomposition_total_uncertainty_baseline/  # outputs
```

`05_gdp_amr_inequality_analysis.R` writes its tables into
`burden/GDP_AMR_inequality_simulation/`.

### Pointing R scripts at your data

Each R script searches upward for this repository's `.here` marker and, by
default, uses a sibling `data/` folder. To use another data location, export
the variable before running R:

```r
Sys.setenv(AU_FIRE_ROOT = "/your/project/data/root")
```

### Pointing the Python notebooks at your data

Run the notebooks from the `deep_learning/` folder; the default data root is
then `data/` at the repository root. Alternatively set the environment
variable `AU_FIRE_ROOT` to an absolute path before starting Jupyter:

```bash
export AU_FIRE_ROOT=/your/project/data/root
```

## Requirements

R scripts require:

```r
install.packages(c(
  "tidyverse", "arrow", "data.table", "lubridate", "tsModel",
  "future", "future.apply", "parallelly"
))
```

Python notebooks require the packages in `requirements.txt` and a PyTorch
installation with CUDA if you want to reproduce GPU inference/training:

```bash
pip install -r requirements.txt
```

## Notes

- Notebook execution outputs have been cleared before sharing.
- The deep-learning notebooks originally run on a Linux/HPC environment; the
  QDM notebook contains an explicit Linux platform check.
- The shared `association/` folder intentionally contains only the analysis
  scripts described above; figure/plot code and preprocessing code are not
  part of this repository.

## License

MIT. See `LICENSE`.
