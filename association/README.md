# Association analysis: numbered analytical scripts

This folder contains only the analysis code shared on GitHub. Plotting code
has been removed.

| File | Purpose |
| --- | --- |
| `00_install_packages.R` | Installs the R packages used by scripts 01-05. |
| `01_calculate_burden_spatial_cdr.R` | Per-SA2 attributable mortality burden under `spatial_cdr`; writes raw Parquet files to `sa2_burden_by_demography/spatial_cdr/`. |
| `02_summarise_burden_by_state.R` | Summarises SA2 burden by Australian state. |
| `03_summarise_burden_by_period.R` | Summarises SA2 burden by multi-year period. |
| `04_gbd_factor_decomposition_baseline.R` | Baseline-to-future GBD/Shapley three-factor decomposition with total uncertainty. |
| `05_gdp_amr_inequality_analysis.R` | GDP-AMR SII/RII inequality analysis and manuscript tables (no plotting code). |

## Data assumptions

- Script 01 reads pre-computed exposure, population and mortality baselines.
  Those preprocessing scripts are not included here.
- Scripts 02-05 read the raw per-SA2 burden files produced by script 01.
- Script 04 reads a factor-input file normally produced by the adjacent-period
  GBD decomposition, which is not included in this repository. Only the
  baseline decomposition script is shared.
- `sa2_GDPpc_2000_2100_ssp123.feather` is required by script 05.

See the repository-level README for the expected `data/` layout and how to
set `AU_FIRE_ROOT`.
