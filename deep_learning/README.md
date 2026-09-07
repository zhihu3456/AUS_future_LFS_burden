# Deep-learning prediction and calibration

1. `01_training_firepm25_rawflux_2000_2020_unet_random_10fold_external_2021_2023.ipynb`
   — trains the fire-related PM2.5 U-Net with a locked random 10-fold split and
   external 2021-2023 validation.
2. `02_predict_firepm25_ssp_17models_4ssp_to_daily_tif.ipynb`
   — applies the trained model to 17 GCMs × 4 SSPs and writes daily TIF files.
3. `03_qdm_grid_calibration_from_yearly_tif.ipynb`
   — fits quantile-delta-mapping (QDM) parameters against historical data and
   calibrates the predicted future daily TIF series.

Run notebooks from this folder. Data paths default to `../data` relative to the
repository root and can be overridden with the `AU_FIRE_ROOT` environment
variable. Execution outputs have been cleared for sharing.
