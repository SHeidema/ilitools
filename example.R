# example.R
#
# Minimal example: fit all seven models on the synthetic sample data and run
# the Shewhart monitoring chart.  Use this as a template for your own data.
#
# NOTE: ili_data_synth is purely synthetic (NegBin draws with sinusoidal
# seasonality).  It has no connection to real GeoSentinel records and is
# provided only to verify that the functions run correctly.

source("R/prepare_data.R")
source("R/fit_models.R")
source("R/shewhart.R")

# ---- 1. Load and prepare synthetic data ----
load("data/ili_data_synth.rda")   # loads ili_data_synth

df <- prepare_data(ili_data_synth)

# ---- 2. Fit models ----
# Fit one model (preferred):
fit <- fit_hybrid_autoregressive(df)

# Or fit all seven at once (clustering models require J to be specified):
models <- fit_all_models(df, J_direct = 3, J_spherical = 3)

# ---- 3. Shewhart monitoring ----
# Here we monitor on the same data; in practice df_monitor would be a
# new period prepared with the same prepare_data() call.
#
# alpha is the per-location significance level (Bonferroni-corrected):
#   alpha = (1 / ARL0_global) / n_locations
#   e.g. one false alarm per 3 years globally across 10 locations:
n_locations <- length(unique(df$location))
alpha       <- (1 / 156) / n_locations

result <- shewhart(fit, df, c = c(1, 2, 3), alpha = alpha)
