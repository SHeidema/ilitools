# ilitools

R functions for fitting baseline epidemiology models and running Shewhart outbreak detection on traveler surveillance data, as described in:

*Global Detection of Respiratory Illness Outbreaks in Travelers: A Statistical Approach using GeoSentinel Data.*

---

## Files

```
ilitools/
├── R/
│   ├── prepare_data.R   — convert raw data into the expected format
│   ├── fit_models.R     — seven baseline models + convenience wrapper
│   └── shewhart.R       — Shewhart control chart monitoring
└── data/
    └── ili_data_synth.rda  — synthetic dataset for testing (not real GeoSentinel data)
```

---

## Quick start

```r
source("R/prepare_data.R")
source("R/fit_models.R")
source("R/shewhart.R")

load("data/ili_data_synth.rda")   # loads ili_data_synth (synthetic — not real GeoSentinel data)

df <- prepare_data(ili_data_synth)

# Fit the preferred model (Eq. 1 in the paper)
fit <- fit_hybrid_autoregressive(df)

# Or fit all seven at once
models <- fit_all_models(df)

# Monitor a new period (here we reuse the same data as an example)
result <- shewhart(fit, df, c = c(1, 2, 3))
```

---

## Sample data

`ili_data_synth` is a **fully synthetic** dataset generated from Negative Binomial draws with sinusoidal seasonality (10 locations × 5 years × 52 weeks). It has no connection to real GeoSentinel records and exists solely so you can verify the functions run correctly before substituting your own data.

## Data format

`prepare_data()` expects a data frame with four columns:

| Column     | Type      | Description                         |
|------------|-----------|-------------------------------------|
| `location` | character | Country / region of exposure        |
| `year`     | integer   | Calendar year                       |
| `week`     | integer   | ISO week number (1–52)              |
| `cases`    | integer   | Aggregated weekly case count        |

---

## The seven models

| Function                       | Description                                              |
|-------------------------------|----------------------------------------------------------|
| `fit_full_fixed(df)`          | Country-specific fixed intercept + fixed harmonics       |
| `fit_full_random(df)`         | Shared mean, random intercept + random harmonics         |
| `fit_hybrid(df)`              | Fixed intercept, random harmonics (no autocorrelation)   |
| `fit_direct_clustering(df, J)`| Fixed model with k-means–clustered seasonal coefficients |
| `fit_spherical_clustering(df, J)` | Fixed model with spherical k-means–clustered phase  |
| `fit_autoregressive(df)`      | Fixed intercept + cyclic GP (no explicit harmonics)      |
| `fit_hybrid_autoregressive(df)` | **Preferred**: fixed intercept, random harmonics, cyclic GP |

---

## Shewhart monitoring

```r
result <- shewhart(fit, df_monitor, c = c(1, 2, 3), alpha = 1e-4)
```

`shewhart()` returns `df_monitor` enriched with, for each value of `c`:

| Column           | Description                             |
|------------------|-----------------------------------------|
| `mu_hat`         | Predicted mean under baseline           |
| `c{c}_mu`        | Scaled mean: `c × mu_hat`              |
| `c{c}_upper_lim` | Upper control limit                     |
| `c{c}_signal`    | `TRUE` when `cases > upper_lim`         |

`alpha` is the **per-location** significance level. Following the paper, set it via a Bonferroni correction:

```r
alpha <- (1 / ARL0_global) / n_locations
# e.g. for ARL0 = 156 weeks and 64 locations:
alpha <- (1 / 156) / 64   # ≈ 1e-4
```

If `alpha` is omitted, this default is applied automatically.

---

## Dependencies

```r
install.packages(c("glmmTMB", "MASS", "skmeans"))
```
