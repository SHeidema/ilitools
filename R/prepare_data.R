# R/prepare_data.R
#
# prepare_data(df)
#
# Converts a raw surveillance data frame into the format expected by all
# model-fitting functions in this package.
#
# Input columns (required)
# ------------------------
#   location  <character>  Country / region of exposure
#   year      <integer>    Calendar year
#   week      <integer>    ISO week number (1-52)
#   cases     <integer>    Aggregated weekly case count
#
# Added columns
# -------------
#   i    <factor>   Integer index for each unique location
#   t    <integer>  Global week index: week + (year - min_year) * 52
#   cos  <numeric>  cos(2*pi*week/52)  — harmonic regressor
#   sin  <numeric>  sin(2*pi*week/52)  — harmonic regressor
#   pos  <factor>   glmmTMB::numFactor(cos, sin) — needed for cyclic GP term
#   Y    <integer>  Alias for cases (keeps model formulas uniform)

prepare_data <- function(df) {

  required <- c("location", "year", "week", "cases")
  missing  <- setdiff(required, names(df))
  if (length(missing) > 0)
    stop("prepare_data: missing column(s): ", paste(missing, collapse = ", "))

  # Location index (alphabetical, stable across calls)
  locs   <- sort(unique(df$location))
  loc_df <- data.frame(location = locs, i = factor(seq_along(locs)),
                       stringsAsFactors = FALSE)

  df <- merge(df, loc_df, by = "location")

  # Global time index
  df$t <- df$week + (df$year - min(df$year)) * 52

  # Harmonic regressors
  df$cos <- cos(2 * pi * df$week / 52)
  df$sin <- sin(2 * pi * df$week / 52)

  # Cyclic GP position (unit-circle coordinates encoded as a numFactor)
  df$pos <- glmmTMB::numFactor(df$cos, df$sin)

  # Response alias
  df$Y <- df$cases

  df[order(df$i, df$year, df$week), ]
}
