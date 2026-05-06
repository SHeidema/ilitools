# R/shewhart.R
#
# shewhart(fit, df, c, alpha)
#
# Applies a Shewhart control chart to a (potentially new) monitoring data
# frame using a previously fitted baseline model.
#
# The null hypothesis is:
#
#   H0 : mu_it  <=  c * mu_hat_it
#
# A signal is raised when the observed count exceeds the (1 - alpha) quantile
# of NegBin(c * mu_hat_it, theta).
#
# Arguments
# ---------
#   fit    Fitted model object (output of any fit_*() function).
#   df     Prepared data frame for the monitoring period (output of
#          prepare_data()).  Must contain the same locations as the training
#          data, or new locations will be handled via allow.new.levels = TRUE.
#   c      Multiplicative travel-volume allowance.  The paper uses c = 1, 2, 3.
#          Pass a single value or a numeric vector to get results for multiple
#          values at once.
#   alpha  Per-location significance level.  Typically set via a Bonferroni
#          correction: alpha = (1 / ARL0_global) / n_locations.
#          Default is 1/156/n_locations (≈ 1 false alarm per 3 years globally).
#
# Returns
# -------
#   The input data frame df enriched with, for each value of c:
#     c{c}_mu         numeric   c * mu_hat  (scaled predicted mean)
#     c{c}_upper_lim  numeric   upper control limit
#     c{c}_signal     logical   TRUE when cases > upper_lim

shewhart <- function(fit, df, c = 1, alpha = NULL) {

  # ---- predicted means ----
  # For clustering models, attach cluster assignments before predicting
  nd <- df
  cl <- attr(fit, "cluster_assignments")
  if (!is.null(cl)) {
    nd <- merge(nd, cl, by = "i", all.x = TRUE)
    if ("cluster" %in% names(nd))
      nd$cluster[is.na(nd$cluster)] <- "non_seasonal"
  }

  df$mu_hat <- predict(fit, newdata = nd, type = "response",
                       allow.new.levels = TRUE)

  # ---- dispersion parameter ----
  theta <- if (inherits(fit, "glmmTMB")) sigma(fit) else fit$theta

  # ---- default alpha: 1 false alarm per 3 years globally ----
  if (is.null(alpha)) {
    n_locs <- length(unique(df$location))
    alpha  <- (1 / 156) / n_locs
  }

  # ---- add columns for each requested c ----
  for (cv in c) {
    col_mu  <- paste0("c", cv, "_mu")
    col_ul  <- paste0("c", cv, "_upper_lim")
    col_sig <- paste0("c", cv, "_signal")

    df[[col_mu]]  <- cv * df$mu_hat
    df[[col_ul]]  <- qnbinom(p = 1 - alpha, mu = df[[col_mu]], size = theta)
    df[[col_sig]] <- df$cases > df[[col_ul]]
  }

  df
}
