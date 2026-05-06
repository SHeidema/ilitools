# R/fit_models.R
#
# One fitting function per model from the paper, plus a convenience wrapper.
#
# All functions take a prepared data frame (output of prepare_data()) and
# return a fitted model object on which predict(..., type = "response") works.
#
# Dependencies: glmmTMB, MASS, skmeans

library(glmmTMB)
library(MASS)


# ------------------------------------------------------------------------------
# 1. Full fixed-effects
#
#    log(mu_it) = beta_i + beta2_i * cos(2pi*t/52) + beta3_i * sin(2pi*t/52)
#
#    Every location gets its own intercept and seasonal shape.
# ------------------------------------------------------------------------------
fit_full_fixed <- function(df) {
  glm.nb(Y ~ 0 + i + i:cos + i:sin, data = df)
}


# ------------------------------------------------------------------------------
# 2. Full random-effects
#
#    log(mu_it) = b1_i + b2_i * cos + b3_i * sin
#    b1_i ~ N(beta1, sigma1^2),  b2_i ~ N(beta2, sigma2^2),  b3_i ~ N(beta3, sigma3^2)
# ------------------------------------------------------------------------------
fit_full_random <- function(df) {
  glmmTMB(Y ~ cos + sin + (1 + cos + sin | i),
          data = df, family = nbinom2())
}


# ------------------------------------------------------------------------------
# 3. Hybrid  (fixed intercept, random seasonality, no autocorrelation)
#
#    log(mu_it) = beta_i + b2_i * cos + b3_i * sin
# ------------------------------------------------------------------------------
fit_hybrid <- function(df) {
  glmmTMB(Y ~ 0 + i + (0 + cos + sin | i),
          data = df, family = nbinom2())
}


# ------------------------------------------------------------------------------
# 4. Direct clustering
#
#    Step 1: fit full fixed-effects model.
#    Step 2: k-means cluster locations in (beta2, beta3) space with J clusters.
#    Step 3: refit with shared seasonal coefficients within each cluster.
#
#    Arguments
#    ---------
#    df  : prepared data frame
#    J   : number of clusters (user-specified)
# ------------------------------------------------------------------------------
fit_direct_clustering <- function(df, J = 3) {

  # Step 1
  m_sat   <- glm.nb(Y ~ 0 + i + i:cos + i:sin, data = df)
  cf      <- coef(m_sat)
  b2      <- cf[grep(":cos$", names(cf))]
  b3      <- cf[grep(":sin$", names(cf))]

  # Step 2
  cl      <- kmeans(data.frame(b2, b3), centers = J, nstart = 20)$cluster
  cl_df   <- data.frame(i = factor(seq_along(cl)), cluster = factor(cl))

  # Step 3
  df_cl   <- merge(df, cl_df, by = "i")
  fit     <- glm.nb(Y ~ 0 + i + cluster:cos + cluster:sin, data = df_cl)

  attr(fit, "cluster_assignments") <- cl_df
  fit
}


# ------------------------------------------------------------------------------
# 5. Spherical clustering
#
#    Step 1: fit full fixed-effects model.
#    Step 2: test each location for non-zero seasonal amplitude (Hotelling T^2).
#            Locations that fail the test are placed in a "non-seasonal" group.
#    Step 3: spherical k-means on the phase angles of the remaining locations.
#    Step 4: refit with cluster-shared seasonal coefficients.
#
#    Arguments
#    ---------
#    df              : prepared data frame
#    J               : number of seasonal clusters (user-specified; the
#                      non-seasonal group is added on top)
#    alpha_hotelling : significance level for the Hotelling T^2 test
# ------------------------------------------------------------------------------
fit_spherical_clustering <- function(df, J = 3, alpha_hotelling = 0.05) {
  library(skmeans)

  n_t   <- length(unique(df$t))
  m_sat <- glm.nb(Y ~ 0 + i + i:cos + i:sin, data = df)
  cf    <- coef(m_sat)
  vc    <- vcov(m_sat)
  b2v   <- cf[grep(":cos$", names(cf))]
  b3v   <- cf[grep(":sin$", names(cf))]
  I     <- length(b2v)

  # Hotelling T^2 p-value for each location
  pvals <- vapply(seq_len(I), function(k) {
    nm2  <- paste0("i", k, ":cos");  nm3 <- paste0("i", k, ":sin")
    if (!(nm2 %in% rownames(vc))) return(1)
    Sig  <- matrix(c(vc[nm2, nm2], vc[nm2, nm3],
                     vc[nm3, nm2], vc[nm3, nm3]), 2, 2)
    bv   <- c(b2v[k], b3v[k])
    T2   <- n_t * t(bv) %*% solve(Sig) %*% bv
    F_   <- as.numeric(T2) * (n_t - 2) / (2 * (n_t - 1))
    1 - pf(F_, 2, n_t - 2)
  }, numeric(1))

  seas_idx <- which(pvals <= alpha_hotelling)
  cl_vec   <- rep("non_seasonal", I)

  if (length(seas_idx) >= J) {
    mat            <- matrix(c(b2v[seas_idx], b3v[seas_idx]), ncol = 2)
    sk             <- skmeans(mat, k = J, method = "pclust")
    cl_vec[seas_idx] <- paste0("s", sk$cluster)
  } else if (length(seas_idx) > 0) {
    cl_vec[seas_idx] <- "s1"
  }

  cl_df   <- data.frame(i = factor(seq_len(I)), cluster = cl_vec)
  df_cl   <- merge(df, cl_df, by = "i")
  fit     <- glm.nb(Y ~ 0 + i + cluster:cos + cluster:sin, data = df_cl)

  attr(fit, "cluster_assignments") <- cl_df
  fit
}


# ------------------------------------------------------------------------------
# 6. Autoregressive  (cyclic GP only, no explicit harmonic terms)
#
#    log(mu_it) = beta_i + eps_it
#    eps_it ~ GP with cyclic covariance (Eq. 2 in the paper)
# ------------------------------------------------------------------------------
fit_autoregressive <- function(df) {
  glmmTMB(Y ~ 0 + i + exp(pos + 0 | i),
          data = df, family = nbinom2(), REML = TRUE)
}


# ------------------------------------------------------------------------------
# 7. Hybrid autoregressive  (PREFERRED MODEL — Eq. 1 in the paper)
#
#    log(mu_it) = beta_i + b2_i * cos + b3_i * sin + eps_it
#
#    Fixed intercept per location, random seasonal harmonics, cyclic GP.
# ------------------------------------------------------------------------------
fit_hybrid_autoregressive <- function(df) {
  glmmTMB(Y ~ 0 + i + (0 + cos + sin | i) + exp(pos + 0 | i),
          data = df, family = nbinom2(), REML = TRUE)
}


# ------------------------------------------------------------------------------
# Convenience wrapper — fit all seven models at once
#
# Returns a named list of fitted model objects.  Failed models are stored as
# NULL (with a warning) so the others are unaffected.
# ------------------------------------------------------------------------------
fit_all_models <- function(df, J_direct = 3, J_spherical = 3,
                           alpha_hotelling = 0.05) {
  fitters <- list(
    full_fixed            = function(d) fit_full_fixed(d),
    full_random           = function(d) fit_full_random(d),
    hybrid                = function(d) fit_hybrid(d),
    direct_clustering     = function(d) fit_direct_clustering(d, J = J_direct),
    spherical_clustering  = function(d) fit_spherical_clustering(d, J = J_spherical,
                                         alpha_hotelling = alpha_hotelling),
    autoregressive        = function(d) fit_autoregressive(d),
    hybrid_autoregressive = function(d) fit_hybrid_autoregressive(d)
  )

  lapply(names(fitters), function(nm) {
    message("Fitting: ", nm, " ...")
    tryCatch(fitters[[nm]](df), error = function(e) {
      warning("fit_all_models: '", nm, "' failed — ", e$message)
      NULL
    })
  }) |> setNames(names(fitters))
}
