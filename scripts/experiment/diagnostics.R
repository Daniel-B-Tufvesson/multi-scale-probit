library(coda)

#' Compute Cumulative Gelman-Rubin R-hat Statistic
#'
#' Computes the cumulative Gelman-Rubin R-hat convergence diagnostic for multiple MCMC chains.
#' This tracks how the R-hat statistic improves as more samples are collected.
#'
#' @param chains An mcmc.list object, a list of mcmc objects, or a numeric matrix/array of shape 
#'   (n_chains, n_samples) containing MCMC samples.
#' @param min_samples Minimum number of samples before computing R-hat. Default is 20.
#'
#' @return A numeric vector of cumulative R-hat values, with length equal to n_samples.
#'   Elements 1 to (min_samples - 1) are NA. Elements min_samples onward contain R-hat values.
#'
#' @details
#' The Gelman-Rubin R-hat statistic (also called PSRF, potential scale reduction factor) 
#' compares within-chain and between-chain variance to assess MCMC convergence.
#' R-hat < 1.05 is typically considered good convergence.
#'
#' @examples
#' \dontrun{
#' library(coda)
#' 
#' # Example with synthetic data and list of mcmc objects
#' set.seed(42)
#' chain1 <- mcmc(rnorm(1000))
#' chain2 <- mcmc(rnorm(1000))
#' chain3 <- mcmc(rnorm(1000))
#' chain4 <- mcmc(rnorm(1000))
#' chains_list <- list(chain1, chain2, chain3, chain4)
#' 
#' rhat_cumulative <- cumulative_gelman_rubin_rhat(chains_list, min_samples = 20)
#' plot(rhat_cumulative, type = "l", ylim = c(1, 1.2))
#' }
#'
cumulative_gelman_rubin_rhat <- function(chains, for_every = 100) {
  # Convert to matrix depending on input type
  if (inherits(chains, "mcmc.list")) {
    chains <- as.matrix(chains)
  } else if (is.list(chains)) {
    # Handle list of mcmc objects
    chains <- do.call(rbind, lapply(chains, function(x) {
      if (inherits(x, "mcmc")) {
        as.numeric(x)
      } else {
        as.numeric(x)
      }
    }))
  } else {
    chains <- as.matrix(chains)
  }

  chains = t(chains)
  
  # Ensure numeric type
  storage.mode(chains) <- "numeric"
  
  # Handle 1D case by reshaping
  if (is.null(nrow(chains)) || nrow(chains) == 1) {
    chains <- t(chains)
  }
  
  n_chains <- nrow(chains)
  n_samples <- ncol(chains)
  
  # Initialize output vector with NAs
  n_rhats = n_samples / for_every
  rhat_values <- rep(NA_real_, n_rhats)
  
  # Compute R-hat for each cumulative sample size
  for (t in 1:n_rhats) {
    up_to = t * for_every
    chains_subset <- chains[, 1:up_to, drop = FALSE]
    rhat_values[t] <- gelman_rubin_rhat_single_param(chains_subset)
  }
  
  return(rhat_values)
}


#' Gelman-Rubin R-hat for a Single Parameter
#'
#' Internal helper function to compute R-hat for a single parameter
#' across multiple chains.
#'
#' @param chains Numeric matrix of shape (n_chains, n_samples)
#'
#' @return Numeric scalar representing the R-hat value
#'
#' @keywords internal
gelman_rubin_rhat_single_param <- function(chains) {
  n_chains <- nrow(chains)
  n_samples <- ncol(chains)
  
  # Compute mean of each chain
  chain_means <- rowMeans(chains)
  
  # Compute overall mean
  overall_mean <- mean(chain_means)
  
  # Between-chain variance: B
  # B = n / (m - 1) * sum((chain_i_mean - overall_mean)^2)
  B <- (n_samples / (n_chains - 1)) * sum((chain_means - overall_mean)^2)
  
  # Within-chain variance: W
  # W = 1/m * sum(var_i) where var_i is variance within chain i
  chain_vars <- apply(chains, 1, var)
  W <- mean(chain_vars)
  
  # Estimated variance of the stationary distribution
  # var_hat = ((n - 1) / n) * W + (1 / n) * B
  var_hat <- ((n_samples - 1) / n_samples) * W + (1 / n_samples) * B
  
  # R-hat (Potential Scale Reduction Factor)
  # R-hat = sqrt(var_hat / W)
  rhat <- sqrt(var_hat / W)
  
  return(rhat)
}


#' Compute Effective Sample Size (ESS)
#'
#' Computes the effective sample size for each parameter, accounting for autocorrelation
#' in the MCMC samples. ESS quantifies the number of independent samples equivalent to
#' the autocorrelated samples.
#'
#' @param samples A numeric matrix of shape (n_samples, n_parameters) where rows are 
#'   MCMC samples and columns are parameters.
#' @param method The method for computing ESS. Options are:
#'   - "tau" (default): Uses integrated autocorrelation time (Sokal 1997)
#'   - "acf": Uses autocorrelation function with automatic cutoff
#' @param max_lag Maximum lag to consider for autocorrelation. Default is NULL,
#'   which uses min(nrow(samples) / 2, 5 * ceiling(tau)) where tau is the 
#'   autocorrelation time.
#'
#' @return A numeric vector of ESS values for each parameter, with length equal
#'   to the number of columns in the input matrix.
#'
#' @details
#' The effective sample size is computed as:
#'   ESS = N / (1 + 2 * sum of autocorrelations)
#' where N is the total number of samples and the sum is over all lags with 
#' significant autocorrelation.
#'
#' @examples
#' \dontrun{
#' # Example with synthetic autocorrelated data
#' set.seed(42)
#' n_samples <- 10000
#' n_params <- 3
#' 
#' # Create samples with different autocorrelation patterns
#' samples <- matrix(nrow = n_samples, ncol = n_params)
#' samples[, 1] <- arima.sim(model = list(ar = 0.1), n = n_samples)
#' samples[, 2] <- arima.sim(model = list(ar = 0.5), n = n_samples)
#' samples[, 3] <- arima.sim(model = list(ar = 0.9), n = n_samples)
#' 
#' ess <- effective_sample_size(samples)
#' print(ess)
#' # Parameters with high autocorrelation should have lower ESS
#' }
#'
effective_sample_size <- function(samples, method = "tau", max_lag = NULL) {
  # Ensure matrix format
  samples <- as.matrix(samples)
  storage.mode(samples) <- "numeric"
  
  n_samples <- nrow(samples)
  n_params <- ncol(samples)
  
  # Initialize output vector
  ess_values <- numeric(n_params)
  
  if (method == "tau") {
    # Method using integrated autocorrelation time
    for (i in seq_len(n_params)) {
      ess_values[i] <- ess_from_tau(samples[, i], max_lag = max_lag)
    }
  } else if (method == "acf") {
    # Method using ACF with automatic cutoff
    for (i in seq_len(n_params)) {
      ess_values[i] <- ess_from_acf(samples[, i], max_lag = max_lag)
    }
  } else {
    stop("method must be 'tau' or 'acf'")
  }
  
  return(ess_values)
}


#' Compute ESS using Integrated Autocorrelation Time
#'
#' @param x Numeric vector of samples
#' @param max_lag Maximum lag to consider
#'
#' @return Numeric scalar, the effective sample size
#'
#' @keywords internal
ess_from_tau <- function(x, max_lag = NULL) {
  n <- length(x)
  
  # Center the data
  x_centered <- x - mean(x)
  
  # Compute variance
  var_x <- mean(x_centered^2)
  
  if (var_x <= 0) {
    return(n)
  }
  
  # Compute autocorrelation at different lags
  if (is.null(max_lag)) {
    # Start with a generous upper bound
    max_lag_init <- min(n %/% 2, 5000)
  } else {
    max_lag_init <- max_lag
  }
  
  acf_values <- numeric(max_lag_init)
  for (lag in seq_len(max_lag_init)) {
    acf_values[lag] <- mean(x_centered[1:(n - lag)] * x_centered[(lag + 1):n]) / var_x
  }
  
  # Find automatic cutoff point where ACF becomes small
  # Use criterion: cutoff where acf < 0.05 and stays small
  cutoff <- 1
  for (lag in seq_along(acf_values)) {
    if (acf_values[lag] < 0.05) {
      cutoff <- lag
      break
    }
  }
  
  # If no cutoff found, use max_lag
  if (cutoff == 1) {
    cutoff <- max_lag_init
  }
  
  # Integrated autocorrelation time (Sokal 1997)
  # tau = 0.5 + sum of autocorrelations from lag 1 onwards
  tau <- 0.5 + sum(acf_values[1:min(cutoff, max_lag_init)])
  
  # Ensure tau is at least 0.5
  tau <- max(tau, 0.5)
  
  # ESS = N / (2 * tau)
  ess <- n / (2 * tau)
  
  return(ess)
}


#' Compute ESS using Autocorrelation Function
#'
#' @param x Numeric vector of samples
#' @param max_lag Maximum lag to consider
#'
#' @return Numeric scalar, the effective sample size
#'
#' @keywords internal
ess_from_acf <- function(x, max_lag = NULL) {
  n <- length(x)
  
  # Center the data
  x_centered <- x - mean(x)
  
  # Compute variance
  var_x <- mean(x_centered^2)
  
  if (var_x <= 0) {
    return(n)
  }
  
  # Determine max_lag
  if (is.null(max_lag)) {
    max_lag <- min(n %/% 2, 5000)
  }
  
  # Compute ACF using stats::acf (extract values)
  acf_obj <- stats::acf(x, lag.max = max_lag, plot = FALSE)
  acf_vals <- as.numeric(acf_obj$acf)[-1]  # Remove lag 0
  
  # Find cutoff: where ACF becomes negligible
  # Using criterion: first lag where |ACF| < 0.05 and stays small
  cutoff <- length(acf_vals)
  threshold <- 0.05
  
  for (lag in seq_along(acf_vals)) {
    if (abs(acf_vals[lag]) < threshold) {
      # Check if next few values also small (avoid noise)
      next_lags <- acf_vals[lag:min(lag + 2, length(acf_vals))]
      if (all(abs(next_lags) < threshold)) {
        cutoff <- lag - 1
        if (cutoff < 1) cutoff <- 1
        break
      }
    }
  }
  
  # Compute ESS as N / (1 + 2 * sum of autocorrelations)
  sum_acf <- sum(acf_vals[1:cutoff])
  tau <- 0.5 + sum_acf
  tau <- max(tau, 0.5)
  
  ess <- n / (2 * tau)
  
  return(ess)
}