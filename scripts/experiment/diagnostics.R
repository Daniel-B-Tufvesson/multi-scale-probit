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