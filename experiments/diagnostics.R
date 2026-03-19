library(coda)

source("experiments/plotting.R")

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

#   chains = t(chains)
  
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

#' Compute First Index of Convergence
#'
#' Determines the first sample index where MCMC convergence is achieved and maintained.
#' Convergence is defined as R-hat <= 1.01, with no subsequent violations.
#'
#' @param rhat_cumulative A numeric vector of cumulative R-hat values, typically
#'   produced by \code{\link{cumulative_gelman_rubin_rhat}}.
#' @param threshold The R-hat threshold for convergence. Default is 1.01.
#'
#' @return An integer representing the first index where R-hat <= threshold and
#'   all subsequent values remain <= threshold. Returns NA_integer_ if convergence
#'   is never achieved (i.e., some later values exceed the threshold).
#'
#' @details
#' This function searches for the first index where R-hat drops to the convergence
#' threshold and remains below (or equal to) that threshold for all remaining samples.
#' This is stricter than simply finding the first R-hat value below the threshold,
#' as it requires sustained convergence.
#'
#' @examples
#' \dontrun{
#' # Example with synthetic R-hat values showing convergence
#' rhat_values <- c(NA, NA, 1.5, 1.3, 1.15, 1.08, 1.05, 1.03, 1.02, 1.01, 
#'                  1.00, 0.99, 1.00, 0.99, 0.98)
#' convergence_idx <- first_convergence_index(rhat_values)
#' print(convergence_idx)  # Should return 9 (first index with rhat <= 1.01 
#'                         # where all later values also <= 1.01)
#' }
#'
#' @export
first_convergence_index <- function(rhat_cumulative, threshold = 1.01, step_size = 1) {
  # Ensure vector format
  rhat_cumulative <- as.numeric(rhat_cumulative)
  
  n <- length(rhat_cumulative)
  
  if (n == 0) {
    return(NA_integer_)
  }
  
  # Find all indices where R-hat <= threshold
  converged_indices <- which(rhat_cumulative <= threshold)
  
  if (length(converged_indices) == 0) {
    # Never converged
    return(NA_integer_)
  }
  
  # For each potential convergence index, check if all subsequent values are <= threshold
  for (idx in converged_indices) {
    # Check all values from idx onwards
    subsequent_values <- rhat_cumulative[idx:n]
    
    # Check if all subsequent values (excluding NAs) are <= threshold
    subsequent_valid <- subsequent_values[!is.na(subsequent_values)]
    
    if (all(subsequent_valid <= threshold)) {
      return(as.integer(idx * step_size))
    }
  }
  
  # If no index satisfies the condition, return NA
  return(NA_integer_)
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

run_full_diagnostic <- function(all_runs) {
    n_splits<- length(all_runs)
    n_samples <- nrow(all_runs[[1]]$beta)
    n_beta_params <- ncol(all_runs[[1]]$beta)
    n_targets <- length(all_runs[[1]]$gammas)
    n_gammas <- lapply(1:n_targets, function(target_idx) ncol(all_runs[[1]]$gammas[[target_idx]]))
    n_gamma_params <- sum(unlist(n_gammas))

    # Unpack all beta samples into a single matrix for each parameter.
    all_betas <- lapply(1:n_beta_params, function(j){
        mat <- do.call(cbind, lapply(all_runs, function(run) run$beta[, j, drop = FALSE]))
        return(t(mat))
    })

    print("dims")
    print(dim(all_betas[[1]]))

    # Unpack all gamma samples into a single matrix for each parameter and target.
    all_gammas <- lapply(1:n_targets, function(target_idx){
        lapply(1:n_gammas[[target_idx]], function(gamma_idx) {
            mat <- do.call(cbind, lapply(all_runs, function(run) run$gammas[[target_idx]][, gamma_idx, drop = FALSE]))
            return(t(mat))
        })
    })

    # Discard the 50% first burn-in samples.
    all_sampled_betas <- lapply(all_betas, function(beta_mat) beta_mat[, (n_samples / 2 + 1):n_samples, drop = FALSE])
    all_sampled_gammas <- lapply(all_gammas, function(gamma_list) {
        lapply(gamma_list, function(gamma_mat) gamma_mat[, (n_samples / 2 + 1):n_samples, drop = FALSE])
    })


    ###################################################################################################
    # Compute the cumulative R-hat for full samples samples.
    rhat_for_every <- 100
    cum_rhat_burnin_beta <- list()
    for (j in 1:n_beta_params) {
        cum_rhat_burnin_beta[[j]] <- cumulative_gelman_rubin_rhat(all_betas[[j]], for_every = rhat_for_every)
    }

    cum_rhat_burnin_gammas <- list()
    for (j in 1:n_targets) {
        cum_rhat <- list()
        for (i in 1:n_gammas[[j]]) {
            cum_rhat[[i]] <- cumulative_gelman_rubin_rhat(all_gammas[[j]][[i]], for_every = rhat_for_every)
        }
        cum_rhat_burnin_gammas[[j]] <- cum_rhat
    }

    # Compute convergence time for burn-in beta samples.
    burnin_convergence_time_beta <- list()
    n_burnin_beta_convergences <- 0 # Number of betas that converged.
    for (j in 1:n_beta_params) {
        convergence_time <- first_convergence_index(cum_rhat_burnin_beta[[j]], step_size = rhat_for_every)
        burnin_convergence_time_beta[[j]] <- convergence_time

        if (!is.na(convergence_time)) {
            n_burnin_beta_convergences <- n_burnin_beta_convergences + 1
        }
    }
    # Compute convergence time for burn-in beta samples.
    burnin_convergence_time_gammas <- list()
    n_burnin_gamma_convergences = 0 # Number of gammas that converged.
    for (j in 1:n_targets) {
        convergence_times <- list()
        for (i in 1:n_gammas[[j]]) {
            convergence_time <- first_convergence_index(cum_rhat_burnin_gammas[[j]][[i]], step_size = rhat_for_every)
            convergence_times[[i]] <- convergence_time

            if (!is.na(convergence_time)) {
                n_burnin_gamma_convergences <- n_burnin_gamma_convergences + 1
            }
        }
        burnin_convergence_time_gammas[[j]] <- convergence_times
    }

    # Compute median convergence time for burn-in betas and gammas.
    median_convergence_time_beta <- median(unlist(burnin_convergence_time_beta), na.rm = TRUE)
    median_convergence_time_gammas <- median(unlist(burnin_convergence_time_gammas), na.rm = TRUE)

    ###################################################################################################
    # Compute regular R-hat for sampled values.
    rhat_sampled_betas <- list()
    n_sampled_betas_convergences = 0 # Number of sampled betas that converged.
    for (j in 1:n_beta_params) {
        rhat_sampled_betas[[j]] <- gelman_rubin_rhat_single_param(all_sampled_betas[[j]])

        if (rhat_sampled_betas[[j]] <= 1.01) {
            n_sampled_betas_convergences <- n_sampled_betas_convergences + 1
        }
    }
    rhat_sampled_gammas <- list()
    n_sampled_gammas_convergences = 0 # Number of sampled gammas that converged.
    for (j in 1:n_targets) {
        rhat_gammas <- list()
        for (i in 1:n_gammas[[j]]) {
            rhat_gammas[[i]] <- gelman_rubin_rhat_single_param(all_sampled_gammas[[j]][[i]])

            if (rhat_gammas[[i]] <= 1.01) {
                n_sampled_gammas_convergences <- n_sampled_gammas_convergences + 1
            }
        }
        rhat_sampled_gammas[[j]] <- rhat_gammas
    }

    ##################################################################################################
    # Compute min ESS
    beta_min_ess <- list() # One for each split.
    beta_min_ess_per_second <- list() 
    beta_min_ess_per_iteration <- list()
    for (j in 1:n_splits) {
        ess_values <- effective_sample_size(all_runs[[j]]$beta, method = "acf")
        time_seconds <- all_runs[[j]]$samplingTime
        # ess_values <- effective_sample_size(cv_res$allBetas[[j]], method = "acf")
        # time_seconds <- cv_res$samplingTime[j]
        min_ess <- min(ess_values)
        beta_min_ess[[j]] <- min_ess
        beta_min_ess_per_second[[j]] <- min_ess / time_seconds
        beta_min_ess_per_iteration[[j]] <- min_ess / n_samples
    }

    gamma_min_ess <- list() # One for each split.
    gamma_min_ess_per_second <- list()
    gamma_min_ess_per_iteration <- list()
    for (j in 1:n_splits) {
        gammas_split <- all_runs[[j]]$gammas
        # gammas_split <- cv_res$allGammas[[j]]
        ess_values <- unlist(lapply(gammas_split, function(gamma_mat) effective_sample_size(gamma_mat, method = "acf")))
        # time_seconds <- cv_res$samplingTime[j]
        time_seconds <- all_runs[[j]]$samplingTime
        min_ess <- min(ess_values)
        gamma_min_ess[[j]] <- min_ess
        gamma_min_ess_per_second[[j]] <- min_ess / time_seconds
        gamma_min_ess_per_iteration[[j]] <- min_ess / n_samples
    }

    # Cpmpute median min ESS.
    median_min_ess_beta <- median(unlist(beta_min_ess))
    median_min_ess_gamma <- median(unlist(gamma_min_ess))
    median_min_ess_per_second_beta <- median(unlist(beta_min_ess_per_second))
    median_min_ess_per_second_gamma <- median(unlist(gamma_min_ess_per_second))
    median_min_ess_per_iteration_beta <- median(unlist(beta_min_ess_per_iteration))
    median_min_ess_per_iteration_gamma <- median(unlist(gamma_min_ess_per_iteration))

    return(list(
        n_targets = n_targets,
        n_beta_params = n_beta_params,
        n_gamma_params = n_gamma_params,
        n_gammas = n_gammas,

        n_burnin_beta_convergences = n_burnin_beta_convergences,
        n_burnin_gamma_convergences = n_burnin_gamma_convergences,
        n_sampled_betas_convergences = n_sampled_betas_convergences,
        n_sampled_gammas_convergences = n_sampled_gammas_convergences,
        median_convergence_time_beta = median_convergence_time_beta,
        median_convergence_time_gammas = median_convergence_time_gammas,

        median_min_ess_beta = median_min_ess_beta,
        median_min_ess_gamma = median_min_ess_gamma,
        median_min_ess_per_second_beta = median_min_ess_per_second_beta,
        median_min_ess_per_second_gamma = median_min_ess_per_second_gamma,
        median_min_ess_per_iteration_beta = median_min_ess_per_iteration_beta,
        median_min_ess_per_iteration_gamma = median_min_ess_per_iteration_gamma,

        cum_rhat_burnin_beta = cum_rhat_burnin_beta,
        cum_rhat_burnin_gammas = cum_rhat_burnin_gammas,
        burnin_convergence_time_beta = burnin_convergence_time_beta,
        burnin_convergence_time_gammas = burnin_convergence_time_gammas,
        rhat_sampled_betas = rhat_sampled_betas,
        rhat_sampled_gammas = rhat_sampled_gammas,
        beta_min_ess = beta_min_ess,
        gamma_min_ess = gamma_min_ess,
        beta_min_ess_per_second = beta_min_ess_per_second,
        gamma_min_ess_per_second = gamma_min_ess_per_second,
        beta_min_ess_per_iteration = beta_min_ess_per_iteration,
        gamma_min_ess_per_iteration = gamma_min_ess_per_iteration
    ))
}

# Some helper functions.
is_rhat_ok <- function(rhat) {
    if (rhat <= 1.01) {
        return("✅")
    } else if (rhat <= 1.05) {
        return("🟡")
    } else {
        return("❌")
    }
}

to_percent <- function(value) {
    return(paste0(round(value * 100, 2), "%"))
}

print_diagnostic_report <- function(diagnostic_results) {
    dr <- diagnostic_results
    cat("Diagnostic Report for Gibbs Sampling:\n")
    cat("1. Cumulative R-hat Analysis for Burn-in Beta:\n")
    cat("   - ", dr$n_burnin_beta_convergences, "/", dr$n_beta_params, "(",
        to_percent(dr$n_burnin_beta_convergences / dr$n_beta_params), ") burn-in betas converged.\n")
    for (j in 1:dr$n_beta_params) {
        final_rhat <- tail(dr$cum_rhat_burnin_beta[[j]], n = 1)
        converged_at <- dr$burnin_convergence_time_beta[[j]]
        rhat_ok <- is_rhat_ok(final_rhat)
        cat(paste0("   - ", rhat_ok, " Beta ", j, " converged at ", converged_at, " r-hat = ", final_rhat, "\n"))
    }

    cat("2. Cumulative R-hat Analysis for Burn-in Gammas:\n")
    cat(paste0("      - ", dr$n_burnin_gamma_convergences, "/", dr$n_gamma_params, " (",
            to_percent(dr$n_burnin_gamma_convergences / dr$n_gamma_params), ") burn-in gammas converged.\n"))
    for (j in 1:dr$n_targets) {
        cat(paste0("   - Gammas for target ", j, ":\n"))
        for (i in 1:dr$n_gammas[[j]]) {
            final_rhat <- tail(dr$cum_rhat_burnin_gammas[[j]][[i]], n = 1)
            converged_at <- dr$burnin_convergence_time_gammas[[j]][[i]]
            rhat_ok <- is_rhat_ok(final_rhat)
            cat(paste0("      - ", rhat_ok, " Gamma ", i, " converged at ", converged_at, " r-hat = ", final_rhat, "\n"))
        }
    }
    cat("3. R-hat Analysis for Sampled Betas\n")
    for (j in 1:dr$n_beta_params) {
        final_rhat <- dr$rhat_sampled_betas[[j]]
        rhat_ok <- is_rhat_ok(final_rhat)
        cat(paste0("   - ", rhat_ok, " Beta ", j, ": R-hat = ", round(final_rhat, 4), "\n"))
    }

    cat("4. R-hat Analysis for Sampled Gammas\n")
    for (j in 1:dr$n_targets) {
        cat(paste0("   - Gammas for target ", j, ":\n"))
        for (i in 1:dr$n_gammas[[j]]) {
            final_rhat <- dr$rhat_sampled_gammas[[j]][[i]]
            rhat_ok <- is_rhat_ok(final_rhat)
            cat(paste0("      - ", rhat_ok, " Gamma ", i, ": R-hat = ", round(final_rhat, 4), "\n"))
        }
    }

    cat("5. ESS Analysis for joint Betas and Gammas\n")
    cat("   - Median minimum ESS for sampled betas: ", round(dr$median_min_ess_beta, 4), " ESS\n")
    cat("   - Median minimum ESS for sampled gammas: ", round(dr$median_min_ess_gamma, 4), " ESS\n")
    cat("   - Median minimum ESS/s for sampled betas: ", round(dr$median_min_ess_per_second_beta, 4), " ESS/s\n")
    cat("   - Median minimum ESS/s for sampled gammas: ", round(dr$median_min_ess_per_second_gamma, 4), " ESS/s\n")
    cat("   - Median minimum ESS/iteration for sampled betas: ", round(dr$median_min_ess_per_iteration_beta, 4), " ESS/iteration\n")
    cat("   - Median minimum ESS/iteration for sampled gammas: ", round(dr$median_min_ess_per_iteration_gamma, 4), " ESS/iteration\n")

    cat("6. Summary\n")
    cat("   - ", dr$n_burnin_beta_convergences, "/", dr$n_beta_params, "(",
        to_percent(dr$n_burnin_beta_convergences / dr$n_beta_params), ") burn-in betas converged.\n")
    cat("   - ", dr$n_burnin_gamma_convergences, "/", dr$n_gamma_params, "(",
        to_percent(dr$n_burnin_gamma_convergences / dr$n_gamma_params), ") burn-in gammas converged.\n")
    cat("   - ", dr$n_sampled_betas_convergences, "/", dr$n_beta_params, "(",
        to_percent(dr$n_sampled_betas_convergences / dr$n_beta_params), ") sampled betas converged.\n")
    cat("   - ", dr$n_sampled_gammas_convergences, "/", dr$n_gamma_params, "(",
        to_percent(dr$n_sampled_gammas_convergences / dr$n_gamma_params), ") sampled gammas converged.\n")
    cat("   - Median convergence time for burn-in betas: ", dr$median_convergence_time_beta, " samples\n")
    cat("   - Median convergence time for burn-in gammas: ", dr$median_convergence_time_gammas, " samples\n")
}

plot_diagnostic_plots <- function(diagnostic_results) {
    dr <- diagnostic_results
    # Plot convergence time distribution for burn-in betas and gammas.
    # print(density_plot(unlist(dr$burnin_convergence_time_beta),
    #             title = "Density of Convergence Times for Burn-in Betas", 
    #             xlabel = "Convergence Time (samples)"))
    # print(density_plot(unlist(dr$burnin_convergence_time_gammas),
    #             title = "Density of Convergence Times for Burn-in Gammas", 
    #             xlabel = "Convergence Time (samples)"))
    
    # Plot distributions of ESS/per second.
    print(density_plot(unlist(dia_res$beta_min_ess_per_second),
                title = "Density of Minimum ESS/s for Sampled Betas", 
                xlabel = "Minimum ESS/s"))

    print(density_plot(unlist(dia_res$gamma_min_ess_per_second),
                title = "Density of Minimum ESS/s for Sampled Gammas", 
                xlabel = "Minimum ESS/s"))

    # Plot distributions of ESS/iteration.
    print(density_plot(unlist(dia_res$beta_min_ess_per_iteration),
                title = "Density of Minimum ESS/iteration for Sampled Betas", 
                xlabel = "Minimum ESS/iteration"))

    print(density_plot(unlist(dia_res$gamma_min_ess_per_iteration),
                title = "Density of Minimum ESS/iteration for Sampled Gammas", 
                xlabel = "Minimum ESS/iteration"))
}


# Plot rhat values.
plot_cum_rhat <- function(cum_rhat, rhat_for_every){
    x_vals <- seq_along(cum_rhat) * rhat_for_every
    plot(x_vals, cum_rhat, type = "l", ylim = c(1, 1.2), 
        xlab = "Cumulative Samples", ylab = "R-hat", main = "Cumulative R-hat for Burn-in Beta")
    abline(h = 1.05, col = "red", lty = 2, lwd = 2)  # Not converged.
    abline(h = 1.01, col = "orange", lty = 2, lwd = 2)  # Maybe converged.
}