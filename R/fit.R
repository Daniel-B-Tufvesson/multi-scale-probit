library(coda)
# library(callr)

tune_mspm <- function(
    data,
    iterations,
    ...,
    target_acceptance_rate = 0.234,
    target_epsilon = 0.05,
    stop_early = FALSE,
    window_size = 25,
    seed = NA,
    mean_prior = NULL,
    prec_prior = NULL,
    proposal_variance_initial = NULL,
    beta_initial = NULL,
    gamma_initial = NULL,
    verbose = 0
) {
    .validate_data(data)

    nlevels <- get_n_levels(data)
    ntargets <- get_n_targets(data)
    npredictors <- get_n_predictors(data)

    # Set the seed.
    if (is.na(seed)) {
        seed <- sample.int(.Machine$integer.max, 1)
    }
    set.seed(seed)

    # Set default priors.
    if (is.null(mean_prior)) {
        mean_prior <- rep(0, npredictors)
    }
    if (is.null(prec_prior)) {
        prec_prior <- .create_prec_prior(npredictors)
    }

    # Set starting values for gamma and beta if not provided
    if (is.null(gamma_initial)) {
        gamma_initial <- .create_inital_gammas(ntargets, nlevels)
    }
    if (is.null(beta_initial)) {
        beta_initial <- rep(0, npredictors)
    }

    # Set tuning parameter for the sampler
    proposal_variance_initial <- .reshape_proposal_variance_gibbs(proposal_variance_initial, ntargets)

    tune_results <- cpp_hprobit_tune(
        data$Xlist,
        data$ylist,
        mean_prior,
        prec_prior,
        nlevels,
        gamma_initial,
        beta_initial,
        proposal_variance_initial,
        target_acceptance_rate,
        target_epsilon,
        stop_early,
        iterations,
        window_size,
        seed,
        verbose
    )

    # Tmp: start as subprocess for more robust development.
    # tune_results <- tryCatch({callr::r(
    #     function(data, mean_prior, prec_prior, nlevels, gamma_initial, beta_initial, start_tune, 
    #              target_acceptance_rate, target_epsilon, stop_early, iterations, window_size, seed, 
    #              verbose) {
    #         devtools::load_all()
    #         cpp_hprobit_tune(
    #             data$Xlist,
    #             data$ylist,
    #             mean_prior,
    #             prec_prior,
    #             nlevels,
    #             gamma_initial,
    #             beta_initial,
    #             start_tune,
    #             target_acceptance_rate,
    #             target_epsilon,
    #             stop_early,
    #             iterations,
    #             window_size,
    #             seed,
    #             verbose
    #         )
    #     },
    #     args = list(data, mean_prior, prec_prior, nlevels, gamma_initial, beta_initial, start_tune, 
    #                 target_acceptance_rate, target_epsilon, stop_early, iterations, window_size, 
    #                 seed, verbose),
    #     show = verbose > 0
    # )}, error = function(e) {
    #     message("Error in cpp_hprobit: ", e$message)
    #     if (!is.null(e$stdout)) {
    #         cat("---- STDOUT ----\n")
    #         cat(e$stdout, sep = "\n")
    #     }
    #     if (!is.null(e$stderr)) {
    #         cat("---- STDERR ----\n")
    #         cat(e$stderr, sep = "\n")
    #     }
    #     stop(e)
    # })

    # Return tuning results.
    return(new_mspm_tune_results(
        data_spec = get_data_spec(data),
        proposal_variance = tune_results$proposal_variance,
        target_acceptance_rate = target_acceptance_rate,
        acceptance_rates = tune_results$acceptance_rates,
        max_iterations = iterations,
        final_iteration = tune_results$final_iteration,
        seed = seed,
        call = match.call()
    ))
}

# Fit a multi-scale probit model (MSPM) to the data.  
# 
# Arguments:
# data: An MSPM data structure holding the data.
# burnin: Number of burn-in iterations.
# ndraws: Number of posterior draws to collect.
# thin: Thinning interval for MCMC sampling.
# meanPrior: Prior mean for regression coefficients.
# precPrior: Prior precision for regression coefficients.
# tune: Tuning parameter for the sampler.
# adapt_tune: Whether to adapt the tuning parameter during burn-in.
# tuneWindowSize: Window size for computing acceptance rates for tuning adaptation.
# targetAcceptanceRate: Target acceptance rate for tuning adaptation. 
# seed: Random seed for reproducibility.
# beta.initial: Initial values for regression coefficients.
# gamma.initial: Initial values for threshold parameters.
# verbose: Verbosity level for output.
# 
# Returns:
# An object of class 'mspm' containing the fitted model.
fit_mspm <- function(
    data, 
    burnin,
    ndraws,
    thin,
    ...,
    mean_prior = NULL,
    prec_prior = NULL,
    proposal_variance = NULL,
    seed = NA,
    beta_start = NULL,
    gamma_start = NULL,
    verbose = 0
) {
    .validate_data(data)

    nlevels <- get_n_levels(data)
    ntargets <- get_n_targets(data)
    npredictors <- get_n_predictors(data)

    # Set the seed.
    if (is.na(seed)) {
        seed <- sample.int(.Machine$integer.max, 1)
    }
    set.seed(seed)

    # Set default priors.
    if (is.null(mean_prior)) {
        mean_prior <- rep(0, npredictors)
    }
    if (is.null(prec_prior)) {
        prec_prior <- .create_prec_prior(npredictors)
    }

    # Set starting values for gamma and beta if not provided
    if (is.null(gamma_start)) {
        gamma_start <- .create_inital_gammas(ntargets, nlevels)
    }
    else if (length(gamma_start) == ntargets) {
        # Append edge gammas if only inner gammas are provided.
        gamma_start <- .add_edge_gammas(gamma_start)
    }
    else if (length(gamma_start) != ntargets + 2) {
        # No edge gammas provided, but length does not match number of targets.
        stop("Length of gamma_start must match number of targets.")
    }

    # Set starting values for beta if not provided.
    if (is.null(beta_start)) {
        beta_start <- rep(0, npredictors)
    }

    # Set proposal variance parameter for the sampler
    proposal_variance <- .reshape_proposal_variance_gibbs(proposal_variance, ntargets)

    # Run CPP backend sampler.
    sim <- cpp_hprobit(
        data$Xlist,
        data$ylist,
        mean_prior,
        prec_prior,
        nlevels,
        gamma_start,
        beta_start,
        proposal_variance,
        ndraws,
        burnin,
        thin,
        seed,
        verbose
    )

    # Tmp: start as subprocess for more robust development.
    # sim <- tryCatch({callr::r(
    #     function(data, meanPrior, precPrior, nlevels, gamma.initial, beta.initial, tune, 
    #              ndraws, burnin, thin, seed, verbose, saveBurninSamples, adapt_tune, tuneWindowSize, 
    #              targetAcceptanceRate) {
    #         devtools::load_all()
    #         cpp_hprobit(
    #             data$Xlist,
    #             data$ylist,
    #             meanPrior,
    #             precPrior,
    #             nlevels,
    #             gamma.initial,
    #             beta.initial,
    #             tune,
    #             adapt_tune,
    #             tuneWindowSize,
    #             targetAcceptanceRate,
    #             ndraws,
    #             burnin,
    #             thin,
    #             saveBurninSamples,
    #             seed,
    #             verbose
    #         )
    #     },
    #     args = list(data, meanPrior, precPrior, nlevels, gamma.initial, beta.initial, 
    #                 tune, ndraws, burnin, thin, seed, verbose, saveBurninSamples, adapt_tune,
    #                 tuneWindowSize, targetAcceptanceRate),
    #     show = verbose > 0
    # )}, error = function(e) {
    #     message("Error in cpp_hprobit: ", e$message)
    #     if (!is.null(e$stdout)) {
    #         cat("---- STDOUT ----\n")
    #         cat(e$stdout, sep = "\n")
    #     }
    #     if (!is.null(e$stderr)) {
    #         cat("---- STDERR ----\n")
    #         cat(e$stderr, sep = "\n")
    #     }
    #     stop(e)
    # })

    # Format and store MCMC results
    colnames(sim$storebeta) <- c(data$predictorNames)
    beta <- mcmc(sim$storebeta, start = burnin + 1, thin = thin)
    gammas <- list()
    for (i in 1:ntargets) {
        gammas[[i]] <- mcmc(sim$storegamma[[i]], start = burnin + 1, thin = thin)
    }

    # Return fitted model.
    new_mspm(
        data_spec = get_data_spec(data),
        beta = beta,
        gammas = gammas,
        meanPrior = mean_prior,
        precPrior = prec_prior,
        proposal_variance = proposal_variance,
        acceptanceRate = sim$acceptance_rate,
        seed = seed,
        ndraws = ndraws / thin,
        ndrawsNoThin = ndraws,
        thin = thin,
        burnin = burnin,
        samplingTime = sim$sampling_time,
        burninTime = sim$burnin_time,
        nlikelihood_calls = sim$nlikelihood_calls,
        call = match.call()
    )  
}



tune_mspm_pt <- function(
    data,
    iterations,
    ntemperatures,
    ...,
    tune_proposal_variance = TRUE,
    tune_temperature_ladder = TRUE,
    target_acceptance_rate = 0.234,
    target_acceptance_epsilon = 0.05,
    target_temp_swap_accept_rate = 0.3,
    target_temp_swap_accept_epsilon = 0.05,
    stop_early = FALSE,
    proposal_window_size = 25,
    temperature_window_size = 100,
    temperature_window_growth_factor = 1.2,
    temperature_ladder_learning_rate = 0.01,
    seed = NA,
    mean_prior = NULL,
    prec_prior = NULL,
    proposal_variance_initial = NULL,
    beta_initial = NULL,
    gamma_initial = NULL,
    inv_temperature_ladder = NULL,
    complete_param_swapping = TRUE,
    verbose = 0
) {
    .validate_data(data)

    nlevels <- get_n_levels(data)
    ntargets <- get_n_targets(data)
    npredictors <- get_n_predictors(data)

    # Set the seed.
    if (is.na(seed)) {
        seed <- sample.int(.Machine$integer.max, 1)
    }
    set.seed(seed)

    # Set default priors.
    if (is.null(mean_prior)) {
        mean_prior <- rep(0, npredictors)
    }
    if (is.null(prec_prior)) {
        prec_prior <- .create_prec_prior(npredictors)
    }

    # Set starting values for gamma and beta if not provided
    if (is.null(gamma_initial)) {
        gamma_initial <- .create_inital_gammas(ntargets, nlevels)
    }
    if (is.null(beta_initial)) {
        beta_initial <- rep(0, npredictors)
    }

    # Set tuning parameter for the sampler
    proposal_variance_initial <- .reshape_proposal_variance_pt(proposal_variance_initial, 
        ntargets, ntemperatures)

    inv_temperature_ladder <- .prepare_temperature_ladder(inv_temperature_ladder, ntemperatures)

    # Call backend tuner.
    tune_results <- cpp_hprobit_tune_pt(
        data$Xlist,
        data$ylist,
        mean_prior,
        prec_prior,
        nlevels,
        gamma_initial,
        beta_initial,
        proposal_variance_initial,
        tune_proposal_variance,
        target_acceptance_rate,
        target_acceptance_epsilon,
        proposal_window_size,
        inv_temperature_ladder,
        tune_temperature_ladder,
        target_temp_swap_accept_rate,
        target_temp_swap_accept_epsilon,
        temperature_window_size,
        temperature_window_growth_factor,
        temperature_ladder_learning_rate,
        iterations,
        stop_early,
        seed,
        complete_param_swapping,
        verbose
    )

    # Tmp: start as subprocess for more robust development.
    # tune_results <- tryCatch({callr::r(
    #     function(data, mean_prior, prec_prior, nlevels, gamma_initial, beta_initial, proposal_variance_initial, 
    #              target_acceptance_rate, target_acceptance_epsilon, target_temp_swap_accept_rate,
    #              target_temp_swap_accept_epsilon, proposal_window_size, inv_temperature_ladder, 
    #              temperature_window_size, temperature_window_growth_factor, temperature_ladder_learning_rate,
    #              tune_proposal_variance, tune_temperature_ladder,
    #              stop_early, iterations, seed, complete_param_swapping, verbose) {
    #         devtools::load_all()
    #         cpp_hprobit_tune_pt(
    #             data$Xlist,
    #             data$ylist,
    #             mean_prior,
    #             prec_prior,
    #             nlevels,
    #             gamma_initial,
    #             beta_initial,
    #             proposal_variance_initial,
    #             tune_proposal_variance,
    #             target_acceptance_rate,
    #             target_acceptance_epsilon,
    #             proposal_window_size,
    #             inv_temperature_ladder,
    #             tune_temperature_ladder,
    #             target_temp_swap_accept_rate,
    #             target_temp_swap_accept_epsilon,
    #             temperature_window_size,
    #             temperature_window_growth_factor,
    #             temperature_ladder_learning_rate,
    #             iterations,
    #             stop_early,
    #             seed,
    #             complete_param_swapping,
    #             verbose
    #         )
    #     },
    #     args = list(data, mean_prior, prec_prior, nlevels, gamma_initial, beta_initial, proposal_variance_initial, 
    #                 target_acceptance_rate, target_acceptance_epsilon, target_temp_swap_accept_rate,
    #                 target_temp_swap_accept_epsilon, proposal_window_size, inv_temperature_ladder, temperature_window_size,
    #                 temperature_window_growth_factor, temperature_ladder_learning_rate,
    #                 tune_proposal_variance, tune_temperature_ladder,
    #                 stop_early, iterations, seed, complete_param_swapping, verbose),
    #     show = verbose > 0
    # )}, error = function(e) {
    #     message("Error in cpp_hprobit: ", e$message)
    #     if (!is.null(e$stdout)) {
    #         cat("---- STDOUT ----\n")
    #         cat(e$stdout, sep = "\n")
    #     }
    #     if (!is.null(e$stderr)) {
    #         cat("---- STDERR ----\n")
    #         cat(e$stderr, sep = "\n")
    #     }
    #     stop(e)
    # })

    # Return tuning results.
    return(new_mspm_tune_results_pt(
        data_spec = get_data_spec(data),
        proposal_variance = tune_results$proposal_variance,
        proposal_acceptance_rates = tune_results$proposal_acceptance_rates,
        target_acceptance_rate = target_acceptance_rate,
        inv_temperature_ladder = tune_results$inv_temperature_ladder,
        target_temp_swap_rate = target_temp_swap_accept_rate,
        temp_swap_rates = tune_results$temp_swap_rates,
        max_iterations = iterations,
        final_iteration = tune_results$final_iteration,
        seed = seed,
        call = match.call()
    ))
}

#' Fit a multi-scale probit model (MSPM) using parallel tempering.
#'
#' @param data An MSPM data structure holding the data.
#' @param burnin Number of burn-in iterations.
#' @param ndraws Number of posterior draws to collect.
#' @param thin Thinning interval for MCMC sampling.
#' @param ntemperatures Number of temperatures to use in parallel tempering.
#' @param meanPrior Prior mean for regression coefficients.
#' @param precPrior Prior precision for regression coefficients.
#' @param fix.zero Index of the threshold to fix at zero.
#' @param tune Tuning parameter for the sampler.
#' @param temperatureLadder Optional initial temperature ladder. If NULL, a default ladder 
#' will be created.
#' @param target_temp_swap_accept_ratio Target acceptance ratio for temperature swaps, used for 
#' adaptive temperature ladder. If set to -1, the temperature ladder will not be adapted.
#' @param temperature_window_size Window size for computing swap acceptance ratios for adaptive 
#' temperature ladder.
#' @param tempperature_window_growth_factor Growth factor for the window size used in adaptive
#' temperature ladder.
#' @param temperature_ladder_learning_rate Learning rate for updating the temperature ladder in 
#' the adaptive temperature ladder algorithm.
#' @param complete_param_swapping Whether to perform complete swapping of the parameters (beta and 
#' gammas) during temperature swaps, or just swap the gammas.
#' @param seed Random seed for reproducibility.
#' @param beta.initial Initial values for regression coefficients.
#' @param gamma.initial Initial values for threshold parameters.
#' @param verbose Verbosity level for output.
#' @param saveBurninSamples Whether to save the burn-in samples in the returned model object.
#' @return An object of class 'mspm' containing the fitted model.
fit_mspm_pt <- function(
    data, 
    burnin,
    ndraws,
    thin,
    ntemperatures,
    ...,
    mean_prior = NULL,
    prec_prior = NULL,
    proposal_variance = NULL,
    inv_temperature_ladder = NULL,
    complete_param_swapping = TRUE,
    seed = NA,
    beta_start = NULL,
    gamma_start = NULL,
    verbose = 0
) {
    .validate_data(data)

    nlevels <- get_n_levels(data)
    ntargets <- get_n_targets(data)
    npredictors <- get_n_predictors(data)

    if (!is.numeric(ntemperatures) || ntemperatures < 2) {
        stop("ntemperatures must be a numeric value greater than or equal to 2.")
    }

    # Set the seed.
    if (is.na(seed)) {
        seed <- sample.int(.Machine$integer.max, 1)
    }
    set.seed(seed)

    # Set default priors.
    if (is.null(mean_prior)) {
        mean_prior <- rep(0, npredictors)
    }
    if (is.null(prec_prior)) {
        prec_prior <- .create_prec_prior(npredictors)
    }

    # Set starting values for gamma and beta if not provided
    if (is.null(gamma_start)) {
        gamma_start <- .create_inital_gammas(ntargets, nlevels)
    }
    else if (length(gamma_start) == ntargets) {
        # Append edge gammas if only inner gammas are provided.
        gamma_start <- .add_edge_gammas(gamma_start)
    }
    else if (length(gamma_start) != ntargets + 2) {
        # No edge gammas provided, but length does not match number of targets.
        stop("Length of gamma_start must match number of targets.")
    }

    
    if (is.null(beta_start)) {
        beta_start <- rep(0, npredictors)
    }

    # Set tuning parameter for the sampler
    proposal_variance <- .reshape_proposal_variance_pt(proposal_variance, ntargets, ntemperatures)

    # Set default temperature ladder if not provided.
    inv_temperature_ladder <- .prepare_temperature_ladder(inv_temperature_ladder, ntemperatures)


    # Call the backend parallel tempering sampler.
    sim <- cpp_hprobit_pt(
        data$Xlist,
        data$ylist,
        mean_prior,
        prec_prior,
        nlevels,
        gamma_start,
        beta_start,
        proposal_variance,
        inv_temperature_ladder,
        ndraws,
        burnin,
        thin,
        seed,
        complete_param_swapping,
        verbose
    )

    # Call the backend parallel tempering sampler.
    # Tmp: start as subprocess for more robust development.
    # sim <- tryCatch({callr::r(
    #     function(
    #         data, 
    #         mean_prior,
    #         prec_prior,
    #         nlevels,
    #         gamma_start,
    #         beta_start,
    #         proposal_variance,
    #         inv_temperature_ladder,
    #         ndraws,
    #         burnin,
    #         thin,
    #         seed,
    #         complete_param_swapping,
    #         verbose
            
    #     ) {
    #         devtools::load_all()
    #         cpp_hprobit_pt(
    #             data$Xlist,
    #             data$ylist,
    #             mean_prior,
    #             prec_prior,
    #             nlevels,
    #             gamma_start,
    #             beta_start,
    #             proposal_variance,
    #             inv_temperature_ladder,
    #             ndraws,
    #             burnin,
    #             thin,
    #             seed,
    #             complete_param_swapping,
    #             verbose
    #         )
    #     },
    #     args = list(
    #         data, 
    #         mean_prior,
    #         prec_prior,
    #         nlevels,
    #         gamma_start,
    #         beta_start,
    #         proposal_variance,
    #         inv_temperature_ladder,
    #         ndraws,
    #         burnin,
    #         thin,
    #         seed,
    #         complete_param_swapping,
    #         verbose
    #     ),
    #     show = verbose > 0
    # )}, error = function(e) {
    #     message("Error in cpp_hprobit_pt: ", e$message)
    #     if (!is.null(e$stdout)) {
    #         cat("---- STDOUT ----\n")
    #         cat(e$stdout, sep = "\n")
    #     }
    #     if (!is.null(e$stderr)) {
    #         cat("---- STDERR ----\n")
    #         cat(e$stderr, sep = "\n")
    #     }
    #     stop(e)
    # })

    # Format and store MCMC results
    colnames(sim$storebeta) <- c(data$predictorNames)
    beta <- mcmc(sim$storebeta, start = burnin + 1, thin = thin)
    gammas <- list()
    for (i in 1:ntargets) {
        gammas[[i]] <- mcmc(sim$storegamma[[i]], start = burnin + 1, thin = thin)
    }

    # Return fitted model.
    new_mspm_pt(
        data_spec = get_data_spec(data),
        beta = beta,
        gammas = gammas,
        mean_prior = mean_prior,
        prec_prior = prec_prior,
        proposal_variance = proposal_variance,
        seed = seed,
        ndraws = ndraws / thin,
        ndrawsNoThin = ndraws,
        thin = thin,
        inv_temperatures = inv_temperature_ladder,
        burnin = burnin,
        completeSwapping = complete_param_swapping,
        samplingTime = sim$sampling_time,
        burninTime = sim$burnin_time,
        nlikelihood_calls = sim$nlikelihood_calls,
        round_trip_times = sim$round_trip_times,
        call = match.call()
    )  
}



# Helper functions. -------------------------------------------------------------------------------


.validate_data <- function(data) {
    nlevels <- get_n_levels(data)
    ntargets <- get_n_targets(data)
    npredictors <- get_n_predictors(data)

    # Check for empty matrices or vectors in data
    if (is.null(data$Xlist) || length(data$Xlist) == 0) {
        stop("data$Xlist is empty or NULL.")
    }
    if (is.null(data$ylist) || length(data$ylist) == 0) {
        stop("data$ylist is empty or NULL.")
    }
    for (i in seq_along(data$Xlist)) {
        if (is.null(data$Xlist[[i]]) || any(dim(data$Xlist[[i]]) == 0)) {
            stop(sprintf("data$Xlist[[%d]] is empty or has zero dimension.", i))
        }
    }
    for (i in seq_along(data$ylist)) {
        if (is.null(data$ylist[[i]]) || length(data$ylist[[i]]) == 0) {
            stop(sprintf("data$ylist[[%d]] is empty.", i))
        }
    }
    if (npredictors == 0) {
        stop("No predictors found: data$predictorNames is empty.")
    }
    if (is.null(nlevels) || length(nlevels) == 0) {
        stop("data$nlevels is empty or NULL.")
    }
    if (is.null(ntargets) || ntargets == 0) {
        stop("data$ntargets is empty or zero.")
    }
}

.create_inital_gammas <- function(ntargets, nlevels) {
    gamma.initial <- list()
    for (i in 1:ntargets) {
        gamma.initial[[i]] <- rep(0, nlevels[i]+1)
        gamma.initial[[i]][1] <- -300
        gamma.initial[[i]][2:nlevels[i]] <- 0:(nlevels[i]-2)
        gamma.initial[[i]][nlevels[i] + 1] <- 300
    }
    gamma.initial
}

.add_edge_gammas <- function(gammas) {
    # Add edge gammas for the backend function. This should really not be necessary, 
    # but it is what it is for now. :/

    gammas_with_edges <- list()
    for (i in seq_along(gammas)) {
        gamma <- gammas[[i]]
        nlevels <- length(gamma) + 1
        gamma_with_edges <- c(-300, gamma, 300)
        gammas_with_edges[[i]] <- gamma_with_edges
    }
    return (gammas_with_edges)
}

.create_prec_prior <- function(npredictors) {
    if (npredictors > 1) {
        precPrior <- diag(rep(0.1, npredictors))
    }
    else {
        # Special case for single predictor. diag() returns a 0x0 matrix when npredictors=1.
        precPrior <- matrix(0.1, nrow = 1, ncol = 1)
    }
    precPrior
}

.reshape_proposal_variance_gibbs <- function(proposal_variance, ntargets) {
    if (is.null(proposal_variance)) {
        proposal_variance <- 0.05
    }
    if (is.numeric(proposal_variance) && length(proposal_variance) == 1) {
        proposal_variance <- rep(proposal_variance, ntargets)
    }
    return(proposal_variance)
}

.reshape_proposal_variance_pt <- function(proposal_variance, ntargets, ntemperatures) {
    if (is.null(proposal_variance)) {
        proposal_variance <- 0.05
    } 
    if (length(proposal_variance) != ntemperatures) {
        proposal_variance_list <- list()
        for (i in 1:ntemperatures) {
            proposal_variance_list[[i]] <- rep(proposal_variance, ntargets)
        }
        proposal_variance <- proposal_variance_list
    }
    return(proposal_variance)
}

.prepare_temperature_ladder <- function(inv_temperature_ladder, ntemperatures) {
    if (is.null(inv_temperature_ladder)) {
        inv_temperature_ladder <- 1 / (2^((1:ntemperatures) - 1))
    }
    if (length(inv_temperature_ladder) != ntemperatures) {
        stop("Length of inv_temperature_ladder must match ntemperatures.")
    }
    return(inv_temperature_ladder)
}
