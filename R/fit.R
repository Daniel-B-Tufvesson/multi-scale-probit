library(coda)
library(callr)

source("R/internal.R")

# Fit a multi-scale probit model (MSPM) to the data.  
# 
# Arguments:
# data: An MSPM data structure holding the data.
# burnin: Number of burn-in iterations.
# ndraws: Number of posterior draws to collect.
# thin: Thinning interval for MCMC sampling.
# meanPrior: Prior mean for regression coefficients.
# precPrior: Prior precision for regression coefficients.
# fix.zero: Index of the threshold to fix at zero.
# tune: Tuning parameter for the sampler.
# seed: Random seed for reproducibility.
# beta.initial: Initial values for regression coefficients.
# gamma.initial: Initial values for threshold parameters.
# verbose: Verbosity level for output.
# computeDiagnostics: Whether to compute diagnostics for the fitted model.
# 
# Returns:
# An object of class 'mspm' containing the fitted model.
fit_mspm <- function(
    data, 
    burnin,
    ndraws,
    thin,
    ...,
    meanPrior = NULL,
    precPrior = NULL,
    fix.zero = 1,
    tune = NULL,
    seed = NA,
    beta.initial = NULL,
    gamma.initial = NULL,
    verbose = 0,
    computeDiagnostics = TRUE
) {
    .validate_data(data)

    nlevels = nlevels(data)
    ntargets = ntargets(data)
    npredictors = length(predictorNames(data))

    # Set the seed.
    if (is.na(seed)) {
        seed <- sample.int(.Machine$integer.max, 1)
    }
    set.seed(seed)

    # Set default priors.
    if (is.null(meanPrior)) {
        meanPrior <- rep(0, npredictors)
    }
    if (is.null(precPrior)) {
        precPrior <- .create_prec_prior(npredictors)
    }

    # Set starting values for gamma and beta if not provided
    if (is.null(gamma.initial)) {
        gamma.initial <- .create_inital_gammas(ntargets, nlevels)
    }
    if (is.null(beta.initial)) {
        beta.initial <- rep(0, npredictors)
    }

    # Set tuning parameter for the sampler
    if (is.null(tune)) {
        tune <- 0.05 / nlevels
    } else if (length(tune) != ntargets) {
        tune <- rep(tune, ntargets)
    }

    # Run CPP backend sampler.
    # sim <- cpp_hprobit(
    #     data$Xlist,
    #     data$ylist,
    #     meanPrior,
    #     precPrior,
    #     fix.zero,
    #     nlevels,
    #     gamma.initial,
    #     beta.initial,
    #     tune,
    #     ndraws,
    #     burnin,
    #     thin,
    #     seed,
    #     verbose
    # )

    # Tmp: start as subprocess for more robust development.
    sim <- tryCatch({callr::r(
        function(data, meanPrior, precPrior, fix.zero, nlevels, gamma.initial, beta.initial, tune, ndraws, burnin, thin, seed, verbose) {
            devtools::load_all()
            cpp_hprobit(
                data$Xlist,
                data$ylist,
                meanPrior,
                precPrior,
                fix.zero,
                nlevels,
                gamma.initial,
                beta.initial,
                tune,
                ndraws,
                burnin,
                thin,
                seed,
                verbose
            )
        },
        args = list(data, meanPrior, precPrior, fix.zero, nlevels, gamma.initial, beta.initial, tune, ndraws, burnin, thin, seed, verbose),
        show = FALSE # set to TRUE for debugging
    )}, error = function(e) {
        message("Error in cpp_hprobit: ", e$message)
        if (!is.null(e$stdout)) {
            cat("---- STDOUT ----\n")
            cat(e$stdout, sep = "\n")
        }
        if (!is.null(e$stderr)) {
            cat("---- STDERR ----\n")
            cat(e$stderr, sep = "\n")
        }
        stop(e)
    })

    # Format and store MCMC results
    colnames(sim$storebeta) <- c(data$predictorNames)
    beta <- mcmc(sim$storebeta, start = burnin + 1, thin = thin)
    gammas <- list()
    for (i in 1:ntargets) {
        gammas[[i]] <- mcmc(sim$storegamma[[i]], start = burnin + 1, thin = thin)
    }

    # Compute diagnostics.
    diagnostics <- NULL
    if (computeDiagnostics) {
        diagnostics <- .run_diagnostics(beta, gammas)
    }

    # Return fitted model.
    new_mspm(
        data_spec = data_spec(data),
        beta = beta,
        gammas = gammas,
        meanPrior = meanPrior,
        precPrior = precPrior,
        seed = seed,
        ndraws = ndraws / thin,
        ndrawsNoThin = ndraws,
        thin = thin,
        burnin = burnin,
        diagnostics = diagnostics,
        call = match.call()
    )  
}

.validate_data <- function(data) {
    nlevels = nlevels(data)
    ntargets = ntargets(data)
    npredictors = length(predictorNames(data))

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

.run_diagnostics <- function(beta, gammas) {
    # Compute Geweke diagnostic.
    gewekeBeta <- geweke.diag(beta)
    gewekeGammas <- lapply(gammas, geweke.diag)

    # Compute ESS.
    essBeta <- effectiveSize(beta)
    essGammas <- lapply(gammas, effectiveSize)

    new_mspm_single_chain_diag(
        gewekeBeta = gewekeBeta,
        gewekeGammas = gewekeGammas,
        essBeta = essBeta,
        essGammas = essGammas
    )
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
#' @param computeDiagnostics Whether to compute diagnostics for the fitted model.
#' @return An object of class 'mspm' containing the fitted model.
fit_mspm_pt <- function(
    data, 
    burnin,
    ndraws,
    thin,
    ntemperatures,
    ...,
    meanPrior = NULL,
    precPrior = NULL,
    fix.zero = 1,
    tune = NULL,
    temperatureLadder = NULL,
    target_temp_swap_accept_ratio = 0.3,
    temperature_window_size = 100,
    tempperature_window_growth_factor = 2,
    temperature_ladder_learning_rate = 0.01,
    complete_param_swapping = TRUE,
    seed = NA,
    beta.initial = NULL,
    gamma.initial = NULL,
    verbose = 0,
    computeDiagnostics = TRUE
) {
    .validate_data(data)

    nlevels = nlevels(data)
    ntargets = ntargets(data)
    npredictors = length(predictorNames(data))

    # Set the seed.
    if (is.na(seed)) {
        seed <- sample.int(.Machine$integer.max, 1)
    }
    set.seed(seed)

    # Set default priors.
    if (is.null(meanPrior)) {
        meanPrior <- rep(0, npredictors)
    }
    if (is.null(precPrior)) {
        precPrior <- .create_prec_prior(npredictors)
    }

    # Set starting values for gamma and beta if not provided
    if (is.null(gamma.initial)) {
        gamma.initial <- .create_inital_gammas(ntargets, nlevels)
    }
    if (is.null(beta.initial)) {
        beta.initial <- rep(0, npredictors)
    }

    # Set tuning parameter for the sampler
    if (is.null(tune)) {
        tune <- 0.05 / nlevels
    } else if (length(tune) != ntargets) {
        tune <- rep(tune, ntargets)
    }

    # Set default temperature ladder if not provided.
    if (is.null(temperatureLadder)) {
        temperatureLadder <- 2^(1 - (1:ntemperatures))
    } else if (length(temperatureLadder) != ntemperatures) {
        stop("Length of temperatureLadder must match ntemperatures.")
    }

    # Call the backend parallel tempering sampler.
    # Tmp: start as subprocess for more robust development.
    sim <- tryCatch({callr::r(
        function(
            data, 
            meanPrior, 
            precPrior, 
            fix.zero, 
            nlevels, 
            gamma.initial, 
            beta.initial, 
            tune, 
            ntemperatures,
            temperatureLadder,
            target_temp_swap_accept_ratio,
            temperature_window_size,
            tempperature_window_growth_factor,
            temperature_ladder_learning_rate,
            ndraws,
            burnin,
            thin,
            seed,
            complete_param_swapping,
            verbose
        ) {
            devtools::load_all()
            cpp_hprobit_pt(
                data$Xlist,
                data$ylist,
                meanPrior,
                precPrior,
                fix.zero,
                nlevels,
                gamma.initial,
                beta.initial,
                tune,
                ntemperatures,
                temperatureLadder,
                target_temp_swap_accept_ratio,
                temperature_window_size,
                tempperature_window_growth_factor,
                temperature_ladder_learning_rate,
                ndraws,
                burnin,
                thin,
                seed,
                complete_param_swapping,
                verbose
            )
        },
        args = list(
            data, 
            meanPrior, 
            precPrior, 
            fix.zero, 
            nlevels, 
            gamma.initial, 
            beta.initial, 
            tune, 
            ntemperatures,
            temperatureLadder,
            target_temp_swap_accept_ratio,
            temperature_window_size,
            tempperature_window_growth_factor,
            temperature_ladder_learning_rate,
            ndraws,
            burnin,
            thin,
            seed,
            complete_param_swapping,
            verbose
        ),
        show = TRUE # set to TRUE for debugging
    )}, error = function(e) {
        message("Error in cpp_hprobit_pt: ", e$message)
        if (!is.null(e$stdout)) {
            cat("---- STDOUT ----\n")
            cat(e$stdout, sep = "\n")
        }
        if (!is.null(e$stderr)) {
            cat("---- STDERR ----\n")
            cat(e$stderr, sep = "\n")
        }
        stop(e)
    })

    # Format and store MCMC results
    colnames(sim$storebeta) <- c(data$predictorNames)
    beta <- mcmc(sim$storebeta, start = burnin + 1, thin = thin)
    gammas <- list()
    for (i in 1:ntargets) {
        gammas[[i]] <- mcmc(sim$storegamma[[i]], start = burnin + 1, thin = thin)
    }

    # Compute diagnostics.
    diagnostics <- NULL
    if (computeDiagnostics) {
        diagnostics <- .run_diagnostics(beta, gammas)
    }

    # Return fitted model.
    new_mspm_pt(
        data_spec = data_spec(data),
        beta = beta,
        gammas = gammas,
        meanPrior = meanPrior,
        precPrior = precPrior,
        seed = seed,
        ndraws = ndraws / thin,
        ndrawsNoThin = ndraws,
        thin = thin,
        ntemperatures = ntemperatures,
        burnin = burnin,
        diagnostics = diagnostics,
        initialTemperatureLadder = temperatureLadder,
        adjustedTemperatureLadder = sim$adapted_temps,
        targetSwapAcceptRatio = target_temp_swap_accept_ratio,
        actualSwapAcceptRatio = sim$nswap_accepts / sim$nswap_proposals,
        nacceptedSwaps = sim$nswap_accepts,
        nproposedSwaps = sim$nswap_proposals,
        ladderLearningRate = temperature_ladder_learning_rate,
        initialWindowSize = temperature_window_size,
        windowGrowthFactor = tempperature_window_growth_factor,
        completeSwapping = complete_param_swapping,
        call = match.call()
    )  
}