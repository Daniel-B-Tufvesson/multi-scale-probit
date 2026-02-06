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
    verbose = 0
) {

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
        if (npredictors > 1) {
            precPrior <- diag(rep(0.1, npredictors))
        }
        else {
            # Special case for single predictor. diag() returns a 0x0 matrix when 
            # npredictors=1.
            precPrior <- matrix(0.1, nrow = 1, ncol = 1)
        }
    }

    # Set starting values for gamma and beta if not provided
    if (is.null(gamma.initial)) {
        gamma.initial <- list()
        for (i in 1:ntargets) {
            gamma.initial[[i]] <- rep(0, nlevels[i]+1)
            gamma.initial[[i]][1] <- -300
            gamma.initial[[i]][2:nlevels[i]] <- 0:(nlevels[i]-2)
            gamma.initial[[i]][nlevels[i] + 1] <- 300
        }
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
        call = match.call()
    )  
}
