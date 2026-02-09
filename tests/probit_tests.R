source("R/fit.R")
source("R/util.R")
source("R/data.R")

run_all_probit_tests <-function() {
    .test_fit()
    .test_fit_compare()
    .test_accessors()
}


# Test case 1: Fit a probit model and verify that it can estimate the original parameters.
.test_fit <- function() {
    # Initial parameter tuning code
    burnin <- 2000
    ndraws <- 2000
    thin <- 10
    verbose <- 0
    beta = c(0.5, -1.0, 1.5)
    gammas = list(
        c(0), 
        c(-1, 1, 2), 
        c(-0.5, 0.5, 2.5)
    )

    # Generate the data.
    mspm.data <- generate_synthetic_data(
        nobs = 100,
        ncov = length(beta),
        ngamma = c(1, 3, 3),
        seed = 1234,
        beta = beta,
        gammas = gammas
    )

    # Fit using all the data.
    tst_fit = mspm.fit <- fit_mspm(
        data = mspm.data,
        ndraws = ndraws,
        burnin = burnin,
        thin = thin,
        tune = 0.1,
        seed = 1234,
        verbose = verbose
    )

    failed <- FALSE

    # Check beta.
    beta_estimates <- colMeans(as.matrix(tst_fit$beta))
    if (max(abs(beta_estimates - beta)) > 0.2) {
        stop("Test failed: Estimated beta coefficients deviate significantly from true values.")
    }

    # Check gammas.
    for (i in 1:length(gammas)) {
        gamma_estimates <- colMeans(as.matrix(tst_fit$gammas[[i]]))
        if (max(abs(gamma_estimates - gammas[[i]])) > 0.5) {
            stop(paste("Test failed: Estimated gammas for dataset", i, 
                       "deviate significantly from true values."))
        }
    }

    # Check diagnostics.
    if (is.null(diagnostics(tst_fit))) {
        stop("Test failed: Diagnostics are NULL.")
    }

    # Check essBeta.
    if(is.null(essBeta(tst_fit))) {
        stop("Test failed: essBeta is NULL.")
    }
    if(length(essBeta(tst_fit)) != length(beta)) {
        stop("Test failed: essBeta length does not match number of beta coefficients.")
    }

    # Check essGammas.
    if(is.null(essGammas(tst_fit))) {
        stop("Test failed: essGammas is NULL.")
    }
    for (i in 1:length(gammas)) {
        if(is.null(essGammas(tst_fit)[[i]])) {
            stop(paste("Test failed: essGammas for dataset", i, "is NULL."))
        }
        if(length(essGammas(tst_fit)[[i]]) != length(gammas[[i]])) {
            stop(paste("Test failed: essGammas length for dataset", i, 
                       "does not match number of gamma thresholds."))
        }
    }

    # Check gewekeBeta.
    if(is.null(gewekeBeta(tst_fit))) {
        stop("Test failed: gewekeBeta is NULL.")
    }
    if (length(gewekeBeta(tst_fit)$z) != length(beta)) {
        stop("Test failed: gewekeBeta length does not match number of beta coefficients.")
    }
    
    # Check gewekeGammas.
    if(is.null(gewekeGammas(tst_fit))) {
        stop("Test failed: gewekeGammas is NULL.")
    }
    for (i in 1:length(gammas)) {
        if(is.null(gewekeGammas(tst_fit)[[i]])) {
            stop(paste("Test failed: gewekeGammas for dataset", i, "is NULL."))
        }
        if(length(gewekeGammas(tst_fit)[[i]]$z) != length(gammas[[i]])) {
            stop(paste("Test failed: gewekeGammas length for dataset", i, 
                       "does not match number of gamma thresholds."))
        }
    }

    cat("Test passed: Fitted model estimates parameters accurately and returns valid diagnostics.\n")
}


# Test case 2: Fit two probit models on the same data and verify that the 
# results are identical.
.test_fit_compare <- function() {
    # Initial parameter tuning code
    burnin <- 2000
    ndraws <- 2000
    thin <- 10
    verbose <- 0

    # Generate the data.
    mspm.data <- generate_synthetic_data(
        nobs = 100,
        ncov = 5,
        ngamma = c(2, 3),
        seed = 42
    )

    # Fit using all the data.
    fit1 = fit_mspm(
        data = mspm.data,
        ndraws = ndraws,
        burnin = burnin,
        thin = thin,
        tune = 0.1,
        seed = 1234,
        verbose = verbose
    )

    fit2 = fit_mspm(
        data = mspm.data,
        ndraws = ndraws,
        burnin = burnin,
        thin = thin,
        tune = 0.1,
        seed = 1234,
        verbose = verbose
    )

    # Check that the fits are identical.
    if (!all.equal(fit1$beta, fit2$beta)) {
        stop("Test failed: Fitted beta coefficients are not identical.")
    }

    for (i in 1:length(fit1$gammas)) {
        if (!all.equal(fit1$gammas[[i]], fit2$gammas[[i]])) {
            stop(paste("Test failed: Fitted gammas for dataset", i, "are not identical."))
        }
    }

    cat("Test passed: Fitted models are identical.\n")
}

# Test case 3: Test if the accessor functions work.
.test_accessors <- function() {
    # Generate the data.
    mspm.data <- generate_synthetic_data(
        nobs = 100,
        ncov = 5,
        ngamma = c(2, 3),
        seed = 42
    )

    # Fit using all the data.
    fit = fit_mspm(
        data = mspm.data,
        ndraws = 500,
        burnin = 500,
        thin = 1,
        tune = 0.1,
        seed = 1234,
        verbose = 0
    )

    # Test ntargets accessor.
    if (ntargets(fit) != mspm.data$ntargets) {
        stop("Test failed: ntargets accessor returned incorrect value.")
    }

    # Test nlevels accessor.
    if (nlevels(fit) != mspm.data$nlevels) {
        stop("Test failed: nlevels accessor returned incorrect value.")
    }

    # Test predictorNames accessor.
    if (!all.equal(predictorNames(fit), mspm.data$predictorNames)) {
        stop("Test failed: predictorNames accessor returned incorrect value.")
    }

    # Test responseNames accessor.
    if (!all.equal(responseNames(fit), mspm.data$responseNames)) {
        stop("Test failed: responseNames accessor returned incorrect value.")
    }

    # Test levelNames accessor.
    if (!all.equal(levelNames(fit), mspm.data$levelNames)) {
        stop("Test failed: levelNames accessor returned incorrect value.")
    }

    # Test beta accessor.
    if (!all.equal(beta(fit), fit$beta)) {
        stop("Test failed: beta accessor returned incorrect value.")
    }

    # Test gammas accessor.
    if (length(gammas(fit)) != length(fit$gammas)) {
        stop("Test failed: gammas accessor returned incorrect value.")
    }

    # Test meanPrior accessor.
    if (!all.equal(meanPrior(fit), fit$meanPrior)) {
        stop("Test failed: meanPrior accessor returned incorrect value.")
    }

    # Test precPrior accessor.
    if (!all.equal(precPrior(fit), fit$precPrior)) {
        stop("Test failed: precPrior accessor returned incorrect value.")
    }

    # Test ndraws accessor.
    if (ndraws(fit) != fit$ndraws) {
        stop("Test failed: ndraws accessor returned incorrect value.")
    }

    # Test ndraws with thinning accessor.
    if (ndraws(fit, withoutThinning = TRUE) != fit$ndrawsNoThin) {
        stop("Test failed: ndraws without thinning accessor returned incorrect value.")
    }

    # Test burnin accessor.
    if (burnin(fit) != fit$burnin) {
        stop("Test failed: burnin accessor returned incorrect value.")
    }

    # Test thin accessor.
    if (thin(fit) != fit$thin) {
        stop("Test failed: thin accessor returned incorrect value.")
    }

    cat("Test passed: Accessor functions work correctly.\n")
}