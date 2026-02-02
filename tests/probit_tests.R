source("R/fit.R")
source("R/util.R")
source("R/data.R")

run_all_probit_tests <-function() {
    .test_generate_data()
    .test_fit()
    .test_fit_compare()
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
        failed <- TRUE
        warning("Test failed: Estimated beta coefficients deviate significantly from true values.")
    }

    # Check gammas.
    for (i in 1:length(gammas)) {
        gamma_estimates <- colMeans(as.matrix(tst_fit$gammas[[i]]))
        if (max(abs(gamma_estimates - gammas[[i]])) > 0.5) {
            failed <- TRUE
            warning(paste("Test failed: Estimated gammas for dataset", i, "deviate significantly from true values."))
        }
    }

    if (!failed) {
        cat("Test passed: Fitted model estimates parameters accurately.\n")
    }
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
