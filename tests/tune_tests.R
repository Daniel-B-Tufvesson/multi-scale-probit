# Some tests for automatic tuning of the samplers.

# Test case 1: Test that the tuning function runs for the MH-gibbs sampler without errors and 
# returns results in the expected format.
run_all_tune_tests <- function() {
    .test_tune_proposal_variance_gibbs()
    .test_tune_proposal_variance_pt()
}

.test_tune_proposal_variance_gibbs <- function() {
    # Generate the data.
    data <- generate_synthetic_data(
        nobs = 60,
        ncov = 5,
        ngamma = c(1, 3, 3),
        seed = 1234
    )

    # Find the optimal tuning parameters for the sampler.
    tune_results <- tune_mspm(
        data = data,
        iterations = 100,
        target_epsilon = 0.01,
        stop_early = TRUE
    )

    # Validate accessors.
    proposal_variance <- proposal_variance(tune_results)
    if (is.null(proposal_variance) || !is.numeric(proposal_variance)) {
        stop("Test failed: proposal_variance accessor did not return a numeric value.")
    }

    # Fit using all the data and the tuned parameters.
    fit <- fit_mspm(
        data = data,
        ndraws = 10,
        burnin = 10,
        thin = 1,
        proposal_variance = proposal_variance,
        seed = 1234
    )

    # Validate accessor.
    proposal_variance_fit <- proposal_variance(fit)
    if (is.null(proposal_variance_fit) || !is.numeric(proposal_variance_fit)) {
        stop("Test failed: proposal_variance accessor did not return a numeric value for the fit object.")
    }

    if (!identical(proposal_variance, proposal_variance_fit)) {
        stop("Test failed: proposal_variance from tuning results and fit object do not match.")
    }

    cat("Test passed: Tuning of proposal variance for Gibbs sampler completed successfully.\n")
}

# Test case 2: Test that the tuning function runs for parallel tempering without errors and returns 
# results in the expected format.
.test_tune_proposal_variance_pt <- function() {
    # Generate the data.
    data <- generate_synthetic_data(
        nobs = 60,
        ncov = 5,
        ngamma = c(1, 3, 3),
        seed = 1234
    )

    # Find the optimal tuning parameters for the sampler.
    tune_results <- tune_mspm_pt(
        data = data,
        iterations = 100,
        target_epsilon = 0.01,
        stop_early = TRUE
    )

    # Validate accessors.
    proposal_variance <- proposal_variance(tune_results)
    if (is.null(proposal_variance) || !is.list(proposal_variance) || !all(sapply(proposal_variance, is.numeric))) {
        stop("Test failed: proposal_variance accessor did not return a list of numeric vectors.")
    }

    # Fit using all the data and the tuned parameters.
    fit <- fit_mspm_pt(
        data = data,
        ndraws = 10,
        burnin = 10,
        thin = 1,
        proposal_variance = proposal_variance,
        seed = 1234
    )

    # Validate accessor.
    proposal_variance_fit <- proposal_variance(fit)
    if (is.null(proposal_variance_fit) || !is.list(proposal_variance_fit) || !all(sapply(proposal_variance_fit, is.numeric))) {
        stop("Test failed: proposal_variance accessor did not return a list of numeric vectors.")
    }

    if (!identical(proposal_variance, proposal_variance_fit)) {
        stop("Test failed: proposal_variance from tuning results and fit object do not match.")
    }

    cat("Test passed: Tuning of proposal variance for PT sampler completed successfully.\n")
}