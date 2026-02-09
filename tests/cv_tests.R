source("R/data.R")
source("R/cv.R")


run_all_cv_tests <- function() {
    .test_cv()
    .test_replication()

    # Fix: parallelization fails stochastically, likely due to random number generation 
    # issues. Need to investigate and fix before enabling these tests.
    #.test_parallelization()
    #.test_parallelization_diff_nworkers()
}

# Test case 1: Test that cross-validation runs without errors and returns results in the expected format.
.test_cv <- function() {
    # Genera te synthetic data for testing.
    data <- generate_synthetic_data(
        nobs = 30,
        ncov = 6,
        ngamma = c(1, 3, 3),
        seed = 1234
    )

    # Run cross-validation with a small number of splits and draws for testing.
    cv_res <- cross_validate(
        data = data,
        nsplits = 3,
        prop = 0.7,
        ndraws = 10,
        burnin = 10,
        thin = 1,
        seed = 42,
        nworkers = 1,
        meansOnly = FALSE
    )

    if (class(cv_res) != "mspm_cv_result") {
        stop("Expected cv_res to be an mspm_cv_result.")
    }

    # Check nsplits.
    if (nsplits(cv_res) != 3) {
        stop("Expected nsplits to be 3.")
    }

    # Check metrics.
    expectedMetrics <- c("f1", "kendall")
    if (!identical(expectedMetrics, evalMetrics(cv_res))) {
        stop("Expected evalMetrics to match the specified metrics. ", 
             "Expected: ", paste(expectedMetrics, collapse = ", "), 
             " Got: ", paste(evalMetrics(cv_res), collapse = ", ")
        )
    }

    # Check data_spec.
    if (!identical(data_spec(cv_res), data_spec(data))) {
        stop("Expected data_spec in cv_res to match the original data_spec.")
    }

    # Check allEvals.
    allEvals <- cvAllEvaluations(cv_res)
    if (length(allEvals) != 3) {
        stop("Expected allEvals to have length equal to nsplits.")
    }

    # Check means.
    means <- cvMeans(cv_res)
    if (length(means) != 3) {
        stop("Expected means to have length equal to ntargets.")
    }

    # Check metrics in means.
    for (i in 1:3) {
        mean_i <- means[[i]]

        # Check metrics.
        if (!all(expectedMetrics %in% colnames(mean_i))) {
            stop(paste("Expected means[[", i, "]] to contain all specified metrics."))
        }

        # Check number of columns.
        expectedNCols <- if (length(expectedMetrics) > 1) length(expectedMetrics) + 1 else length(expectedMetrics)
        if (ncol(mean_i) != expectedNCols) {
            stop(paste("Expected means[[", i, "]] to have number of columns equal to number of metrics + 1 for harmonic mean."))
        }

        # Check column names.
        expectedColNames <- if (length(expectedMetrics) > 1) c(expectedMetrics, "HarmonicMean") else expectedMetrics
        if (!identical(colnames(mean_i), expectedColNames)) {
            stop(paste("Expected means[[", i, "]] to have column names matching the specified metrics and harmonic mean."))
        }

        # Check number of rows. One for each split.
        if (nrow(mean_i) != 3) {
            stop(paste("Expected means[[", i, "]] to have number of rows equal to nsplits."))
        }
    }

    cat("Test passed: cross-validation runs without errors and returns results in the expected format.\n")
}

# Test case 2: Test that cross-validation results are reproducible with the same seed.
.test_replication <- function() {
    # Generate synthetic data for testing.
    data <- generate_synthetic_data(
        nobs = 30,
        ncov = 3,
        ngamma = c(1, 3, 3),
        seed = 1234
    )

    # Run cross-validation with a fixed seed and check that results are the same across runs.
    cv_res1 <- cross_validate(
        data = data,
        nsplits = 3,
        prop = 0.7,
        ndraws = 10,
        burnin = 10,
        thin = 1,
        seed = 42,
        nworkers = 1,
        meansOnly = FALSE
    )

    cv_res2 <- cross_validate(
        data = data,
        nsplits = 3,
        prop = 0.7,
        ndraws = 10,
        burnin = 10,
        thin = 1,
        seed = 42,
        nworkers = 1,
        meansOnly = FALSE
    )

    # Check that allEvals are identical.
    if (!identical(cvAllEvaluations(cv_res1), cvAllEvaluations(cv_res2))) {
        stop("Expected allEvals to be identical across runs with the same seed.")
    }

    # Check that means are identical.
    if (!identical(cvMeans(cv_res1), cvMeans(cv_res2))) {
        stop("Expected means to be identical across runs with the same seed.")
    }

    cat("Test passed: cross-validation results are reproducible with the same seed.\n")
}

# Test case 3: Test that cross-validation replicates for parallelization.
.test_parallelization <- function() {
    # Generate synthetic data for testing.
    data <- generate_synthetic_data(
        nobs = 30,
        ncov = 3,
        ngamma = c(1, 3, 3),
        seed = 1234
    )

    cv_res1 <- cross_validate(
        data = data,
        nsplits = 3,
        prop = 0.7,
        ndraws = 10,
        burnin = 10,
        thin = 1,
        seed = 42,
        nworkers = 4,
        meansOnly = FALSE
    )

    cv_res2 <- cross_validate(
        data = data,
        nsplits = 3,
        prop = 0.7,
        ndraws = 10,
        burnin = 10,
        thin = 1,
        seed = 42,
        nworkers = 4,
        meansOnly = FALSE
    )

    # Check that allEvals are identical.
    if (!identical(cvAllEvaluations(cv_res1), cvAllEvaluations(cv_res2))) {
        stop("Expected allEvals to be identical across runs with the same seed and parallelization.")
    }

    # Check that means are identical.
    if (!identical(cvMeans(cv_res1), cvMeans(cv_res2))) {
        stop("Means do not replicate across runs with parallelization.")
    }

    cat("Test passed: cross-validation results replicate with parallelization.\n")
}

# Test case 4: 
.test_parallelization_diff_nworkers <- function() {
    # Generate synthetic data for testing.
    data <- generate_synthetic_data(
        nobs = 30,
        ncov = 3,
        ngamma = c(1, 3, 3),
        seed = 1234
    )

    cv_res1 <- cross_validate(
        data = data,
        nsplits = 3,
        prop = 0.7,
        ndraws = 10,
        burnin = 10,
        thin = 1,
        seed = 42,
        nworkers = 2,
        meansOnly = FALSE
    )

    cv_res2 <- cross_validate(
        data = data,
        nsplits = 3,
        prop = 0.7,
        ndraws = 10,
        burnin = 10,
        thin = 1,
        seed = 42,
        nworkers = 4,
        meansOnly = FALSE
    )

    # Check that allEvals are identical.
    if (!identical(cvAllEvaluations(cv_res1), cvAllEvaluations(cv_res2))) {
        stop("Expected allEvals to be identical across runs with the same seed and parallelization",
             " with different number of workers."
        )
    }

    # Check that means are identical.
    if (!identical(cvMeans(cv_res1), cvMeans(cv_res2))) {
        stop("Means do not replicate across runs with parallelization with different number of workers.")
    }

    cat("Test passed: cross-validation results replicate with parallelization with different number of workers.\n")
}