source("R/persistence.R")
source("R/data.R")
source("R/fit.R")
source("R/internal.R")


run_all_persistence_tests <- function() {
    .test_save_and_load()
    .test_predict()
}

#' Test case 1: Save a fitted model to disk and load it back, verifying key components are 
#' identical.
.test_save_and_load <- function() {
    # Generate synthetic data and fit model.
    data <- generate_synthetic_data(
        nobs = 100,
        ncov = 5,
        ngamma = c(2, 3),
        seed = 42
    )

    # Fit model.
    fit = fit_mspm(
        data = mspm.data,
        ndraws = 100,
        burnin = 100,
        thin = 2,
        tune = 0.1,
        seed = 1234,
        verbose = 0
    )

    # Save model to temporary file.
    temp_file <- tempfile(fileext = ".rds")
    save_model(fit, temp_file)

    # Load model from file.
    loaded_fit <- load_model(temp_file)

    # Check beta.
    if (!identical(beta(fit), beta(loaded_fit))) {
        stop("Test failed: beta are not identical after loading.")
    }
    
    # Check gammas.
    if (!identical(gammas(fit), gammas(loaded_fit))) {
        stop("Test failed: gammas are not identical after loading.")
    }

    # Check ntargets.
    if (ntargets(fit) != ntargets(loaded_fit)) {
        stop("Test failed: ntargets are not identical after loading.")
    }

    # Check nlevels.
    if (!identical(nlevels(fit), nlevels(loaded_fit))) {
        stop("Test failed: nlevels are not identical after loading.")
    }

    # Check meanPrior.
    if (!identical(meanPrior(fit), meanPrior(loaded_fit))) {
        stop("Test failed: meanPrior are not identical after loading.")
    }

    # Check precPrior.
    if (!identical(precPrior(fit), precPrior(loaded_fit))) {
        stop("Test failed: precPrior are not identical after loading.")
    }

    # Check ndraws.
    if (ndraws(fit) != ndraws(loaded_fit)) {
        stop("Test failed: ndraws are not identical after loading.")
    }

    # Check ndraws without thinning.
    if (ndraws(fit, withoutThinning = FALSE) != ndraws(loaded_fit, withoutThinning = FALSE)) {
        stop("Test failed: ndraws without thinning are not identical after loading.")
    }

    # Check burnin.
    if (burnin(fit) != burnin(loaded_fit)) {
        stop("Test failed: burnin are not identical after loading.")
    }

    # Check thin.
    if (thin(fit) != thin(loaded_fit)) {
        stop("Test failed: thin are not identical after loading.")
    }

    # Check training data. Should be null.
    if (!is.null(loaded_fit$data)) {
        stop("Test failed: training data should be null after loading.")
    }

    cat("Test passed: Model saved and loaded successfully.\n")
}

#' Test case 2: Predict using the loaded model and verify predictions are consistent with the 
#' original model.
.test_predict <- function() {
    # Generate synthetic data and fit model.
    data <- generate_synthetic_data(
        nobs = 100,
        ncov = 5,
        ngamma = c(2, 3),
        seed = 42
    )

    # Fit model.
    fit = fit_mspm(
        data = data,
        ndraws = 100,
        burnin = 100,
        thin = 2,
        tune = 0.1,
        seed = 1234,
        verbose = 0
    )

    # Save model to temporary file.
    temp_file <- tempfile(fileext = ".rds")
    save_model(fit, temp_file)

    # Load model from file.
    loaded_fit <- load_model(temp_file)

    #data <- fit$data
    predict1 <- predict_mspm(fit, newdata = data, type = "draws")
    predict2 <- predict_mspm(loaded_fit, newdata = data, type = "draws")

    # Check that latent predictions are identical.
    if (!identical(latent(predict1), latent(predict2))) {
        stop("Test failed: latent predictions are not identical after loading.")
    }

    # Check predictedLabels.
    if (!identical(predictedLabels(predict1), predictedLabels(predict2))) {
        cat("diff: ", all.equal(predictedLabels(predict1), predictedLabels(predict2)), "\n")
        cat("len: ", length(predictedLabels(predict1)), ", ", length(predictedLabels(predict2)), "\n")
        stop("Test failed: predicted labels are not identical after loading.")
    }

    cat("Test passed: Predictions from loaded model are consistent with original model.\n")
}