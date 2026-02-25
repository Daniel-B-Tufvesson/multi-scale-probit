
source("R/predict.R")

run_all_pt_tests <- function() {
    .test_pt_predict()
    .test_predict_reproducibility()
}

#' Test case 1: Fit a modeel with parallel tempering and test the predict function. 
.test_pt_predict <- function() {
    # Generate the data.
    data <- generate_synthetic_data(
        nobs = 60,
        ncov = 6,
        ngamma = c(1, 3, 3),
        seed = 1234
    )

    # Fit using all the data.
    fit <- fit_mspm_pt(
        data = data,
        ndraws = 100,
        burnin = 100,
        thin = 2,
        tune = 0.1,
        seed = 1234,
        ntemperatures = 10,
        verbose = 0
    )

    # Validate accessors.
    # Check beta.
    beta <- beta(fit);
    if (is.null(beta)) {
        stop("Beta is NULL")
    }
    # Check beta size.
    expected_beta_size <- 6
    if (ncol(beta) != expected_beta_size) {
        stop(paste("Wrong column size for beta. Expected ",
            expected_beta_size, " but got ", length(beta)))
    }
    # Check draws in beta.
    expected_draws <- 50
    if (nrow(beta) != expected_draws) {
        stop(paste("Wrong row size for beta. Expected ", 
            expected_draws, " but got ", nrow(beta)))
    }
    
    # Check gamma.
    gammas <- gammas(fit);
    if (is.null(gammas)) {
        stop("Gammas is NULL")
    }
    # Check gamma groups.
    if (length(gammas) != 3) {
        stop(paste("Wrong number of gamma groups. Expected 3 but got ", length(gammas)))
    }
    # Check gamma group sizes.
    expected_gamma_group_sizes <- c(1, 3, 3)
    for (i in 1:length(gammas)) {
        if (ncol(gammas[[i]]) != expected_gamma_group_sizes[i]) {
            stop(paste("Wrong column size for gamma group ", i, ". Expected ", 
                expected_gamma_group_sizes[i], " but got ", ncol(gammas[[i]])))
        }
        # Check draws in gamma group.
        if (nrow(gammas[[i]]) != expected_draws) {
            stop(paste("Wrong row size for gamma group ", i, ". Expected ", 
                expected_draws, " but got ", nrow(gammas[[i]])))
        }
    }

    # Check data_spec.
    data_spec <- data_spec(fit)
    if (is.null(data_spec)) {
        stop("Data spec is NULL")
    }
    if (!identical(data_spec, data_spec(data))) {
        stop("Data spec in fit does not match original data spec")
    }

    # Chek ndraws.
    if (ndraws(fit) != expected_draws) {
        stop(paste("Wrong number of draws. Expected ", 
            expected_draws, " but got ", ndraws(fit)))
    }

    # Check ndrawsNoThin.
    expected_draws_no_thin <- 100
    if (ndrawsNoThin(fit) != expected_draws_no_thin) {
        stop(paste("Wrong number of draws without thinning. Expected ", 
            expected_draws_no_thin, " but got ", ndrawsNoThin(fit)))
    }

    # Check thin.
    expected_thin <- 2
    if (thin(fit) != expected_thin) {
        stop(paste("Wrong thin. Expected ", expected_thin, " but got ", thin(fit)))
    }

    # Check burnin.
    expected_burnin <- 100
    if (burnin(fit) != expected_burnin) {
        stop(paste("Wrong burnin. Expected ", expected_burnin, " but got ", burnin(fit)))
    }

    # Check diagnostics.
    diagnostics <- diagnostics(fit)
    if (is.null(diagnostics)) {
        stop("Diagnostics is NULL")
    }

    # Check ntemperatures.
    expected_ntemperatures <- 10
    if (ntemperatures(fit) != expected_ntemperatures) {
        stop(paste("Wrong number of temperatures. Expected ", 
            expected_ntemperatures, " but got ", ntemperatures(fit)))
    }

    # Make prediction.
    predictions <- predict_mspm(
        fit = fit,
        newdata = data
    )
    labels <- predictedLabels(predictions)
    labelIndexes <- predictedLabelIndexes(predictions)
    levelNames <- levelNames(data)
    nlevels <- nlevels(data)

    # Check that number of ylabels elements match number of targets.
    if (length(labels) != ntargets(data)) {
        stop("Test failed: Number of predicted label sets does not match number of targets.")
    }

    # Check that the predicted labels are factors and have correct lengths.
    for (i in 1:ntargets(data)) {

        # Check that all predicted label strings are in the allowed level names
        if (!all(labels[[i]] %in% levelNames[[i]])) {
            stop(paste("Test failed: Predicted labels for target", i, "are not in allowed level names."))
        }

        # Check number of labels match number of observations
        if (nrow(labels[[i]]) != nrow(data$Xlist[[i]])) {
            stop(paste("Test failed: Number of predicted labels does not match number of observations for target", i))
        }

        # Check number of labels match number of draws.
        if (ncol(labels[[i]]) != ndraws(fit)) {
            stop(paste("Test failed: Number of predicted label draws does not match number of posterior draws for target", i))
        }

        # Check that predicted label indexes are a matrix.
        if (!is.matrix(labelIndexes[[i]])) {
            stop(paste("Test failed: Predicted label indexes for target", i, "are not a matrix."))
        }

        # Check that all predicted label indexes are valid (within 1:length(levelNames))
        if (!all(labelIndexes[[i]] %in% 1:nlevels[i])) {
            stop(paste("Test failed: Predicted label indexes for target", i, "are not valid indexes."))
        }

        # Check number of label indexes match number of observations
        if (nrow(labelIndexes[[i]]) != nrow(data$Xlist[[i]])) {
            stop(paste("Test failed: Number of predicted label indexes does not match number of observations for target", i))
        }

        # Check number of label indexes match number of draws.
        if (ncol(labelIndexes[[i]]) != ndraws(fit)) {
            stop(paste("Test failed: Number of predicted label index draws does not match number of posterior draws for target", i))
        }
    }

    cat("Test passed: predict_mspm returns predicted labels and label indexes in the expected format with correct dimensions and valid values.\n")
}


#' Test case 2: Test that predictions are reproducible with the same seed.
.test_predict_reproducibility <- function() {
    # Generate the data.
    data <- generate_synthetic_data(
        nobs = 60,
        ncov = 6,
        ngamma = c(1, 3, 3),
        seed = 1234
    )

    # Fit using all the data.
    fit1 <- fit_mspm_pt(
        data = data,
        ndraws = 100,
        burnin = 100,
        thin = 2,
        tune = 0.1,
        seed = 1234,
        ntemperatures = 10,
        verbose = 0
    )
    fit2 <- fit_mspm_pt(
        data = data,
        ndraws = 100,
        burnin = 100,
        thin = 2,
        tune = 0.1,
        seed = 1234,
        ntemperatures = 10,
        verbose = 0
    )

    # Make predictions.
    pred1 <- predict_mspm(
        fit = fit1,
        newdata = data
    )
    pred2 <- predict_mspm(
        fit = fit2,
        newdata = data
    )

    pred1_gl <<- pred1
    pred2_gl <<- pred2

    # Check that predicted labels are the same for both fits.
    if (!identical(predictedLabels(pred1), predictedLabels(pred2))) {
        stop("Test failed: Predictions are not reproducible with same seed and data.")
    }

    # Compare label indexes.
    if (!identical(predictedLabelIndexes(pred1), predictedLabelIndexes(pred2))) {
        stop("Test failed: Predicted label indexes are not reproducible with same seed and data.")
    }

    cat("Test passed: Predictions are reproducible with same seed and data.\n")
}