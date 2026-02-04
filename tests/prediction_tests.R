
source("R/predict.R")

run_all_prediction_tests <- function() {
    .test_predict_latent_draws()
    .test_predict_labels_draws()
    .test_predict_on_new_data()
    .test_predict_reproducibility()
    .test_predict_minimal_data()
    .test_accessors_mspm_latent_prediction()
    .test_accessors_mspm_labeled_prediction()
}

# Test case 1: Generate synthetic data, fit a model, and verify latent variable predictions.
.test_predict_latent_draws <- function() {
    # Generate the data.
    mspm.data <- generate_synthetic_data(
        nobs = 100,
        ncov = 48,
        ngamma = c(1, 3, 3),
        seed = 1234
    )

    # Fit model.
    mspm.fit = fit_mspm(
        data = mspm.data,
        ndraws = 500,
        burnin = 500,
        thin = 1,
        tune = 0.1,
        seed = 1234,
        verbose = 0
    )

    # Predict latent ystar values.
    mspm.latent <- predict_mspm(
        fit = mspm.fit,
        latentOnly = TRUE
    )

    # Check that number of ystar elements match number of targets.
    if (length(mspm.latent$ystars) != mspm.data$ntargets) {
        stop("Test failed: Number of predicted latent variable sets does not match number of targets.")
    }

    # Check that the dimensions of the predicted latent values match the data.
    for (i in 1:mspm.data$ntargets) {

        # Check number of predictions match number of observations
        if (nrow(mspm.latent$ystars[[i]]) != nrow(mspm.data$Xlist[[i]])) {
            stop(paste("Test failed: Number of predicted latent values does not match number of observations for target", i))
        }

        # Check number of predictions match number of draws
        if (ncol(mspm.latent$ystars[[i]]) != mspm.fit$ndraws) {
            stop(paste("Test failed: Number of predicted latent value draws does not match number of posterior draws for target", i))
        }
    }

    cat("Test passed: Predicted latent values have correct dimensions.\n")
}

# Test case 2: Generate synthetic data, fit a model, and verify label predictions.
.test_predict_labels_draws <- function() {
    # Generate the data.
    mspm.data <- generate_synthetic_data(
        nobs = 100,
        ncov = 48,
        ngamma = c(1, 3, 3),
        seed = 1234
    )

    # Fit model.
    mspm.fit = fit_mspm(
        data = mspm.data,
        ndraws = 500,
        burnin = 500,
        thin = 1,
        tune = 0.1,
        seed = 1234,
        verbose = 0
    )

    # Predict labels.
    mspm.labels <- predict_mspm(
        fit = mspm.fit
    )

    # Check that number of ylabels elements match number of targets.
    if (length(mspm.labels$ylabels) != mspm.data$ntargets) {
        stop("Test failed: Number of predicted label sets does not match number of targets.")
    }

    # Check that the predicted labels are factors and have correct lengths.
    for (i in 1:mspm.data$ntargets) {

        # Check that all predicted label strings are in the allowed level names
        if (!all(mspm.labels$ylabels[[i]] %in% mspm.data$levelNames[[i]])) {
            stop(paste("Test failed: Predicted labels for target", i, "are not in allowed level names."))
        }

        # Check number of labels match number of observations
        if (nrow(mspm.labels$ylabels[[i]]) != nrow(mspm.data$Xlist[[i]])) {
            stop(paste("Test failed: Number of predicted labels does not match number of observations for target", i))
        }

        # Check number of labels match number of draws.
        if (ncol(mspm.labels$ylabels[[i]]) != mspm.fit$ndraws) {
            stop(paste("Test failed: Number of predicted label draws does not match number of posterior draws for target", i))
        }

        # Check that predicted label indexes are a matrix.
        if (!is.matrix(mspm.labels$ylabelIndexes[[i]])) {
            stop(paste("Test failed: Predicted label indexes for target", i, "are not a matrix."))
        }

        # Check that all predicted label indexes are valid (within 1:length(levelNames))
        if (!all(mspm.labels$ylabelIndexes[[i]] %in% 1:mspm.data$nlevels[i])) {
            stop(paste("Test failed: Predicted label indexes for target", i, "are not valid indexes."))
        }

        # Check number of label indexes match number of observations
        if (nrow(mspm.labels$ylabelIndexes[[i]]) != nrow(mspm.data$Xlist[[i]])) {
            stop(paste("Test failed: Number of predicted label indexes does not match number of observations for target", i))
        }

        # Check number of label indexes match number of draws.
        if (ncol(mspm.labels$ylabelIndexes[[i]]) != mspm.fit$ndraws) {
            stop(paste("Test failed: Number of predicted label index draws does not match number of posterior draws for target", i))
        }
    }

    cat("Test passed: Predicted labels are factors with correct lengths.\n")
}

# Test case 3: Predict on new data (not used in fitting) and check output dimensions/types.
.test_predict_on_new_data <- function() {
    # Generate training data.
    train.data <- generate_synthetic_data(
        nobs = 100,
        ncov = 48,
        ngamma = c(1, 3, 3),
        seed = 1234
    )

    # Fit model on training data.
    mspm.fit = fit_mspm(
        data = train.data,
        ndraws = 500,
        burnin = 500,
        thin = 1,
        tune = 0.1,
        seed = 1234,
        verbose = 0
    )

    # Generate new test data (same structure, different seed).
    test.data <- generate_synthetic_data(
        nobs = 50,
        ncov = 48,
        ngamma = c(1, 3, 3),
        seed = 5678
    )

    # Predict on new data.
    mspm.pred <- predict_mspm(
        fit = mspm.fit,
        newdata = test.data
    )

    # Check that number of ylabels elements match number of targets.
    if (length(mspm.pred$ylabels) != test.data$ntargets) {
        stop("Test failed: Number of predicted label sets does not match number of targets (new data).")
    }

    # Check that the predicted labels are factors and have correct lengths.
    for (i in 1:test.data$ntargets) {
        # Check that all predicted label strings are in the allowed level names
        if (!all(mspm.pred$ylabels[[i]] %in% test.data$levelNames[[i]])) {
            stop(paste("Test failed: Predicted labels for target", i, "(new data) are not in allowed level names."))
        }
        # Check number of labels match number of observations
        if (nrow(mspm.pred$ylabels[[i]]) != nrow(test.data$Xlist[[i]])) {
            stop(paste("Test failed: Number of predicted labels does not match number of observations for target", i, "(new data)"))
        }
        # Check number of labels match number of draws.
        if (ncol(mspm.pred$ylabels[[i]]) != mspm.fit$ndraws) {
            stop(paste("Test failed: Number of predicted label draws does not match number of posterior draws for target", i, "(new data)"))
        }
    }

    cat("Test passed: Prediction on new data returns correct dimensions and label levels.\n")
}

# Test case 4: Check that predictions are reproducible with the same seed and data.
.test_predict_reproducibility <- function() {
    # Generate the data.
    mspm.data <- generate_synthetic_data(
        nobs = 100,
        ncov = 48,
        ngamma = c(1, 3, 3),
        seed = 4321
    )

    # Fit model.
    mspm.fit = fit_mspm(
        data = mspm.data,
        ndraws = 200,
        burnin = 200,
        thin = 1,
        tune = 0.1,
        seed = 4321,
        verbose = 0
    )

    # Predict twice with the same fit and seed.
    set.seed(999)
    pred1 <- predict_mspm(fit = mspm.fit)
    set.seed(999)
    pred2 <- predict_mspm(fit = mspm.fit)

    # Compare predictions for all targets.
    for (i in seq_along(pred1$ylabels)) {
        if (!identical(pred1$ylabels[[i]], pred2$ylabels[[i]])) {
            stop(paste("Test failed: Predictions are not reproducible for target", i))
        }
    }

    cat("Test passed: Predictions are reproducible with same seed and data.\n")
}

# Test case 5: Predict with minimal data (6 observation, 1 covariate, 3x1 gammas)
.test_predict_minimal_data <- function() {
    cat("Generating minimal data test...\n")
    # Generate minimal data
    min.data <- generate_synthetic_data(
        nobs = 6,
        ncov = 1,
        ngamma = c(1, 1, 1),
        seed = 2026
    )

    cat("Fitting model on minimal data...\n")
    # Fit model
    min.fit = fit_mspm(
        data = min.data,
        ndraws = 10,
        burnin = 10,
        thin = 1,
        tune = 0.1,
        seed = 2026,
        verbose = 0
    )

    cat("Predicting on minimal data...\n")
    # Predict
    min.pred <- predict_mspm(fit = min.fit)

    cat("Checking predictions on minimal data...\n")

    # Check that number of ylabels elements match number of targets
    if (length(min.pred$ylabels) != min.data$ntargets) {
        stop("Test failed: Number of predicted label sets does not match number of targets (minimal data).")
    }

    # Check that the predicted labels have correct shape
    for (i in 1:min.data$ntargets) {
        if (nrow(min.pred$ylabels[[i]]) != 6) {
            stop(paste("Test failed: Number of predicted labels does not match 6 observation for target", i, "(minimal data)"))
        }
        if (ncol(min.pred$ylabels[[i]]) != min.fit$ndraws) {
            stop(paste("Test failed: Number of predicted label draws does not match number of posterior draws for target", i, "(minimal data)"))
        }
    }

    cat("Test passed: Prediction works for minimal data (6 obs, 1 covariate).\n")
}

# Test case 6: Test if the accessor functions work for mspm_latent_prediction object.
.test_accessors_mspm_latent_prediction <- function() {
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
        ndraws = 1000,
        burnin = 500,
        thin = 2,
        tune = 0.1,
        seed = 1234,
        verbose = 0
    )

    # Predict latent ystar values.
    latent <- predict_mspm(
        fit = fit,
        latentOnly = TRUE
    )

    # Test ntargets accessor
    if (ntargets(latent) != mspm.data$ntargets) {
        stop("Test failed: ntargets accessor does not return correct value. Got ", 
             ntargets(latent), " expected ", mspm.data$ntargets, ".")
    }

    # Test nlevels accessor
    if (!all(nlevels(latent) == mspm.data$nlevels)) {
        stop("Test failed: nlevels accessor does not return correct value.")
    }

    # Test predictorNames accessor.
    if (!all.equal(predictorNames(latent), mspm.data$predictorNames)) {
        stop("Test failed: predictorNames accessor returned incorrect value.")
    }

    # Test responseNames accessor.
    if (!all.equal(responseNames(latent), mspm.data$responseNames)) {
        stop("Test failed: responseNames accessor returned incorrect value.")
    }

    # Test levelNames accessor.
    if (!all.equal(levelNames(latent), mspm.data$levelNames)) {
        stop("Test failed: levelNames accessor returned incorrect value.")
    }

    # Test ndraws accessor.
    if (ndraws(latent) != fit$ndraws) {
        stop("Test failed: ndraws accessor does not return correct value.")
    }

    # Test ndraws accessor without thinning.
    if (ndraws(latent, withoutThinning = TRUE) != fit$ndrawsNoThin) {
        stop("Test failed: ndraws accessor with thinning does not return correct value.")
    }

    # Test model accessor.
    if (!identical(model(latent), fit)) {
        stop("Test failed: fit accessor does not return correct model.")
    }

    # Test latent accessor.
    if (!identical(latent(latent), latent$ystars)) {
        stop("Test failed: latent accessor does not return correct latent values.")
    }
    
    cat("Test passed: Accessor functions for mspm_latent_prediction work correctly.\n")
}

# Test case 7: Test if accessor functions work for mspm_labeled_prediction object.
.test_accessors_mspm_labeled_prediction <- function() {
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
        ndraws = 1000,
        burnin = 500,
        thin = 2,
        tune = 0.1,
        seed = 1234,
        verbose = 0
    )

    # Predict labels.
    labeled <- predict_mspm(
        fit = fit
    )

    # Test ntargets accessor
    if (ntargets(labeled) != mspm.data$ntargets) {
        stop("Test failed: ntargets accessor does not return correct value. Got ", 
             ntargets(labeled), " expected ", mspm.data$ntargets, ".")
    }

    # Test nlevels accessor
    if (!all(nlevels(labeled) == mspm.data$nlevels)) {
        stop("Test failed: nlevels accessor does not return correct value.")
    }

    # Test predictorNames accessor.
    if (!all.equal(predictorNames(labeled), mspm.data$predictorNames)) {
        stop("Test failed: predictorNames accessor returned incorrect value.")
    }

    # Test responseNames accessor.
    if (!all.equal(responseNames(labeled), mspm.data$responseNames)) {
        stop("Test failed: responseNames accessor returned incorrect value.")
    }

    # Test levelNames accessor.
    if (!all.equal(levelNames(labeled), mspm.data$levelNames)) {
        stop("Test failed: levelNames accessor returned incorrect value.")
    }

    # Test ndraws accessor.
    if (ndraws(labeled) != fit$ndraws) {
        stop("Test failed: ndraws accessor does not return correct value.")
    }

    # Test ndraws accessor without thinning.
    if (ndraws(labeled, withoutThinning = TRUE) != fit$ndrawsNoThin) {
        stop("Test failed: ndraws accessor with thinning does not return correct value.")
    }

    # Test model accessor.
    if (!identical(model(labeled), fit)) {
        stop("Test failed: fit accessor does not return correct model.")
    }

    # Test predictedLabels accessor.
    if (!identical(predictedLabels(labeled), labeled$ylabels)) {
        stop("Test failed: predictedLabels accessor does not return correct labels.")
    }

    # Test predictedLabelIndexes accessor.
    if (!identical(predictedLabelIndexes(labeled), labeled$ylabelIndexes)) {
        stop("Test failed: predictedLabelIndexes accessor does not return correct label indexes.")
    }

    cat("Test passed: Accessor functions for mspm_labeled_prediction work correctly.\n")
}