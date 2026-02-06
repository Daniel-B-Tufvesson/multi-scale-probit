
source("R/predict.R")

run_all_prediction_tests <- function() {
    .test_predict_latent_draws()
    .test_predict_labels_draws()
    .test_predict_on_new_data()
    .test_predict_reproducibility()
    .test_predict_minimal_data()
    .test_accessors_mspm_labeled_prediction()
}

# Test case 1: Generate synthetic data, fit a model, and verify latent variable predictions.
.test_predict_latent_draws <- function() {
    # Generate the data.
    data <- generate_synthetic_data(
        nobs = 100,
        ncov = 48,
        ngamma = c(1, 3, 3),
        seed = 1234
    )

    # Fit model.
    fit = fit_mspm(
        data = data,
        ndraws = 500,
        burnin = 500,
        thin = 1,
        tune = 0.1,
        seed = 1234,
        verbose = 0
    )

    # Predict latent ystar values.
    latent <- predict_mspm(
        fit = fit,
        newdata = data,
        latentOnly = TRUE
    )

    # Check that number of ystar elements match number of targets.
    if (length(latent) != ntargets(data)) {
        stop("Test failed: Number of predicted latent variable sets does not match number of targets.")
    }

    # Check that the dimensions of the predicted latent values match the data.
    for (i in 1:ntargets(data)) {

        # Check number of predictions match number of observations
        if (nrow(latent[[i]]) != nrow(data$Xlist[[i]])) {
            stop(paste("Test failed: Number of predicted latent values does not match number of observations for target", i))
        }

        # Check number of predictions match number of draws
        if (ncol(latent[[i]]) != ndraws(fit)) {
            stop(paste("Test failed: Number of predicted latent value draws does not match number of posterior draws for target", i))
        }
    }

    cat("Test passed: Predicted latent values have correct dimensions.\n")
}

# Test case 2: Generate synthetic data, fit a model, and verify label predictions.
.test_predict_labels_draws <- function() {
    # Generate the data.
    data <- generate_synthetic_data(
        nobs = 100,
        ncov = 48,
        ngamma = c(1, 3, 3),
        seed = 1234
    )

    # Fit model.
    fit = fit_mspm(
        data = data,
        ndraws = 500,
        burnin = 500,
        thin = 1,
        tune = 0.1,
        seed = 1234,
        verbose = 0
    )

    # Predict labels.
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
    fit = fit_mspm(
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
    ntargets <- ntargets(test.data)
    levelNames <- levelNames(test.data)

    # Predict on new data.
    pred <- predict_mspm(
        fit = fit,
        newdata = test.data
    )
    labels <- predictedLabels(pred)

    # Check that number of ylabels elements match number of targets.
    if (length(labels) != ntargets) {
        stop("Test failed: Number of predicted label sets does not match number of targets (new data).")
    }

    # Check that the predicted labels are factors and have correct lengths.
    for (i in 1:ntargets) {
        # Check that all predicted label strings are in the allowed level names
        if (!all(labels[[i]] %in% levelNames[[i]])) {
            stop(paste("Test failed: Predicted labels for target", i, "(new data) are not in allowed level names."))
        }
        # Check number of labels match number of observations
        if (nrow(labels[[i]]) != nrow(test.data$Xlist[[i]])) {
            stop(paste("Test failed: Number of predicted labels does not match number of observations for target", i, "(new data)"))
        }
        # Check number of labels match number of draws.
        if (ncol(labels[[i]]) != ndraws(fit)) {
            stop(paste("Test failed: Number of predicted label draws does not match number of posterior draws for target", i, "(new data)"))
        }
    }

    cat("Test passed: Prediction on new data returns correct dimensions and label levels.\n")
}

# Test case 4: Check that predictions are reproducible with the same seed and data.
.test_predict_reproducibility <- function() {
    # Generate the data.
    data <- generate_synthetic_data(
        nobs = 100,
        ncov = 48,
        ngamma = c(1, 3, 3),
        seed = 4321
    )

    # Fit model.
    fit = fit_mspm(
        data = data,
        ndraws = 200,
        burnin = 200,
        thin = 1,
        tune = 0.1,
        seed = 4321,
        verbose = 0
    )

    # Predict twice with the same fit and seed.
    set.seed(999)
    pred1 <- predict_mspm(fit = fit, newdata = data)
    set.seed(999)
    pred2 <- predict_mspm(fit = fit, newdata = data)

    # Compare labels.
    if (!identical(predictedLabels(pred1), predictedLabels(pred2))) {
        stop("Test failed: Predictions are not reproducible with same seed and data.")
    }

    # Compare label indexes.
    if (!identical(predictedLabelIndexes(pred1), predictedLabelIndexes(pred2))) {
        stop("Test failed: Predicted label indexes are not reproducible with same seed and data.")
    }

    cat("Test passed: Predictions are reproducible with same seed and data.\n")
}

# Test case 5: Predict with minimal data (6 observation, 1 covariate, 3x1 gammas)
.test_predict_minimal_data <- function() {
    # Generate minimal data
    data <- generate_synthetic_data(
        nobs = 6,
        ncov = 1,
        ngamma = c(1, 1, 1),
        seed = 2026
    )

    # Fit model
    fit = fit_mspm(
        data = data,
        ndraws = 10,
        burnin = 10,
        thin = 1,
        tune = 0.1,
        seed = 2026,
        verbose = 0
    )

    # Predict
    pred <- predict_mspm(fit = fit, newdata = data)
    labels <- predictedLabels(pred)

    # Check that number of ylabels elements match number of targets
    if (length(labels) != ntargets(data)) {
        stop("Test failed: Number of predicted label sets does not match number of targets (minimal data).")
    }

    # Check that the predicted labels have correct shape
    for (i in 1:ntargets(data)) {
        if (nrow(labels[[i]]) != 6) {
            stop(paste("Test failed: Number of predicted labels does not match 6 observation for target", i, "(minimal data)"))
        }
        if (ncol(labels[[i]]) != ndraws(fit)) {
            stop(paste("Test failed: Number of predicted label draws does not match number of posterior draws for target", i, "(minimal data)"))
        }
    }

    cat("Test passed: Prediction works for minimal data (6 obs, 1 covariate).\n")
}


# Test case 6: Test if accessor functions work for mspm_labeled_prediction object.
.test_accessors_mspm_labeled_prediction <- function() {
    # Generate the data.
    data <- generate_synthetic_data(
        nobs = 100,
        ncov = 5,
        ngamma = c(2, 3),
        seed = 42
    )

    # Fit using all the data.
    fit = fit_mspm(
        data = data,
        ndraws = 1000,
        burnin = 500,
        thin = 2,
        tune = 0.1,
        seed = 1234,
        verbose = 0
    )

    # Predict labels.
    prediction <- predict_mspm(
        fit = fit,
        newdata = data
    )
    labels = prediction$ylabels
    labelIndexes = prediction$ylabelIndexes

    # Test ntargets accessor
    if (ntargets(prediction) != ntargets(data)) {
        stop("Test failed: ntargets accessor does not return correct value. Got ", 
             ntargets(prediction), " expected ", ntargets(data), ".")
    }

    # Test nlevels accessor
    if (!all(nlevels(prediction) == nlevels(data))) {
        stop("Test failed: nlevels accessor does not return correct value.")
    }

    # Test predictorNames accessor.
    if (!all.equal(predictorNames(prediction), predictorNames(data))) {
        stop("Test failed: predictorNames accessor returned incorrect value.")
    }

    # Test responseNames accessor.
    if (!all.equal(responseNames(prediction), responseNames(data))) {
        stop("Test failed: responseNames accessor returned incorrect value.")
    }

    # Test levelNames accessor.
    if (!all.equal(levelNames(prediction), levelNames(data))) {
        stop("Test failed: levelNames accessor returned incorrect value.")
    }

    # Test ndraws accessor.
    if (ndraws(prediction) != ndraws(fit)) {
        stop("Test failed: ndraws accessor does not return correct value.")
    }

    # Test predictedLabels accessor.
    if (!identical(predictedLabels(prediction), labels)) {
        stop("Test failed: predictedLabels accessor does not return correct labels.")
    }

    # Test predictedLabelIndexes accessor.
    if (!identical(predictedLabelIndexes(prediction), labelIndexes)) {
        stop("Test failed: predictedLabelIndexes accessor does not return correct label indexes.")
    }

    cat("Test passed: Accessor functions for mspm_labeled_prediction work correctly.\n")
}