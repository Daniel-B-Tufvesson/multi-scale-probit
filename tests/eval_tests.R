
source("R/predict.R")
source("R/eval.R")

run_all_eval_tests <- function() {
    .test_predict_and_eval_draws()
    .test_eval_reproducibility()
    .test_only_f1_metric()
    .test_only_kendall()
    .test_no_metric()
    .test_unsupported_metric()
    .test_accessors()
}


# Test case 1: Generate synthetic data, fit a model, make predictions, and evaluate 
# predictions.
.test_predict_and_eval_draws <- function() {
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
    mspm.pred <- predict_mspm(
        fit = mspm.fit
    )

    # Evaluate predictions.
    metrics = c("f1", "kendall")
    mspm.eval <- eval_mspm_prediction_draws(
        predictions = mspm.pred,
        metrics = metrics
    )

    # Check that prediction object is identical.
    if (!identical(mspm.eval$prediction, mspm.pred)) {
        stop("Test failed: Prediction object in evaluation does not match original predictions.")
    }
    # Check that metrics are identical.
    if (!identical(mspm.eval$metrics, metrics)) {
        stop("Test failed: Metrics in evaluation do not match requested metrics. Requested metric: ", 
             paste(metrics, collapse = ", "), ". Found metrics: ", paste(mspm.eval$metrics, collapse = ", "), ".")
    }

    # Check that there is a draw result for each metric.
    if (length(mspm.eval$drawResults) != length(metrics)) {
        stop("Test failed: Number of draw results does not match number of metrics. Expected number of metrics: ", 
             length(metrics), ". Found number of draw results: ", length(mspm.eval$drawResults), ".")
    }

    # Check that the number of draws in result match the number of posterior draws, and 
    # check that number of targets match.
    for (metric in metrics) {
        actualDraws = nrow(mspm.eval$drawResults[[metric]])
        if (actualDraws != mspm.fit$ndraws) {
            stop(paste("Test failed: Number of draws in result for metric", metric, 
                       "does not match number of posterior draws. Expected:", 
                       mspm.fit$ndraws, 
                       "Got:", 
                       actualDraws))
        }

        # Need to accomodate harmonic mean if more than 1 target.
        if (mspm.data$ntargets > 1) {
            expected_ncols <- mspm.data$ntargets + 1
        } else {
            expected_ncols <- mspm.data$ntargets
        }
        # One column for each target + one for harmonic mean.
        if (ncol(mspm.eval$drawResults[[metric]]) != expected_ncols) {
            stop(paste("Test failed: Number of targets in result for metric", metric, "does not match number of targets in data."))
        }
    }

    # Check that there is a mean for each target.
    for (metric in metrics) {
        actualTargets = length(mspm.eval$targetMeans[[metric]])
        # Need to accomodate harmonic mean if more than 1 target.
        if (mspm.data$ntargets > 1) {
            expectedTargets <- mspm.data$ntargets + 1
        } else {
            expectedTargets <- mspm.data$ntargets
        }
        if (actualTargets != expectedTargets) {
            stop(paste("Test failed: Number of target means for metric", metric, 
                       "does not match number of targets. Expected:", 
                       expectedTargets, 
                       "Got:", 
                       actualTargets))
        }
    }

    # Check that there is a mean for each draw.
    for (metric in metrics) {
        actualDraws = length(mspm.eval$drawMeans[[metric]])
        if (actualDraws != mspm.fit$ndraws) {
            stop(paste("Test failed: Number of draw means for metric", metric, 
                       "does not match number of posterior draws. Expected:", 
                       mspm.fit$ndraws, 
                       "Got:", 
                       actualDraws))
        }
    }

    # Check that there is a mean for each draw across the metrics.
    actualMetricMeansRows = nrow(mspm.eval$metricMeans)
    actualMetricMeansCols = ncol(mspm.eval$metricMeans)
    if (actualMetricMeansRows != mspm.fit$ndraws) {
        stop(paste("Test failed: Number of rows in metric means does not match number of posterior draws. Expected:", 
                   mspm.fit$ndraws, 
                   "Got:", 
                   actualMetricMeansRows))
    }
    # Need to accomodate harmonic mean if more than 1 target.
    if (mspm.data$ntargets > 1) {
        expectedMetricMeansCols <- mspm.data$ntargets + 1
    } else {
        expectedMetricMeansCols <- mspm.data$ntargets
    }
    if (actualMetricMeansCols != expectedMetricMeansCols) {
        stop(paste("Test failed: Number of columns in metric means does not match number of targets. Expected:", 
                   expectedMetricMeansCols, 
                   "Got:", 
                   actualMetricMeansCols))
    }

    cat("Test passed: Evaluation metrics computed correctly for all targets.\n")
}

# Test case 2: Make two evaluations on the same predictions and ensure results are identical.
.test_eval_reproducibility <- function() {
    # Generate the data.
    mspm.data <- generate_synthetic_data(
        nobs = 50,
        ncov = 20,
        ngamma = c(2, 2),
        seed = 5678
    )

    # Fit model.
    mspm.fit = fit_mspm(
        data = mspm.data,
        ndraws = 300,
        burnin = 300,
        thin = 1,
        tune = 0.1,
        seed = 5678,
        verbose = 0
    )

    # Predict labels.
    mspm.pred <- predict_mspm(
        fit = mspm.fit
    )

    # Evaluate predictions twice.
    metrics = c("f1", "kendall")
    mspm.eval1 <- eval_mspm_prediction_draws(
        predictions = mspm.pred,
        metrics = metrics
    )
    mspm.eval2 <- eval_mspm_prediction_draws(
        predictions = mspm.pred,
        metrics = metrics
    )

    # Check that evaluations are identical.
    if (!identical(mspm.eval1, mspm.eval2)) {
        stop("Test failed: Evaluations from repeated calls are not identical.")
    }

    cat("Test passed: Evaluations are reproducible across multiple calls.\n")
}

# Test case 3: Evaluate only for F1 metric.
.test_only_f1_metric <- function() {
    # Generate the data.
    mspm.data <- generate_synthetic_data(
        nobs = 80,
        ncov = 30,
        ngamma = c(3, 2),
        seed = 91011
    )

    # Fit model.
    mspm.fit = fit_mspm(
        data = mspm.data,
        ndraws = 400,
        burnin = 400,
        thin = 1,
        tune = 0.1,
        seed = 91011,
        verbose = 0
    )

    # Predict labels.
    mspm.pred <- predict_mspm(
        fit = mspm.fit
    )

    # Evaluate predictions for only F1 metric.
    metrics = c("f1")
    mspm.eval <- eval_mspm_prediction_draws(
        predictions = mspm.pred,
        metrics = metrics
    )

    # Check that only F1 metric is present in results.
    if (!identical(names(mspm.eval$drawResults), metrics)) {
        stop("Test failed: Evaluation results contain metrics other than F1.")
    }

    # Check that there is a draw result for F1 metric.
    if (length(mspm.eval$drawResults) != 1 || is.null(mspm.eval$drawResults[["f1"]])) {
        stop("Test failed: F1 metric results are missing in evaluation.")
    }

    cat("Test passed: Evaluation correctly computes only the specified F1 metric.\n")
}

# Test case 4: Evaluate only for Kendall metric.
.test_only_kendall <- function() {
    # Generate the data.
    mspm.data <- generate_synthetic_data(
        nobs = 80,
        ncov = 30,
        ngamma = c(3, 2),
        seed = 91011
    )

    # Fit model.
    mspm.fit = fit_mspm(
        data = mspm.data,
        ndraws = 400,
        burnin = 400,
        thin = 1,
        tune = 0.1,
        seed = 91011,
        verbose = 0
    )

    # Predict labels.
    mspm.pred <- predict_mspm(
        fit = mspm.fit
    )

    # Evaluate predictions for only Kendall metric.
    metrics = c("kendall")
    mspm.eval <- eval_mspm_prediction_draws(
        predictions = mspm.pred,
        metrics = metrics
    )

    # Check that only Kendall metric is present in results.
    if (!identical(names(mspm.eval$drawResults), metrics)) {
        stop("Test failed: Evaluation results contain metrics other than Kendall.")
    }

    # Check that there is a draw result for Kendall metric.
    if (length(mspm.eval$drawResults) != 1 || is.null(mspm.eval$drawResults[["kendall"]])) {
        stop("Test failed: Kendall metric results are missing in evaluation.")
    }

    cat("Test passed: Evaluation correctly computes only the specified Kendall metric.\n")
}

# Test case 5: Evaluate with no metrics specified. We expect an error.
.test_no_metric <- function() {
    # Generate the data.
    mspm.data <- generate_synthetic_data(
        nobs = 80,
        ncov = 30,
        ngamma = c(3, 2),
        seed = 91011
    )

    # Fit model.
    mspm.fit = fit_mspm(
        data = mspm.data,
        ndraws = 400,
        burnin = 400,
        thin = 1,
        tune = 0.1,
        seed = 91011,
        verbose = 0
    )

    # Predict labels.
    mspm.pred <- predict_mspm(
        fit = mspm.fit
    )

    # Evaluate predictions with no metrics specified.
    error_caught <- FALSE
    error_message <- NULL
    tryCatch(
        {
            mspm.eval <- eval_mspm_prediction_draws(
                predictions = mspm.pred,
                metrics = c()
            )
        },
        error = function(e) {
            error_caught <<- TRUE
            error_message <<- e$message
        }
    )

    if (!error_caught) {
        stop("Test failed: No error thrown when evaluating with no metrics specified.")
    }
    if (!grepl("At least one metric must be specified.", error_message, ignore.case = TRUE)) {
        stop("Test failed: Incorrect error message. Got: ", error_message)
    }

    cat("Test passed: Error correctly thrown when no metrics are specified for evaluation.\n")
}

# Test case 6: Evaluate with an unsupported metric. We expect an error.
.test_unsupported_metric <- function () {
    # Generate the data.
    mspm.data <- generate_synthetic_data(
        nobs = 80,
        ncov = 30,
        ngamma = c(3, 2),
        seed = 91011
    )

    # Fit model.
    mspm.fit = fit_mspm(
        data = mspm.data,
        ndraws = 400,
        burnin = 400,
        thin = 1,
        tune = 0.1,
        seed = 91011,
        verbose = 0
    )

    # Predict labels.
    mspm.pred <- predict_mspm(
        fit = mspm.fit
    )

    # Evaluate predictions with an unsupported metric.
    error_caught <- FALSE
    error_message <- NULL
    tryCatch(
        {
            mspm.eval <- eval_mspm_prediction_draws(
                predictions = mspm.pred,
                metrics = c("unsupported_metric")
            )
        },
        error = function(e) {
            error_caught <<- TRUE
            error_message <<- e$message
        }
    )

    if (!error_caught) {
        stop("Test failed: No error thrown when evaluating with an unsupported metric.")
    }
    if (!grepl("Unsupported metric: unsupported_metric", error_message, ignore.case = TRUE)) {
        stop("Test failed: Incorrect error message. Got: ", error_message)
    }

    cat("Test passed: Error correctly thrown when unsupported metric is specified for evaluation.\n")
}

.test_accessors <- function() {
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
    mspm.pred <- predict_mspm(
        fit = mspm.fit
    )

    # Evaluate predictions.
    metrics = c("f1", "kendall")
    mspm.eval <- eval_mspm_prediction_draws(
        predictions = mspm.pred,
        metrics = metrics
    )

    # Test ndraws accessor.
    if (ndraws(mspm.eval) != mspm.fit$ndraws) {
        stop("Test failed: ndraws accessor returned incorrect value.")
    }

    # Test ndraws accessor without thinning.
    if (ndraws(mspm.eval, withoutThinning = TRUE) != mspm.fit$ndrawsNoThin) {
        stop("Test failed: ndraws accessor without thinning returned incorrect value.")
    }

    # Test nlevels accessor.
    if (!all(nlevels(mspm.eval) == mspm.data$nlevels)) {
        stop("Test failed: nlevels accessor returned incorrect value.")
    }
    
    # Test predictorNames accessor.
    if (!identical(predictorNames(mspm.eval), mspm.data$predictorNames)) {
        stop("Test failed: predictorNames accessor returned incorrect value.")
    }

    # Test responseNames accessor.
    if (!identical(responseNames(mspm.eval), mspm.data$responseNames)) {
        stop("Test failed: responseNames accessor returned incorrect value.")
    }

    # Test levelNames accessor.
    if (!identical(levelNames(mspm.eval), mspm.data$levelNames)) {
        stop("Test failed: levelNames accessor returned incorrect value.")
    }

    # Test model accessor.
    if (!identical(model(mspm.eval), mspm.fit)) {
        stop("Test failed: model accessor returned incorrect value.")
    }

    # Test predictedLabels accessor.
    if (!identical(predictedLabels(mspm.eval), mspm.pred$ylabels)) {
        stop("Test failed: predictedLabels accessor returned incorrect value.")
    }

    # Test latent accessor.
    if (!identical(latent(mspm.eval), mspm.pred$latentPredictions$ystars)) {
        stop("Test failed: latent accessor returned incorrect value.")
    }

    # Test predictedLabelIndexes accessor.
    if (!identical(predictedLabelIndexes(mspm.eval), mspm.pred$ylabelIndexes)) {
        stop("Test failed: predictedLabelIndexes accessor returned incorrect value.")
    }

    # Test evalMetrics accessor.
    if (!identical(evalMetrics(mspm.eval), metrics)) {
        stop("Test failed: metrics accessor returned incorrect value.")
    }

    # Test evalDrawResults accessor.
    if (!identical(evalDrawResults(mspm.eval), mspm.eval$drawResults)) {
        stop("Test failed: drawResults accessor returned incorrect value.")
    }

    # Test evalTargetMeans accessor.
    if (!identical(evalTargetMeans(mspm.eval), mspm.eval$targetMeans)) {
        stop("Test failed: targetMeans accessor returned incorrect value.")
    }

    # Test evalDrawMeans accessor.
    if (!identical(evalDrawMeans(mspm.eval), mspm.eval$drawMeans)) {
        stop("Test failed: drawMeans accessor returned incorrect value.")
    }

    # Test evalMetricMeans accessor.
    if (!identical(evalMetricMeans(mspm.eval), mspm.eval$metricMeans)) {
        stop("Test failed: metricMeans accessor returned incorrect value.")
    }

    cat("Test passed: Accessor functions return correct values.\n")
}