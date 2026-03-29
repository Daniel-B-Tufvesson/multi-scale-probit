
run_all_eval_tests <- function() {
    .test_predict_and_eval_draws()
    .test_eval_reproducibility()
    .test_only_f1_metric()
    .test_only_kendall()
    .test_no_metric()
    .test_unsupported_metric()
    .test_eval_accessors()
}


# Test case 1: Generate synthetic data, fit a model, make predictions, and evaluate predictions.
.test_predict_and_eval_draws <- function() {
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
    pred <- predict_mspm(
        fit = fit,
        newdata = data
    )

    # Evaluate predictions.
    metrics = c("f1", "kendall")
    eval <- eval_mspm_prediction_draws(
        predictions = pred,
        test_data = data,
        metrics = metrics
    )
    drawResults <- get_eval_draw_results(eval)
    targetMeans <- get_eval_target_means(eval)
    drawMeans <- get_eval_draw_means(eval)
    metricMeans <- get_eval_metric_means(eval)

    # Check data_spec.
    if (!identical(get_data_spec(eval), get_data_spec(data))) {
        stop("Test failed: Data specification in evaluation does not match original data specification.")
    }

    # Check that metrics are identical.
    if (!identical(get_eval_metrics(eval), metrics)) {
        stop("Test failed: Metrics in evaluation do not match requested metrics. Requested metric: ", 
             paste(metrics, collapse = ", "), ". Found metrics: ", paste(get_eval_metrics(eval), collapse = ", "), ".")
    }

    # Check that there is a draw result for each metric.
    if (length(drawResults) != length(metrics)) {
        stop("Test failed: Number of draw results does not match number of metrics. Expected number of metrics: ", 
             length(metrics), ". Found number of draw results: ", length(drawResults), ".")
    }

    # Check correct number of draws.
    if (get_n_draws(eval) != get_n_draws(fit)) {
        stop("Test failed: Number of draws in evaluation does not match number of posterior draws. Expected:", 
             get_n_draws(fit), 
             "Got:",
             get_n_draws(eval))
    }

    # Check that the number of draws in result match the number of posterior draws, and 
    # check that number of targets match.
    for (metric in metrics) {
        actualDraws = nrow(drawResults[[metric]])
        if (actualDraws != get_n_draws(fit)) {
            stop(paste("Test failed: Number of draws in result for metric", metric, 
                       "does not match number of posterior draws. Expected:", 
                       get_n_draws(fit), 
                       "Got:", 
                       actualDraws))
        }

        # Need to accomodate harmonic mean if more than 1 target.
        if (get_n_targets(data) > 1) {
            expected_ncols <- get_n_targets(data) + 1
        } else {
            expected_ncols <- get_n_targets(data)
        }
        # One column for each target + one for harmonic mean.
        if (ncol(drawResults[[metric]]) != expected_ncols) {
            stop(paste("Test failed: Number of targets in result for metric", metric, "does not match number of targets in data."))
        }
    }

    # Check that there is a mean for each target.
    for (metric in metrics) {
        actualTargets = length(targetMeans[[metric]])
        # Need to accomodate harmonic mean if more than 1 target.
        if (get_n_targets(data) > 1) {
            expectedTargets <- get_n_targets(data) + 1
        } else {
            expectedTargets <- get_n_targets(data)
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
        actualDraws = length(drawMeans[[metric]])
        if (actualDraws != get_n_draws(fit)) {
            stop(paste("Test failed: Number of draw means for metric", metric, 
                       "does not match number of posterior draws. Expected:", 
                       get_n_draws(fit), 
                       "Got:", 
                       actualDraws))
        }
    }

    # Check that there is a mean for each draw across the metrics.
    actualMetricMeansRows = nrow(metricMeans)
    actualMetricMeansCols = ncol(metricMeans)
    if (actualMetricMeansRows != get_n_draws(fit)) {
        stop(paste("Test failed: Number of rows in metric means does not match number of posterior draws. Expected:", 
                   get_n_draws(fit), 
                   "Got:", 
                   actualMetricMeansRows))
    }
    # Need to accomodate harmonic mean if more than 1 target.
    if (get_n_targets(data) > 1) {
        expectedMetricMeansCols <- get_n_targets(data) + 1
    } else {
        expectedMetricMeansCols <- get_n_targets(data)
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
    data <- generate_synthetic_data(
        nobs = 50,
        ncov = 20,
        ngamma = c(2, 2),
        seed = 5678
    )

    # Fit model.
    fit = fit_mspm(
        data = data,
        ndraws = 300,
        burnin = 300,
        thin = 1,
        tune = 0.1,
        seed = 5678,
        verbose = 0
    )

    # Predict labels.
    pred <- predict_mspm(
        fit = fit,
        newdata = data
    )

    # Evaluate predictions twice.
    metrics = c("f1", "kendall")
    eval1 <- eval_mspm_prediction_draws(
        predictions = pred,
        test_data = data,
        metrics = metrics
    )
    eval2 <- eval_mspm_prediction_draws(
        predictions = pred,
        test_data = data,
        metrics = metrics
    )

    # Check that evaluations are identical.
    if (!identical(eval1, eval2)) {
        stop("Test failed: Evaluations from repeated calls are not identical.")
    }

    cat("Test passed: Evaluations are reproducible across multiple calls.\n")
}

# Test case 3: Evaluate only for F1 metric.
.test_only_f1_metric <- function() {
    # Generate the data.
    data <- generate_synthetic_data(
        nobs = 80,
        ncov = 30,
        ngamma = c(3, 2),
        seed = 91011
    )

    # Fit model.
    fit = fit_mspm(
        data = data,
        ndraws = 400,
        burnin = 400,
        thin = 1,
        tune = 0.1,
        seed = 91011,
        verbose = 0
    )

    # Predict labels.
    pred <- predict_mspm(
        fit = fit,
        data
    )

    # Evaluate predictions for only F1 metric.
    metrics = c("f1")
    eval <- eval_mspm_prediction_draws(
        predictions = pred,
        test_data = data,
        metrics = metrics
    )
    drawResults <- get_eval_draw_results(eval)

    # Check that only F1 metric is present in results.
    if (!identical(names(drawResults), metrics)) {
        stop("Test failed: Evaluation results contain metrics other than F1.")
    }

    # Check that there is a draw result for F1 metric.
    if (length(drawResults) != 1 || is.null(drawResults[["f1"]])) {
        stop("Test failed: F1 metric results are missing in evaluation.")
    }

    cat("Test passed: Evaluation correctly computes only the specified F1 metric.\n")
}

# Test case 4: Evaluate only for Kendall metric.
.test_only_kendall <- function() {
    # Generate the data.
    data <- generate_synthetic_data(
        nobs = 80,
        ncov = 30,
        ngamma = c(3, 2),
        seed = 91011
    )

    # Fit model.
    fit = fit_mspm(
        data = data,
        ndraws = 400,
        burnin = 400,
        thin = 1,
        tune = 0.1,
        seed = 91011,
        verbose = 0
    )

    # Predict labels.
    pred <- predict_mspm(
        fit = fit,
        newdata = data
    )

    # Evaluate predictions for only Kendall metric.
    metrics = c("kendall")
    eval <- eval_mspm_prediction_draws(
        predictions = pred,
        test_data = data,
        metrics = metrics
    )
    drawResults <- get_eval_draw_results(eval)

    # Check that only Kendall metric is present in results.
    if (!identical(names(drawResults), metrics)) {
        stop("Test failed: Evaluation results contain metrics other than Kendall.")
    }

    # Check that there is a draw result for Kendall metric.
    if (length(drawResults) != 1 || is.null(drawResults[["kendall"]])) {
        stop("Test failed: Kendall metric results are missing in evaluation.")
    }

    cat("Test passed: Evaluation correctly computes only the specified Kendall metric.\n")
}

# Test case 5: Evaluate with no metrics specified. We expect an error.
.test_no_metric <- function() {
    # Generate the data.
    data <- generate_synthetic_data(
        nobs = 80,
        ncov = 30,
        ngamma = c(3, 2),
        seed = 91011
    )

    # Fit model.
    fit = fit_mspm(
        data = data,
        ndraws = 400,
        burnin = 400,
        thin = 1,
        tune = 0.1,
        seed = 91011,
        verbose = 0
    )

    # Predict labels.
    pred <- predict_mspm(
        fit = fit,
        newdata = data
    )

    # Evaluate predictions with no metrics specified.
    error_caught <- FALSE
    error_message <- NULL
    tryCatch(
        {
            eval <- eval_mspm_prediction_draws(
                predictions = pred,
                test_data = data,
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
    data <- generate_synthetic_data(
        nobs = 80,
        ncov = 30,
        ngamma = c(3, 2),
        seed = 91011
    )

    # Fit model.
    fit = fit_mspm(
        data = data,
        ndraws = 400,
        burnin = 400,
        thin = 1,
        tune = 0.1,
        seed = 91011,
        verbose = 0
    )

    # Predict labels.
    pred <- predict_mspm(
        fit = fit,
        newdata = data
    )

    # Evaluate predictions with an unsupported metric.
    error_caught <- FALSE
    error_message <- NULL
    tryCatch(
        {
            eval <- eval_mspm_prediction_draws(
                predictions = pred,
                test_data = data,
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

.test_eval_accessors <- function() {
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
        ndraws = 10,
        burnin = 10,
        thin = 1,
        tune = 0.1,
        seed = 1234,
        verbose = 0
    )

    # Predict labels.
    pred <- predict_mspm(
        fit = fit,
        newdata = data
    )

    # Evaluate predictions.
    metrics = c("f1", "kendall")
    eval <- eval_mspm_prediction_draws(
        predictions = pred,
        test_data = data,
        metrics = metrics
    )

    # Test data_spec
    if (!identical(get_data_spec(eval), get_data_spec(data))) {
        stop("Test failed: data_spec accessor does not return correct value.")
    }

    # Test ndraws accessor.
    if (get_n_draws(eval) != 10) {
        stop("Test failed: ndraws accessor returned incorrect value.")
    }

    # Test ntargets accessor.
    if (get_n_targets(eval) != get_n_targets(data)) {
        stop("Test failed: ntargets accessor returned incorrect value.")
    }

    # Test nlevels accessor.
    if (!all(get_n_levels(eval) == get_n_levels(data))) {
        stop("Test failed: nlevels accessor returned incorrect value.")
    }
    
    # Test predictorNames accessor.
    if (!identical(get_predictor_names(eval), get_predictor_names(data))) {
        stop("Test failed: predictorNames accessor returned incorrect value.")
    }

    # Test responseNames accessor.
    if (!identical(get_response_names(eval), get_response_names(data))) {
        stop("Test failed: responseNames accessor returned incorrect value.")
    }

    # Test levelNames accessor.
    if (!identical(get_level_names(eval), get_level_names(data))) {
        stop("Test failed: levelNames accessor returned incorrect value.")
    }

    # Test evalMetrics accessor.
    if (!identical(get_eval_metrics(eval), metrics)) {
        stop("Test failed: metrics accessor returned incorrect value.")
    }

    # Test evalDrawResults accessor.
    if (!identical(get_eval_draw_results(eval), eval$drawResults)) {
        stop("Test failed: drawResults accessor returned incorrect value.")
    }

    # Test evalTargetMeans accessor.
    if (!identical(get_eval_target_means(eval), eval$targetMeans)) {
        stop("Test failed: targetMeans accessor returned incorrect value.")
    }

    # Test evalDrawMeans accessor.
    if (!identical(get_eval_draw_means(eval), eval$drawMeans)) {
        stop("Test failed: drawMeans accessor returned incorrect value.")
    }

    # Test evalMetricMeans accessor.
    if (!identical(get_eval_metric_means(eval), eval$metricMeans)) {
        stop("Test failed: metricMeans accessor returned incorrect value.")
    }

    cat("Test passed: Accessor functions return correct values.\n")
}