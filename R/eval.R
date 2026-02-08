library(callr)

# Evaluate the predictions from an mspm model. 
#
# Supported metrics are "f1" for F1 score and "kendall" for Kendall's tau correlation.
#
# Arguments:
# predictions: An mspm_labeled_prediction object containing predicted labels.
# test_data: An mspm_data object containing the test data with true labels.
# metrics: A character vector specifying which evaluation metrics to compute. Defaults to 
#          c("f1", "kendall").
# harmonic: Logical indicating whether to compute harmonic mean of metrics when multiple
#           targets are present. Defaults to TRUE.
#
# Returns:
# An mspm_labeled_evaluation object containing the evaluation results.
eval_mspm_prediction_draws <- function(
    predictions,
    test_data,
    ...,
    metrics = c("f1", "kendall")
) {
    .validate_metrics(metrics)

    # Extract necessary components from the labeled predictions.
    ylabels_pred <- predictedLabelIndexes(predictions)
    ylist_true <- test_data$ylist
    ntargets <- ntargets(test_data)
    ndraws <- ndraws(predictions)
    nlevels <- nlevels(test_data)

    # Determine the number of result columns (for multi-scale, add one)
    # Last column is used for harmonic mean accross datasets.
    if(ntargets > 1) {
        nresults <- ntargets + 1
    } else {
        nresults <- ntargets
    }

    # Compute each metric for each target.
    results <- list()
    for (metric in metrics) {

        # Create the result matrix for this metric.
        res <- matrix(0, nrow = ndraws, ncol = nresults)
        colnames(res) <- if (ntargets > 1) {
            c(paste0("D", 1:ntargets), "HarmonicMean")
        } else {
            paste0("D", 1:ntargets)
        }

        for (i in 1:ntargets) {
            y_true <- ylist_true[[i]]
            y_pred_draws <- ylabels_pred[[i]]

            # Compute metric for each draw.
            if (metric == "f1") {
                res[, i] <- compute_f1_score_draws(y_true, y_pred_draws, nlevels[i])
            } 
            else if (metric == "kendall") {
                res[, i] <- apply(y_pred_draws, 2, function(y_pred){
                    return(suppressWarnings(cor(y_true, y_pred, method = "kendall")))
                })
            } 
            else {
                stop(paste("Unsupported metric:", metric))
            }
        }

        # Compute harmonic mean over datasets.
        if (ntargets > 1) {
            res[, ncol(res)] <- apply(res[, 1:(ncol(res)-1)], 1, function(x) {
                 if (any(x == 0, na.rm = TRUE)) {
                    return(0)
                } else {
                    return(length(x) / sum(1/x, na.rm = TRUE))
                }
            })
        } 

        results[[metric]] <- res
    }

    # Compute harmonic mean results for each target.
    targetMeans <- list()
    for (metric in metrics) {
        res <- results[[metric]]
        means <- colMeans(res)
        targetMeans[[metric]] <- means
    }

    # Compute harmonic mean results for each draw.
    drawMeans <- list()
    for (metric in metrics) {
        res <- results[[metric]]
        means <- rowMeans(res)
        drawMeans[[metric]] <- means
    }

    # Compute harmonic mean over the metrics for each draw.
    metricMeans <- matrix(0, nrow = ndraws, ncol = nresults)
    colnames(metricMeans) <- if (ntargets > 1) {
        c(paste0("D", 1:ntargets), "HarmonicMean")
    } else {
        paste0("D", 1:ntargets)
    }
    for (j in 1:nresults) { # for each dataset/column
        # Create a matrix where each row is a draw, each column is a metric
        metric_mat <- sapply(results, function(res) res[, j])
        # Harmonic mean for each draw in this dataset
        metricMeans[, j] <- apply(metric_mat, 1, function(x) {
             if (any(x == 0, na.rm = TRUE)) {
                return(0)
            } else {
                return(length(x) / sum(1/x, na.rm = TRUE))
            }
        })
    }

    return(new_mspm_labeled_evaluation(
        data_spec = data_spec(predictions),
        ndraws = ndraws(predictions),
        metrics = metrics,
        drawResults = results,
        targetMeans = targetMeans,
        drawMeans = drawMeans,
        metricMeans = metricMeans,
        call = match.call()
    ))
}

# Compute the F1 score given true and predicted labels.
#
# Arguments:
# y_true: A vector of true labels.
# y_pred: A matrix of predicted labels, where each row is an observation and each column 
#         is a draw.
# nlabels: The number of unique labels.
# 
# Returns:
# A vector of F1 scores for each draw.
compute_f1_score_draws <- function(y_true, y_pred, nlabels) {
    # Check arguments are valid.
    if (is.matrix(y_pred) == FALSE) {
        stop("y_pred must be a matrix with rows as observations and columns as draws.")
    }
    if (is.vector(y_true) == FALSE) {
        stop("y_true must be a vector of true labels.")
    }
    if (nrow(y_pred) != length(y_true)) {
        stop("Number of rows in y_pred must match length of y_true.")
    }
    if (any(y_pred < 1) || any(y_pred > nlabels)) {
        stop("y_pred contains labels outside the valid range.")
    }

    # Call the cpp backend as a subprocess.
    f1 <- tryCatch({callr::r(
        function(y_pred, y_true, nlabels) {
            devtools::load_all()
            cpp_fmeasure_distribution(y_pred, y_true, nlabels)
        },
        args = list(
            y_pred = y_pred-1, # 0-based for cpp
            y_true = y_true-1, # 0-based for cpp
            nlabels = nlabels
        )
    )}, error = function(e) {
        message("Error in cpp_hprobit: ", e$message)
        if (!is.null(e$stdout)) {
            cat("---- STDOUT ----\n")
            cat(e$stdout, sep = "\n")
        }
        if (!is.null(e$stderr)) {
            cat("---- STDERR ----\n")
            cat(e$stderr, sep = "\n")
        }
        stop(e)
    })
    return(f1)
}


# Validate that the specified metrics are supported.
#
# Arguments:
# metrics: A character vector of metric names to validate.
.validate_metrics <- function(metrics) {
    supported_metrics <- c("f1", "kendall")
    for (metric in metrics) {
        if (!(metric %in% supported_metrics)) {
            stop(paste("Unsupported metric:", metric))
        }
    }

    if (length(metrics) == 0) {
        stop("At least one metric must be specified.")
    }
}

