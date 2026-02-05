source("R/internal.R")

# The main predict method for mspm objects.
#
# Arguments:
# fit: A fitted mspm object.
# newdata: An mspm_data object containing new data for prediction. If NULL, uses data from fit.
# type: The type of prediction to perform. Currently only "draws" is supported.
# latentOnly: If TRUE, only returns latent variable predictions without discretizing to labels.
#
# Returns:
# An mspm_latent_prediction object if latentOnly is TRUE, otherwise an mspm_labeled_prediction 
# object.
predict_mspm <- function(
    fit,
    ...,
    newdata = NULL,
    type = "draws",
    latentOnly = FALSE
) {
    # If newdata is null, use data from fit.
    if (is.null(newdata)) {
        newdata <- fit$data
    }

    if (type == "draws") {

        # Predict latent variable draws.
        ystars <- .predict_mspm_latent_draws(fit, newdata)
        latentPredictions <- new_mspm_latent_prediction(
            data = newdata,
            fit = fit,
            ystars = ystars,
            call = match.call()
        )

        # Return latent predictions only if specified.
        if (latentOnly) {
            return(latentPredictions)
        }

        # Predict labels as well.
        labels <- .predict_mspm_labels_draws(fit, newdata, ystars)
        return (new_mspm_labeled_prediction(
            data = newdata,
            fit = fit,
            ylabels = labels$ylabels,
            ylabelIndexes = labels$ylabelIndexes,
            latentPredictions = latentPredictions,
            call = match.call()
        ))
    }
    else {
        stop("Currently only 'draws' type is supported.")
    }
}

# Predict the latent variable draws for new data. Prediction is done for each
# posterior draw in the fit.
# 
# Arguments:
# fit: A fitted mspm object.
# newdata: An mspm_data object containing new data for prediction.
#
# Returns:
# A list of matrices of latent variable samples for each dataset. Each matrix has rows
# corresponding to observations and columns to posterior draws.
.predict_mspm_latent_draws <- function(
    fit,
    newdata
) {
    ystars = list()
    for (i in 1:newdata$ntargets) {
        Xnew <- newdata$Xlist[[i]]
        nobs_new <- nrow(Xnew)
        ndraws <- nrow(fit$beta)

        # Compute the linear predictor for all draws
        ystar <- tcrossprod(Xnew, fit$beta)

        # Store the latent variable samples in the matrix.
        ystars[[i]] <- ystar
    }
    return (ystars)
}

# Predict the labels from the latent variable draws for new data.
#
# Arguments:
# fit: A fitted mspm object.
# newdata: An mspm_data object containing new data for prediction.
# ystars: A list of matrices of latent variable samples for each dataset. 
#
# Returns:
# A list of predicted labels for each dataset. Each element is a matrix where rows correspond
# to observations and columns to posterior draws. Each element in the matrix is a predicted 
# label.
.predict_mspm_labels_draws <- function(
    fit,
    newdata,
    ystars
) {
    ntargets <- ntargets(fit)
    ndraws <- ndraws(fit)

    # Label each dataset
    ylabels = list()
    ylabelIndexes = list()
    for (i in 1:ntargets) {
        ystar <- ystars[[i]]

        # Discretize using sampled gammas for each draw
        y <- sapply(1:ndraws, function(draw) {
            return(findInterval(ystar[, draw], fit$gammas[[i]][draw, ]))
        })

        # Store the label indexes (1-based)
        ylabelIndexes[[i]] <- y + 1
        
        # Convert to factor labels using level names.
        levelNames <- levelNames(fit)[[i]]
        labels <- apply(y, c(1,2), function(idx) {
            return(levelNames[idx+1])
        })
        ylabels[[i]] <- labels
    }

    return (list(
        ylabels = ylabels,
        ylabelIndexes = ylabelIndexes
    ))
}