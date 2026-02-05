source("R/generics.R")

# Constructor for creating a new multi-scale probit model data structure.
#
# Arguments:
# predictorNames: A vector of predictor variable names.
# responseNames: A vector of response variable names corresponding to each dataset.
# Xlist: A list of design matrices for each dataset.
# ylist: A list of response vectors for each dataset.
# levelNames: A list of vectors containing the names of levels for each response variable.
# nlevels: A vector indicating the number of levels for each response variable.
# ntargets: The number of target datasets.
# seed: An optional random seed for reproducibility. Only used if data is synthetic.
# call: The original function call used to create the data structure.
new_mspm_data <- function(
    predictorNames,
    responseNames,
    Xlist,
    ylist,
    levelNames,
    nlevels,
    ntargets,
    seed,
    call
) {
    structure(
        list(
            predictorNames = predictorNames,
            responseNames = responseNames,
            Xlist = Xlist,
            ylist = ylist,
            levelNames = levelNames,
            nlevels = nlevels,
            ntargets = ntargets,
            seed = seed,
            call = call
        ),
        class = "mspm_data"
    )
}

ntargets.mspm_data <- function(object, ...) {
    object$ntargets
}

nlevels.mspm_data <- function(object, ...) {
    object$nlevels
}

predictorNames.mspm_data <- function(object, ...) {
    object$predictorNames
}

responseNames.mspm_data <- function(object, ...) {
    object$responseNames
}

levelNames.mspm_data <- function(object, ...) {
    object$levelNames
}

# Constructor for creating a new multi-scale probit model (MSPM).
# 
# Arguments: 
# data: An MSPM data structure holding the data used for fitting.
# beta: An MCMC object for the beta (the coefficients) parameters.
# gammas: A list of MCMC gammas (the thresholds) of each target dataset.
# meanPrior: Prior mean for regression coefficients.
# precPrior: Prior precision for regression coefficients.
# seed: The random seed used for reproducibility.
# call: The original function call used to create the model.
new_mspm <- function(
    data,
    beta,
    gammas,
    meanPrior,
    precPrior,
    seed,
    ndraws,
    ndrawsNoThin,
    thin,
    burnin,
    call
) {
    structure(
        list(
            data = data,
            beta = beta,
            gammas = gammas,
            meanPrior = meanPrior,
            precPrior = precPrior,
            call = call,
            seed = seed,
            nobs = nrow(beta), # Why is nobs the number of rows in beta? Shouldn't this be ndraws?
            ndraws = ndraws,
            ndrawsNoThin = ndrawsNoThin,
            thin = thin,
            burnin = burnin,

            training_data_info = list(
                ntargets = ntargets(data),
                nlevels = nlevels(data),
                levelNames = levelNames(data),
                predictorNames = predictorNames(data),
                responseNames = responseNames(data),
                seed = data$seed
            )
            
            # diagnostics = list(
            #     rhat = rhat,
            #     ess = ess,
            #     divergences = divergences
            # ),
        ),
        class = "mspm"
    )
}

ntargets.mspm <- function(object, ...) {
    object$training_data_info$ntargets
}

nlevels.mspm <- function(object, ...) {
    object$training_data_info$nlevels
}

levelNames.mspm <- function(object, ...) {
    object$training_data_info$levelNames
}

predictorNames.mspm <- function(object, ...) {
    object$training_data_info$predictorNames
}

responseNames.mspm <- function(object, ...) {
    object$training_data_info$responseNames
}

beta.mspm <- function(object, ...) {
    object$beta
}

gammas.mspm <- function(object, ...) {
    object$gammas
}

meanPrior.mspm <- function(object, ...) {
    object$meanPrior
}

precPrior.mspm <- function(object, ...) {
    object$precPrior
}

ndraws.mspm <- function(object, withoutThinning = FALSE, ...) {
    if (withoutThinning) {
        object$ndrawsNoThin
    } else {
        object$ndraws
    }
}

burnin.mspm <- function(object, ...) {
    object$burnin
}

thin.mspm <- function(object, ...) {
    object$thin
}


# Constructor for creating a new multi-scale probit model latent prediction object.
# 
# Arguments:
# data: The mspm_data object containing the data used for predicting (the test data).
# fit: The fitted mspm object that was used for predicting.
# ystars: A data frame of latent variable samples for each observation and each draw.
# call: The original function call used to create the prediction.
new_mspm_latent_prediction <- function(
    data,
    fit,
    ystars,
    call
) {
    structure(
        list(
            data = data,
            fit = fit,
            ystars = ystars,
            call = call
        ),
        class = "mspm_latent_prediction"
    )
}

ntargets.mspm_latent_prediction <- function(object, ...) {
    ntargets(object$fit)
}

nlevels.mspm_latent_prediction <- function(object, ...) {
    nlevels(object$fit)
}

predictorNames.mspm_latent_prediction <- function(object, ...) {
    predictorNames(object$fit)
}

responseNames.mspm_latent_prediction <- function(object, ...) {
    responseNames(object$fit)
}

levelNames.mspm_latent_prediction <- function(object, ...) {
    levelNames(object$fit)
}

ndraws.mspm_latent_prediction <- function(object, withoutThinning = FALSE, ...) {
    if (withoutThinning) {
        object$fit$ndrawsNoThin
    } else {
        object$fit$ndraws
    }
}

model.mspm_latent_prediction <- function(object, ...) {
    object$fit
}

latent.mspm_latent_prediction <- function(object, ...) {
    object$ystars
}

# Constructor for creating a new multi-scale probit model labeled prediction object.
# 
# Arguments:
# data: The mspm_data object containing the data used for predicting.
# fit: The fitted mspm object that was used for predicting.
# ylabels: A list of data frames of observed categorical labels for each observation and 
#          each draw. There is one data frame for each scale.
# ylabelIndexes: A list of matrices indicating the indexes of predicted labels for each
#                observation and each draw. There is one matrix for each scale.
# latentPredictions: The mspm_latent_prediction object containing latent variable samples.
# call: The original function call used to create the prediction.
new_mspm_labeled_prediction <- function(
    data,
    fit,
    ylabels,
    ylabelIndexes,
    latentPredictions,
    call
) {
    structure(
        list(
            data = data,
            fit = fit,
            ylabels = ylabels,
            ylabelIndexes = ylabelIndexes,
            latentPredictions = latentPredictions,
            call = call
        ),
        class = "mspm_labeled_prediction"
    )
}

ntargets.mspm_labeled_prediction <- function(object, ...) {
    ntargets(object$fit)
}

nlevels.mspm_labeled_prediction <- function(object, ...) {
    nlevels((object$fit))
}

predictorNames.mspm_labeled_prediction <- function(object, ...) {
    predictorNames(object$fit)
}

responseNames.mspm_labeled_prediction <- function(object, ...) {
    responseNames(object$fit)
}

levelNames.mspm_labeled_prediction <- function(object, ...) {
    levelNames(object$fit)
}

ndraws.mspm_labeled_prediction <- function(object, withoutThinning = FALSE, ...) {
    if (withoutThinning) {
        object$fit$ndrawsNoThin
    } else {
        object$fit$ndraws
    }
}

model.mspm_labeled_prediction <- function(object, ...) {
    object$fit
}

predictedLabels.mspm_labeled_prediction <- function(object, ...) {
    object$ylabels
}

predictedLabelIndexes.mspm_labeled_prediction <- function(object, ...) {
    object$ylabelIndexes
}

latent.mspm_labeled_prediction <- function(object, ...) {
    object$ystars
}

# Constructor for creating a new multi-scale probit model labeled evaluation object.
# 
# Arguments:
# prediction: The mspm_labeled_prediction object containing predicted labels.
# metrics: A character vector specifying which evaluation metrics were computed.
# drawResults: A list containing the computed evaluation results. Each element corresponds 
#              to a metric and contains a matrix of results with rows for draws and columns 
#              for targets (and possibly harmonic mean as the last column).
# targetMeans: A list containing the mean evaluation results for each target (and possibly 
#              harmonic mean) across all draws.
# metricMeans: A matrix containing the mean over each metric for each draw and target (and 
#              possibly harmonic mean).
# drawMeans: A list containing the mean evaluation results for each draw across all targets.
# call: The original function call used to create the evaluation.
new_mspm_labeled_evaluation <- function(
    prediction,
    metrics,
    drawResults,
    targetMeans,
    drawMeans,
    metricMeans,
    call
) {
    structure(
        list(
            prediction = prediction,
            metrics = metrics,
            drawResults = drawResults,
            targetMeans = targetMeans,
            drawMeans = drawMeans,
            metricMeans = metricMeans,
            call = call
        ),
        class = "mspm_labeled_evaluation"
    )
}

ntargets.mspm_labeled_evaluation <- function(object, ...) {
    ntargets(object$prediction$fit)
}

nlevels.mspm_labeled_evaluation <- function(object, ...) {
    nlevels(object$prediction$fit)
}

predictorNames.mspm_labeled_evaluation <- function(object, ...) {
    predictorNames(object$prediction$fit)
}

responseNames.mspm_labeled_evaluation <- function(object, ...) {
    responseNames(object$prediction$fit)
}

levelNames.mspm_labeled_evaluation <- function(object, ...) {
    levelNames(object$prediction$fit)
}

ndraws.mspm_labeled_evaluation <- function(object, withoutThinning = FALSE, ...) {
    if (withoutThinning) {
        object$prediction$fit$ndrawsNoThin
    } else {
        object$prediction$fit$ndraws
    }
}

model.mspm_labeled_evaluation <- function(object, ...) {
    object$prediction$fit
}

latent.mspm_labeled_evaluation <- function(object, ...) {
    object$prediction$latentPredictions$ystars
}

predictedLabels.mspm_labeled_evaluation <- function(object, ...) {
    object$prediction$ylabels
}

predictedLabelIndexes.mspm_labeled_evaluation <- function(object, ...) {
    object$prediction$ylabelIndexes
}

evalMetrics.mspm_labeled_evaluation <- function(object, ...) {
    object$metrics
}

evalDrawResults.mspm_labeled_evaluation <- function(object, ...) {
    object$drawResults
}

evalTargetMeans.mspm_labeled_evaluation <- function(object, ...) {
    object$targetMeans
}

evalDrawMeans.mspm_labeled_evaluation <- function(object, ...) {
    object$drawMeans
}

evalMetricMeans.mspm_labeled_evaluation <- function(object, ...) {
    object$metricMeans
}