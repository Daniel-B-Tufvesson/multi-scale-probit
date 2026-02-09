
#' The data specification for a multi-scale probit data object.
#'
#' @param object One of:
#' \itemize{
#'   \item An mspm object.
#'   \item An mspm_data object.
#'   \item An mspm_labeled_prediction object.
#'   \item An mspm_labeled_evaluation object.
#'   \item An mspm_cv_result object.
#' }
#' @param ... Additional arguments (not used).
#' @return An mspm_data_spec object containing the data specification for the data.
data_spec <- function(object, ...) {
    UseMethod("data_spec")
}

#' The number of targets (datasets) in a multi-scale probit model object.
#'
#' @param object One of:
#' \itemize{
#'   \item An mspm object.
#'   \item An mspm_data object.
#'   \item An mspm_latent_prediction object.
#'   \item An mspm_labeled_prediction object.
#'   \item An mspm_labeled_evaluation object.
#' }
#' @param ... Additional arguments (not used).
#' @return An integer indicating the number of target datasets.
ntargets <- function(object, ...) {
    UseMethod("ntargets")
}

#' The number of levels for each target dataset in a multi-scale probit model object.
#'
#' @param object One of:
#' \itemize{
#'   \item An mspm object.
#'   \item An mspm_data object.
#'   \item An mspm_latent_prediction object.
#'   \item An mspm_labeled_prediction object.
#'   \item An mspm_labeled_evaluation object.
#' }
#' @param ... Additional arguments (not used).
#' @return A vector of integers indicating the number of levels for each target dataset.
nlevels <- function(object, ...) {
    UseMethod("nlevels")
}

#' The names of predictor variables in a multi-scale probit model object.
#'
#' @param object One of:
#' \itemize{
#'   \item An mspm object.
#'   \item An mspm_data object.
#'   \item An mspm_latent_prediction object.
#'   \item An mspm_labeled_prediction object.
#'   \item An mspm_labeled_evaluation object.
#' }
#' @param ... Additional arguments (not used).
#' @return A character vector of predictor variable names.
predictorNames <- function(object, ...) {
    UseMethod("predictorNames")
}

#' The names of response variables in a multi-scale probit model object.
#'
#' @param object One of:
#' \itemize{
#'   \item An mspm object.
#'   \item An mspm_data object.
#'   \item An mspm_latent_prediction object.
#'   \item An mspm_labeled_prediction object.
#'   \item An mspm_labeled_evaluation object.
#' }
#' @param ... Additional arguments (not used).
#' @return A character vector of response variable names.
responseNames <- function(object, ...) {
    UseMethod("responseNames")
}

#' The names of levels for each target dataset in a multi-scale probit model object.
#'
#' @param object One of:
#' \itemize{
#'   \item An mspm object.
#'   \item An mspm_data object.
#'   \item An mspm_latent_prediction object.
#'   \item An mspm_labeled_prediction object.
#'   \item An mspm_labeled_evaluation object.
#' }
#' @param ... Additional arguments (not used).
#' @return A list of character vectors containing level names for each target dataset.
levelNames <- function(object, ...) {
    UseMethod("levelNames")
}

beta <- function(object, ...) {
    UseMethod("beta")
}

gammas <- function(object, ...) {
    UseMethod("gammas")
}

meanPrior <- function(object, ...) {
    UseMethod("meanPrior")
}

precPrior <- function(object, ...) {
    UseMethod("precPrior")
}

#' The number of posterior draws in a multi-scale probit model object.
#'
#' @param object One of:
#' \itemize{
#'   \item An mspm object.
#'   \item An mspm_labeled_prediction object.
#'   \item An mspm_labeled_evaluation object.
#' }
#' @param ... Additional arguments (not used).
#' @return An integer indicating the number of posterior draws.
ndraws <- function(object, ...) {
    UseMethod("ndraws")
}

#' The number of posterior draws before thinning in a multi-scale probit model object.
#'
#' @param object An mspm object.
#' @param ... Additional arguments (not used).
#' @return An integer indicating the number of posterior draws before thinning.
ndrawsNoThin <- function(object, ...) {
    UseMethod("ndrawsNoThin")
}

burnin <- function(object, ...) {
    UseMethod("burnin")
}

thin <- function(object, ...) {
    UseMethod("thin")
}

#' The predicted values of latent variable. This is also known as y* (y-star).
#'
#' @param object An mspm_labeled_prediction object.
#' @param ... Additional arguments (not used).
#' @return A list of matrices containing predicted latent variable values for each target dataset.
latent <- function(object, ...) {
    UseMethod("latent")
}

#' The predicted labels for each target dataset. 
#'
#' @param object One of:
#' \itemize{
#'   \item An mspm_labeled_prediction object.
#' }
#' @param ... Additional arguments (not used).
#' @return A list of matrices containing predicted labels for each target dataset.
predictedLabels <- function(object, ...) {
    UseMethod("predictedLabels")
}

#' The indexes of predicted labels for each target dataset.
#'
#' @param object One of:
#' \itemize{
#'   \item An mspm_labeled_prediction object.
#' }
#' @param ... Additional arguments (not used).
#' @return A list of matrices containing indexes of predicted labels for each target dataset.
predictedLabelIndexes <- function(object, ...) {
    UseMethod("predictedLabelIndexes")
}

#' The evaluation metrics for predictions.
#'
#' @param object One of:
#' \itemize{
#'   \item An mspm_labeled_evaluation object.
#'   \item An mspm_cv_result object.
#' }
#' @param ... Additional arguments (not used).
#' @return A character vector of evaluation metric names.
evalMetrics <- function(object, ...) {
    UseMethod("evalMetrics")
}

#' The results of evaluation for each posterior draw.
#'
#' @param object One of:
#' \itemize{
#'   \item An mspm_labeled_evaluation object.
#' }
#' @param ... Additional arguments (not used).
#' @return A list containing the computed evaluation results. Each element corresponds 
#' to a metric and contains a matrix of results with rows for draws and columns 
#' for targets (and possibly harmonic mean as the last column).
evalDrawResults <- function(object, ...) {
    UseMethod("evalDrawResults")
}

#' The mean evaluation results for each target across all draws.
#'
#' @param object One of:
#' \itemize{
#'   \item An mspm_labeled_evaluation object.
#' }
#' @param ... Additional arguments (not used).
#' @return A list containing the mean evaluation results for each target (and possibly 
#' harmonic mean) across all draws.
evalTargetMeans <- function(object, ...) {
    UseMethod("evalTargetMeans")
}

#' The mean evaluation results for each draw across all targets.
#'
#' @param object One of:
#' \itemize{
#'   \item An mspm_labeled_evaluation object.
#' }
#' @param ... Additional arguments (not used).
#' @return A list containing the mean evaluation results for each draw across all targets.
evalDrawMeans <- function(object, ...) {
    UseMethod("evalDrawMeans")
}

#' The mean evaluation results for each metric for each draw and target.
#'
#' @param object One of:
#' \itemize{
#'   \item An mspm_labeled_evaluation object.
#' }
#' @param ... Additional arguments (not used).
#' @return A matrix containing the mean over each metric for each draw and target (and 
#' possibly harmonic mean).
evalMetricMeans <- function(object, ...) {
    UseMethod("evalMetricMeans")
}

nsplits <- function(object, ...) {
    UseMethod("nsplits")
}

cvAllEvaluations <- function(object, ...) {
    UseMethod("cvAllEvaluations")
}

cvMeans <- function(object, ...) {
    UseMethod("cvMeans")
}

#' All the evaluation draws aggregated across all cross-validation splits.
#'
#' @param object An mspm_cv_result object.
#' @param ... Additional arguments (not used).
#' @return A list containing all the evaluation draws aggregated across all cross-validation splits.
#' Each element corresponds to a metric and contains a matrix of results with rows for draws and 
#' columns for targets (and possibly harmonic mean as the last column). Returns NULL if the cv result 
#' does not contain evaluation draws.
cvAllDraws <- function(object, ...) {
    UseMethod("cvAllDraws")
}