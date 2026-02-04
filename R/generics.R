
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
#'   \item An mspm_latent_prediction object.
#'   \item An mspm_labeled_prediction object.
#' }
#' @param withoutThinning Logical indicating whether to return the number of draws without 
#' thinning (default FALSE).
#' @param ... Additional arguments (not used).
#' @return An integer indicating the number of posterior draws.
ndraws <- function(object, withoutThinning = FALSE, ...) {
    UseMethod("ndraws")
}

burnin <- function(object, ...) {
    UseMethod("burnin")
}

thin <- function(object, ...) {
    UseMethod("thin")
}

#' The underlying model used for fitting or prediction.
#'
#' @param object One of:
#' \itemize{
#'   \item An mspm_latent_prediction object.
#'   \item An mspm_labeled_prediction object.
#' }
#' @param ... Additional arguments (not used).
#' @return The underlying model object.
model <- function(object, ...) {
    UseMethod("model")
}

#' The predicted values of latent variable. This is also known as y* (y-star).
#'
#' @param object One of:
#' \itemize{
#'   \item An mspm_latent_prediction object.
#'   \item An mspm_labeled_prediction object.
#' }
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
