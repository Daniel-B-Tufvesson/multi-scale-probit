
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

predictorNames <- function(object, ...) {
    UseMethod("predictorNames")
}

responseNames <- function(object, ...) {
    UseMethod("responseNames")
}

levelNames <- function(object, ...) {
    UseMethod("levelNames")
}
