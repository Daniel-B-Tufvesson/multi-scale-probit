# Load and save objects.


#' Save a fitted model object to a file.
#'
#' @param fit A fitted model object to save.
#' @param file A file path to save the model to.
#' @param savePostDraws Logical; whether to retain posterior draws in the saved model. Default 
#' is TRUE.
#' @param saveTrainingData Logical; whether to retain training data in the saved model. Default 
#' is FALSE.
save_model <- function(
    fit, 
    file,
    ...,
    savePostDraws = TRUE,
    saveTrainingData = FALSE
) {
  fit <- .strip_ephemeral_state(
    fit,
    savePostDraws = savePostDraws,
    saveTrainingData = saveTrainingData
  )
  saveRDS(fit, file)
}

#' Strip the ephemeral state from a fitted model object.
#'
#' @param fit A fitted model object.
#' @param savePostDraws Logical; whether to retain posterior draws.
#' @param saveTrainingData Logical; whether to retain training data.
#' @return A copy of the fitted model object with ephemeral state removed as specified. 
.strip_ephemeral_state <- function(
    fit,
    savePostDraws,
    saveTrainingData
) {
    # Strip post draws if not saving.
    if (!savePostDraws) {
        fit$postDraws <- NULL
    }

    # Strip training data if not saving.
    if (!saveTrainingData) {
        fit$data <- NULL
    }

    fit
}

#' Load a fitted model object from a file.
#'
#' @param file A file path to load the model from.
#' @return A fitted model object loaded from the specified file path.
load_model <- function(file) {
  fit <- readRDS(file)
  .restore_model(fit)
}

#' Restore any necessary state in a fitted model object after loading from file.
#'
#' @param fit A fitted model object loaded from file.
#' @return The fitted model object with any necessary state restored. Currently this is a no-op, 
#' but it can be extended in the future if needed.
.restore_model <- function(fit) {
    fit
}

#' Save a data object to a file.
#'
#' @param data A data object to save.
#' @param file A file path to save the data to.
save_data <- function(data, file) {
  saveRDS(data, file)
}

#' Load a data object from a file.
#'
#' @param file A file path to load the data from.
#' @return A data object loaded from the specified file path.
load_data <- function(file) {
    readRDS(file)
}

