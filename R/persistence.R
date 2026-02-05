# Load and save model objects.


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

load_model <- function(file) {
  fit <- readRDS(file)
  .restore_model(fit)
}

.restore_model <- function(fit) {
    fit
}