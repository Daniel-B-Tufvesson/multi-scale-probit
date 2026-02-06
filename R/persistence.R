# Load and save objects.


#' Save a fitted model object to a file.
#'
#' @param fit A fitted model object to save.
#' @param file A file path to save the model to.
save_model <- function(
    fit, 
    file
) {
    saveRDS(fit, file)
}

#' Load a fitted model object from a file.
#'
#' @param file A file path to load the model from.
#' @return A fitted model object loaded from the specified file path.
load_model <- function(file) {
    readRDS(file)
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

