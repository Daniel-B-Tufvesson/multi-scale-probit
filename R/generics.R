
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
get_data_spec <- function(object, ...) {
    UseMethod("get_data_spec")
}

#' The values for the predictor variables in the data.
#'
#' @param object An mspm_data object.
#' @param ... Additional arguments (not used).
#' @return A list of matrices containing the predictor variable values for each target 
#' dataset. Each matrix has rows corresponding to observations and columns 
#' corresponding to predictor variables.
get_x_values <- function(object, ...) {
    UseMethod("get_x_values")
}

#' The values for the response variables in the data.
#'
#' @param object An mspm_data object.
#' @param ... Additional arguments (not used).
#' @return A list of vectors containing the response variable values for each target
#' dataset. Each vector has a length corresponding to the number of observations in 
#' the dataset and contains the response variable values as factors with levels 
#' corresponding to the ordinal levels of the response variable.
get_y_values <- function(object, ...) {
    UseMethod("get_y_values")
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
get_n_targets <- function(object, ...) {
    UseMethod("get_n_targets")
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
get_n_levels <- function(object, ...) {
    UseMethod("get_n_levels")
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
get_predictor_names <- function(object, ...) {
    UseMethod("get_predictor_names")
}

get_n_predictors <- function(object, ...) {
    UseMethod("get_n_predictors")
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
get_response_names <- function(object, ...) {
    UseMethod("get_response_names")
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
get_level_names <- function(object, ...) {
    UseMethod("get_level_names")
}

get_train_split <- function(object, ...) {
    UseMethod("get_train_split")
}

get_test_split <- function(object, ...) {
    UseMethod("get_test_split")
}

get_beta <- function(object, ...) {
    UseMethod("get_beta")
}

get_gammas <- function(object, ...) {
    UseMethod("get_gammas")
}


get_mean_prior <- function(object, ...) {
    UseMethod("get_mean_prior")
}

get_prec_prior <- function(object, ...) {
    UseMethod("get_prec_prior")
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
get_n_draws <- function(object, ...) {
    UseMethod("get_n_draws")
}

#' The number of posterior draws before thinning in a multi-scale probit model object.
#'
#' @param object An mspm object.
#' @param ... Additional arguments (not used).
#' @return An integer indicating the number of posterior draws before thinning.
get_n_draws_no_thin <- function(object, ...) {
    UseMethod("get_n_draws_no_thin")
}

get_burnin <- function(object, ...) {
    UseMethod("get_burnin")
}

get_thin <- function(object, ...) {
    UseMethod("get_thin")
}

#' The proposal variance used in the MCMC sampling for a multi-scale probit model object.
#' 
#' @param object One of:
#' \itemize{
#'   \item An mspm object.
#'   \item An mspm_tune_results object.
#'   \item An mspm_tune_results_pt object.
#' }
#' @param ... Additional arguments (not used).
#' @return The proposal variance on the form of:
#' \itemize{
#'   \item A vector of numerical values for each target gamma group if the sampler is the MH-Gibbs
#'   sampler.
#'   \item A list of numerical vectors for each temperature if the sampler is a parallel tempering 
#'   sampler.
#' }
get_proposal_variance <- function(object, ...) {
    UseMethod("get_proposal_variance")
}



get_n_likelihood_calls <- function(object, ...) {
    UseMethod("get_n_likelihood_calls")
}

#' The time taken for the sampling phase of MCMC sampling. This only encompasses the time taken for 
#' sampling and does not include data preprocessing, model fitting setup, burn-in, or post-processing 
#' steps. The time is measured in seconds.
#'
#' @param object An mspm object.
#' @param ... Additional arguments (not used).
#' @return A numeric value indicating the total time taken for the MCMC sampling process in seconds. 
get_sampling_time <- function(object, ...) {
    UseMethod("get_sampling_time")
}

#' The time taken for the burn-in phase of MCMC sampling. This only encompasses the time taken for
#' the burn-in phase and does not include data preprocessing, model fitting setup, sampling, or
#' post-processing steps. The time is measured in seconds.
#'
#' @param object An mspm object.
#' @param ... Additional arguments (not used).
#' @return A numeric value indicating the total time taken for the burn-in phase of MCMC sampling in
#' seconds.
get_burnin_time <- function(object, ...) {
    UseMethod("get_burnin_time")
}

#' The predicted values of latent variable. This is also known as y* (y-star).
#'
#' @param object An mspm_labeled_prediction object.
#' @param ... Additional arguments (not used).
#' @return A list of matrices containing predicted latent variable values for each target dataset.
get_latent_prediction <- function(object, ...) {
    UseMethod("get_latent_prediction")
}

#' The predicted labels for each target dataset. 
#'
#' @param object One of:
#' \itemize{
#'   \item An mspm_labeled_prediction object.
#' }
#' @param ... Additional arguments (not used).
#' @return A list of matrices containing predicted labels for each target dataset.
get_predicted_labels <- function(object, ...) {
    UseMethod("get_predicted_labels")
}

#' The indexes of predicted labels for each target dataset.
#'
#' @param object One of:
#' \itemize{
#'   \item An mspm_labeled_prediction object.
#' }
#' @param ... Additional arguments (not used).
#' @return A list of matrices containing indexes of predicted labels for each target dataset.
get_predicted_label_indexes <- function(object, ...) {
    UseMethod("get_predicted_label_indexes")
}

#' The evaluation metrics for predictions.
#'
#' @param object One of:
#' \itemize{
#'   \item An mspm_labeled_evaluation object.
#' }
#' @param ... Additional arguments (not used).
#' @return A character vector of evaluation metric names.
get_eval_metrics <- function(object, ...) {
    UseMethod("get_eval_metrics")
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
get_eval_draw_results <- function(object, ...) {
    UseMethod("get_eval_draw_results")
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
get_eval_target_means <- function(object, ...) {
    UseMethod("get_eval_target_means")
}

#' The mean evaluation results for each draw across all targets.
#'
#' @param object One of:
#' \itemize{
#'   \item An mspm_labeled_evaluation object.
#' }
#' @param ... Additional arguments (not used).
#' @return A list containing the mean evaluation results for each draw across all targets.
get_eval_draw_means <- function(object, ...) {
    UseMethod("get_eval_draw_means")
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

get_eval_metric_means <- function(object, ...) {
    UseMethod("get_eval_metric_means")
}

#' The tuning results for all cross-validation splits.
#'
#' @param object An mspm_cv_result object.
#' @param ... Additional arguments (not used).
#' @return A list containing the tuning results for all cross-validation splits. Each element of
#' the list corresponds to a cross-validation split and contains the tuning results for that split.
#' Is NULL if no tuning was done.
get_all_tune_results <- function(object, ...) {
    UseMethod("get_all_tune_results")
}

#' The number of likelihood calls made by the MCMC sampler during cross-validation. 
#'
#' @param object An mspm_cv_result object.
#' @param ... Additional arguments (not used).
#' @return An integer vector indicating the number of likelihood calls made by the MCMC sampler 
#' for each cross-validation split.
get_all_nlikelihood_calls <- function(object, ...) {
    UseMethod("get_all_nlikelihood_calls")
}

#' The MCMC diagnostics for a multi-scale probit model object.
#'
#' @param object An mspm object.
#' @param ... Additional arguments (not used).
#' @return A list containing the MCMC diagnostics for the model. The list may include diagnostics 
#' such as effective sample size (ESS) for beta coefficients and gamma thresholds, Geweke 
#' diagnostic results for beta coefficients and gamma thresholds, and any other relevant 
#' diagnostics. Returns NULL if the diagnostics cannot be computed for the given object or if the 
#' object does not contain the necessary information to compute diagnostics.
diagnostics <- function(object, ...) {
    UseMethod("diagnostics")
}

#' The effective sample size (ESS) for the beta coefficients in a multi-scale probit model object.
#'
#' @param object An mspm object.
#' @param ... Additional arguments (not used).
#' @return A vector containing the effective sample size for each beta coefficient. Returns NULL 
#' if the object does not contain beta coefficients or if the effective sample size cannot be 
#' computed.
essBeta <- function(object, ...) {
    UseMethod("essBeta")
}

#' The effective sample size (ESS) for the gamma thresholds in a multi-scale probit model object.
#'
#' @param object An mspm object.
#' @param ... Additional arguments (not used).
#' @return A list containing the effective sample size for each gamma threshold. Each element of 
#' the list corresponds to a target dataset and contains a vector of effective sample sizes for 
#' the gamma thresholds of that target. Returns NULL if the object does not contain gamma thresholds 
#' or if the effective sample size cannot be computed.
essGammas <- function(object, ...) {
    UseMethod("essGammas")
}

#' The Geweke diagnostic for the beta coefficients in a multi-scale probit model object.
#'
#' @param object An mspm object.
#' @param ... Additional arguments (not used).
#' @return A list containing the Geweke diagnostic results for each beta coefficient. Each element of 
#' the list corresponds to a beta coefficient and contains a list with the following components:
#' \itemize{
#'   \item z: The Geweke z-score for the beta coefficient.
#'   \item p: The p-value associated with the Geweke z-score.
#' } Returns NULL if the object does not contain beta coefficients or if the Geweke diagnostic 
#' cannot be computed.
gewekeBeta <- function(object, ...) {
    UseMethod("gewekeBeta")
}

#' The Geweke diagnostic for the gamma thresholds in a multi-scale probit model object.
#'
#' @param object An mspm object.
#' @param ... Additional arguments (not used).
#' @return A list containing the Geweke diagnostic results for each gamma threshold. Each element of 
#' the list corresponds to a target dataset and contains a list of Geweke diagnostic results for 
#' the gamma thresholds of that target. Each Geweke diagnostic result is a list with the following 
#' components:
#' \itemize{
#'   \item z: The Geweke z-score for the gamma threshold.
#'   \item p: The p-value associated with the Geweke z-score.
#' } Returns NULL if the object does not contain gamma thresholds or if the Geweke diagnostic 
#' cannot be computed.
gewekeGammas <- function(object, ...) {
    UseMethod("gewekeGammas")
}

gelmanRhatBeta <- function(object, ...) {
    UseMethod("gelmanRhatBeta")
}

gelmanRhatGammas <- function(object, ...) {
    UseMethod("gelmanRhatGammas")
}

get_n_temperatures <- function(object, ...) {
    UseMethod("get_n_temperatures")
}

get_inv_temperatures <- function(object, ...) {
    UseMethod("get_inv_temperatures")
}