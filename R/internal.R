
#' Constructor for creating a data specification object for model data.
new_msmp_data_spec <- function(
    predictorNames,
    responseNames,
    levelNames,
    nlevels,
    ntargets,
    call
) {
    structure(
        list(
            predictorNames = predictorNames,
            responseNames = responseNames,
            levelNames = levelNames,
            nlevels = nlevels,
            ntargets = ntargets,
            call = call
        ),
        class = "mspm_data_spec"
    )
}

get_predictor_names.mspm_data_spec <- function(object, ...) {
    object$predictorNames
}

get_n_predictors.mspm_data_spec <- function(object, ...) {
    length(get_predictor_names(object))
}

get_response_names.mspm_data_spec <- function(object, ...) {
    object$responseNames
}

get_level_names.mspm_data_spec <- function(object, ...) {
    object$levelNames
}

get_n_levels.mspm_data_spec <- function(object, ...) {
    object$nlevels
}

get_n_targets.mspm_data_spec <- function(object, ...) {
    object$ntargets
}


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
            Xlist = Xlist,
            ylist = ylist,
            call = call,

            data_spec = new_msmp_data_spec(
                predictorNames = predictorNames,
                responseNames = responseNames,
                levelNames = levelNames,
                nlevels = nlevels,
                ntargets = ntargets,
                call = call
            )
        ),
        class = "mspm_data"
    )
}

get_data_spec.mspm_data <- function(object, ...) {
    object$data_spec
}

get_x_values.mspm_data <- function(object, ...) {
    object$Xlist
}

get_y_values.mspm_data <- function(object, ...) {
    object$ylist
}

get_n_targets.mspm_data <- function(object, ...) {
    get_n_targets(object$data_spec)
}

get_n_levels.mspm_data <- function(object, ...) {
    get_n_levels(object$data_spec)
}

get_predictor_names.mspm_data <- function(object, ...) {
    get_predictor_names(object$data_spec)
}

get_n_predictors.mspm_data <- function(object, ...) {
    length(get_predictor_names(object))
}

get_response_names.mspm_data <- function(object, ...) {
    get_response_names(object$data_spec)
}

get_level_names.mspm_data <- function(object, ...) {
    get_level_names(object$data_spec)
}


#' Constructor for creating a new multi-scale probit model data split object for 
#' train/test splits.
new_mspm_data_split <- function(
    train_data,
    test_data,
    call
) {
    structure(
        list(
            train_data = train_data,
            test_data = test_data,
            call = call
        ),
        class = "mspm_data_split"
    )
}

get_train_split.mspm_data_split <- function(object, ...) {
    object$train_data
}

get_test_split.mspm_data_split <- function(object, ...) {
    object$test_data
}

get_data_spec.mspm_data_split <- function(object, ...) {
    get_data_spec(get_train_split(object))
}


# Constructor for creating a new multi-scale probit model (MSPM).
# 
# Arguments: 
# data_spec: The data specification object containing information about the predictors, responses, levels, etc.
# beta: An MCMC object for the beta (the coefficients) parameters.
# gammas: A list of MCMC gammas (the thresholds) of each target dataset.
# meanPrior: Prior mean for regression coefficients.
# precPrior: Prior precision for regression coefficients.
# adaptTune: The tuning parameter used for adaptive tuning of the sampler.
# tune: The final tuning parameter used for the sampler after adaptive tuning.
# acceptanceRate: The acceptance rate of the sampler for each target.
# burninAcceptanceRate: The acceptance rate during burn-in for each target.
# seed: The random seed used for reproducibility.
# ndraws: The number of posterior draws collected after thinning.
# ndrawsNoThin: The number of posterior draws collected before thinning.
# thin: The thinning interval used in MCMC sampling.
# burnin: The number of burn-in iterations used in MCMC sampling.
# burninBeta: An MCMC object containing the burn-in samples for the beta parameters.
# burninGammas: A list of MCMC objects containing the burn-in samples for the gamma parameters of 
# each target dataset.
# samplingTime: The time taken for the sampling phase of MCMC sampling.
# burninTime: The time taken for the burn-in phase of MCMC sampling.
# call: The original function call used to create the model.
new_mspm <- function(
    data_spec,
    beta,
    gammas,
    meanPrior,
    precPrior,
    proposal_variance,
    acceptanceRate,
    seed,
    ndraws,
    ndrawsNoThin,
    thin,
    burnin,
    samplingTime,
    burninTime,
    nlikelihood_calls,
    call
) {
    structure(
        list(
            data_spec = data_spec,
            beta = beta,
            gammas = gammas,
            meanPrior = meanPrior,
            precPrior = precPrior,
            proposal_variance = proposal_variance,
            acceptanceRate = acceptanceRate,
            call = call,
            seed = seed,
            ndraws = ndraws,
            ndrawsNoThin = ndrawsNoThin,
            thin = thin,
            burnin = burnin,
            samplingTime = samplingTime,
            burninTime = burninTime,
            nlikelihood_calls = nlikelihood_calls
        ),
        class = "mspm"
    )
}

get_data_spec.mspm <- function(object, ...) {
    object$data_spec
}

get_n_targets.mspm <- function(object, ...) {
    get_n_targets(object$data_spec)
}

get_n_levels.mspm <- function(object, ...) {
    get_n_levels(object$data_spec)
}

get_level_names.mspm <- function(object, ...) {
    get_level_names(object$data_spec)
}

get_predictor_names.mspm <- function(object, ...) {
    get_predictor_names(object$data_spec)
}

get_n_predictors.mspm <- function(object, ...) {
    length(get_predictor_names(object))
}

get_response_names.mspm <- function(object, ...) {
    get_response_names(object$data_spec)
}

get_beta.mspm <- function(object, ...) {
    object$beta
}

get_gammas.mspm <- function(object, ...) {
    object$gammas
}

get_mean_prior.mspm <- function(object, ...) {
    object$meanPrior
}

get_prec_prior.mspm <- function(object, ...) {
    object$precPrior
}

get_proposal_variance.mspm <- function(object, ...) {
    object$proposal_variance
}

get_n_likelihood_calls.mspm <- function(object, ...) {
    object$nlikelihood_calls
}

get_n_draws.mspm <- function(object, ...) {
    object$ndraws
}

get_n_draws_no_thin.mspm <- function(object, ...) {
    object$ndrawsNoThin
}

get_burnin.mspm <- function(object, ...) {
    object$burnin
}

get_thin.mspm <- function(object, ...) {
    object$thin
}

get_sampling_time.mspm <- function(object, ...) {
    object$samplingTime
}

get_burnin_time.mspm <- function(object, ...) {
    object$burninTime
}

#' Constructor for creating a new multi-scale probit model fitted object for parallel tempering.
#' This inherits from the mspm class and includes additional information specific to parallel 
#' tempering.
#' @param data_spec The data specification object containing information about the predictors, 
#' responses, levels, etc.
#' @param beta An MCMC object for the beta (the coefficients) parameters.
#' @param gammas A list of MCMC gammas (the thresholds) of each target dataset.
#' @param meanPrior Prior mean for regression coefficients.
#' @param precPrior Prior precision for regression coefficients.
#' @param seed The random seed used for reproducibility.
#' @param ndraws The number of posterior draws collected after thinning.
#' @param ndrawsNoThin The number of posterior draws collected before thinning.
#' @param thin The thinning interval used in MCMC sampling.
#' @param inv_temperatures The inverse temperature ladder used in parallel tempering.
#' @param burnin The number of burn-in iterations used in MCMC sampling.
#' @param burninBeta An MCMC object containing the burn-in samples for the beta parameters.
#' @param burninGammas A list of MCMC objects containing the burn-in samples for the gamma 
#' parameters of each target dataset.
#' @param samplingTime The time taken for the sampling phase of MCMC sampling.
#' @param burninTime The time taken for the burn-in phase of MCMC sampling
#' @param completeSwapping A boolean indicating whether complete swapping was used in parallel
#' tempering.
#' @param nlikelihood_calls The total number of likelihood calls made during MCMC sampling.
#' @param round_trip_times A numeric vector containing the round trip times for each completed 
#' round trip.
#' @param call The original function call used to create the model.
new_mspm_pt <- function(
    data_spec,
    beta,
    gammas,
    mean_prior,
    prec_prior,
    proposal_variance,
    seed,
    ndraws,
    ndrawsNoThin,
    thin,
    inv_temperatures,
    burnin,
    samplingTime,
    burninTime,
    completeSwapping,
    nlikelihood_calls,
    round_trip_times,
    call
) {
    structure(
        list(
            data_spec = data_spec,
            beta = beta,
            gammas = gammas,
            mean_prior = mean_prior,
            prec_prior = prec_prior,
            proposal_variance = proposal_variance,
            call = call,
            seed = seed,
            ndraws = ndraws,
            ndrawsNoThin = ndrawsNoThin,
            thin = thin,
            burnin = burnin,

            inv_temperatures = inv_temperatures,
            ntemperatures = length(inv_temperatures),
            completeSwapping = completeSwapping,

            samplingTime = samplingTime,
            burninTime = burninTime,
            nlikelihood_calls = nlikelihood_calls,
            round_trip_times = round_trip_times
        ),
        class = c("mspm_pt", "mspm")
    )
}

get_n_temperatures.mspm_pt <- function(object, ...) {
    object$ntemperatures
}

get_proposal_variance.mspm_pt <- function(object, ...) {
    object$proposal_variance
}

get_inv_temperatures.mspm_pt <- function(object, ...) {
    object$inv_temperatures
}


# Constructor for creating a new multi-scale probit model labeled prediction object.
# 
# Arguments:
# data_spec: The data specification object containing information about the predictors, responses, levels, etc.
# ndraws: The number of posterior draws used for prediction.
# ystars: A list of matrices of latent variable samples for each dataset. Each matrix has rows 
#         corresponding to observations and columns to posterior draws.
# ylabels: A list of data frames of observed categorical labels for each observation and 
#          each draw. There is one data frame for each scale.
# ylabelIndexes: A list of matrices indicating the indexes of predicted labels for each
#                observation and each draw. There is one matrix for each scale.
# call: The original function call used to create the prediction.
new_mspm_labeled_prediction <- function(
    data_spec,
    ndraws,
    ystars,
    ylabels,
    ylabelIndexes,
    call
) {
    structure(
        list(
            data_spec = data_spec,
            ndraws = ndraws,
            ystars = ystars,
            ylabels = ylabels,
            ylabelIndexes = ylabelIndexes,
            call = call
        ),
        class = "mspm_labeled_prediction"
    )
}

get_data_spec.mspm_labeled_prediction <- function(object, ...) {
    object$data_spec
}

get_n_targets.mspm_labeled_prediction <- function(object, ...) {
    get_n_targets(object$data_spec)
}

get_n_levels.mspm_labeled_prediction <- function(object, ...) {
    get_n_levels((object$data_spec))
}

get_predictor_names.mspm_labeled_prediction <- function(object, ...) {
    get_predictor_names(object$data_spec)
}

get_response_names.mspm_labeled_prediction <- function(object, ...) {
    get_response_names(object$data_spec)
}

get_level_names.mspm_labeled_prediction <- function(object, ...) {
    get_level_names(object$data_spec)
}

get_n_draws.mspm_labeled_prediction <- function(object, ...) {
    object$ndraws
}

get_predicted_labels.mspm_labeled_prediction <- function(object, ...) {
    object$ylabels
}

get_predicted_label_indexes.mspm_labeled_prediction <- function(object, ...) {
    object$ylabelIndexes
}

get_latent_prediction.mspm_labeled_prediction <- function(object, ...) {
    object$ystars
}

# Constructor for creating a new multi-scale probit model labeled evaluation object.
# 
# Arguments:
# data_spec: The data specification object containing information about the predictors, responses, levels, etc.
# ndraws: The number of posterior draws used for evaluation.
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
    data_spec,
    ndraws,
    metrics,
    drawResults,
    targetMeans,
    drawMeans,
    metricMeans,
    call
) {
    structure(
        list(
            data_spec = data_spec,
            ndraws = ndraws,
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

get_data_spec.mspm_labeled_evaluation <- function(object, ...) {
    object$data_spec
}

get_n_targets.mspm_labeled_evaluation <- function(object, ...) {
    get_n_targets(object$data_spec)
}

get_n_levels.mspm_labeled_evaluation <- function(object, ...) {
    get_n_levels(object$data_spec)
}

get_predictor_names.mspm_labeled_evaluation <- function(object, ...) {
    get_predictor_names(object$data_spec)
}

get_response_names.mspm_labeled_evaluation <- function(object, ...) {
    get_response_names(object$data_spec)
}

get_level_names.mspm_labeled_evaluation <- function(object, ...) {
    get_level_names(object$data_spec)
}

get_n_draws.mspm_labeled_evaluation <- function(object, ...) {
    object$ndraws
}

get_eval_metrics.mspm_labeled_evaluation <- function(object, ...) {
    object$metrics
}

get_eval_draw_results.mspm_labeled_evaluation <- function(object, ...) {
    object$drawResults
}

get_eval_target_means.mspm_labeled_evaluation <- function(object, ...) {
    object$targetMeans
}

get_eval_draw_means.mspm_labeled_evaluation <- function(object, ...) {
    object$drawMeans
}

get_eval_metric_means.mspm_labeled_evaluation <- function(object, ...) {
    object$metricMeans
}


#' Constructor for creating a new multi-scale probit model tuning results object.
#'
#' @param data_spec The data specification object containing information about the predictors,
#' responses, levels, etc.
#' @param proposal_variance The final approximation of the tuning parameter.
#' @param acceptance_rates A vector of the final acceptance rates for each target dataset.
#' @param target_acceptance_rate The target acceptance rate used for tuning.
#' @param max_iterations The maximum number of iterations used for tuning.
#' @param final_iteration The final iteration number reached during tuning.
#' @param seed The random seed used for reproducibility.
#' @param call The original function call used to create the tuning results object.
new_mspm_tune_results <- function(
    data_spec,
    proposal_variance,
    acceptance_rates,
    target_acceptance_rate,
    max_iterations,
    final_iteration,
    seed,
    call
) {
    structure(
        list(
            data_spec = data_spec,
            proposal_variance = proposal_variance,
            acceptance_rates = acceptance_rates,
            target_acceptance_rate = target_acceptance_rate,
            max_iterations = max_iterations,
            final_iteration = final_iteration,
            seed = seed,
            call = call
        ),
        class = "mspm_tune_results"
    )
}

get_proposal_variance.mspm_tune_results <- function(object, ...) {
    object$proposal_variance
}

new_mspm_tune_results_pt <- function(
    data_spec,
    proposal_variance,
    proposal_acceptance_rates,
    target_acceptance_rate,
    inv_temperature_ladder,
    target_temp_swap_rate,
    temp_swap_rates,
    max_iterations,
    final_iteration,
    seed,
    call
) {
    structure(
        list(
            data_spec = data_spec,

            proposal_variance = proposal_variance,
            proposal_acceptance_rates = proposal_acceptance_rates,
            target_acceptance_rate = target_acceptance_rate,

            inv_temperature_ladder = inv_temperature_ladder,
            target_temp_swap_rate = target_temp_swap_rate,
            temp_swap_rates = temp_swap_rates,

            max_iterations = max_iterations,
            final_iteration = final_iteration,
            seed = seed,
            call = call
        ),
        class = "mspm_tune_results_pt"
    )
}

get_proposal_variance.mspm_tune_results_pt <- function(object, ...) {
    object$proposal_variance
}

get_inv_temperatures.mspm_tune_results_pt <- function(object, ...) {
    object$inv_temperature_ladder
}