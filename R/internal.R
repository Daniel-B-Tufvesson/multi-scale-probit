source("R/generics.R")

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

predictorNames.mspm_data_spec <- function(object, ...) {
    object$predictorNames
}

responseNames.mspm_data_spec <- function(object, ...) {
    object$responseNames
}

levelNames.mspm_data_spec <- function(object, ...) {
    object$levelNames
}

nlevels.mspm_data_spec <- function(object, ...) {
    object$nlevels
}

ntargets.mspm_data_spec <- function(object, ...) {
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

data_spec.mspm_data <- function(object, ...) {
    object$data_spec
}

ntargets.mspm_data <- function(object, ...) {
    ntargets(object$data_spec)
}

nlevels.mspm_data <- function(object, ...) {
    nlevels(object$data_spec)
}

predictorNames.mspm_data <- function(object, ...) {
    predictorNames(object$data_spec)
}

responseNames.mspm_data <- function(object, ...) {
    responseNames(object$data_spec)
}

levelNames.mspm_data <- function(object, ...) {
    levelNames(object$data_spec)
}

#' Constructor for creating a new multi-scale probit model single chain diagnostics object.
#'
#' @param essBeta A numeric vector containing the effective sample size (ESS) for each beta 
#' coefficient.
#' @param essGammas A list of numeric vectors containing the effective sample size (ESS) for each 
#' gamma threshold in each target dataset.
#' @param gewekeBeta A numeric vector containing the Geweke diagnostic z-scores for each beta 
#' coefficient.
#' @param gewekeGammas A list of numeric vectors containing the Geweke diagnostic z-scores for 
#' each gamma threshold in each target dataset.
new_mspm_single_chain_diag <- function (
    essBeta,
    essGammas,
    gewekeBeta,
    gewekeGammas
) {
    structure(
        list(
            essBeta = essBeta,
            essGammas = essGammas,
            gewekeBeta = gewekeBeta,
            gewekeGammas = gewekeGammas
        ),
        class = "mspm_single_chain_diag"
    )
}

essBeta.mspm_single_chain_diag <- function(object, ...) {
    object$essBeta
}

essGammas.mspm_single_chain_diag <- function(object, ...) {
    object$essGammas
}

gewekeBeta.mspm_single_chain_diag <- function(object, ...) {
    object$gewekeBeta
}

gewekeGammas.mspm_single_chain_diag <- function(object, ...) {
    object$gewekeGammas
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
# diagnostics: An mspm_single_chain_diag object containing diagnostic information about the 
# MCMC sampling, such as effective sample size (ess) and Geweke diagnostics for each parameter.
# samplingTime: The time taken for the sampling phase of MCMC sampling.
# burninTime: The time taken for the burn-in phase of MCMC sampling.
# call: The original function call used to create the model.
new_mspm <- function(
    data_spec,
    beta,
    gammas,
    meanPrior,
    precPrior,
    adaptTune,
    tune,
    acceptanceRate,
    burninAcceptanceRate,
    seed,
    ndraws,
    ndrawsNoThin,
    thin,
    burnin,
    burninBeta,
    burninGammas,
    diagnostics,
    samplingTime,
    burninTime,
    call
) {
    structure(
        list(
            data_spec = data_spec,
            beta = beta,
            gammas = gammas,
            meanPrior = meanPrior,
            precPrior = precPrior,
            adaptTune = adaptTune,
            tune = tune,
            acceptanceRate = acceptanceRate,
            burninAcceptanceRate = burninAcceptanceRate,
            call = call,
            seed = seed,
            nobs = nrow(beta), # Why is nobs the number of rows in beta? Shouldn't this be ndraws?
            ndraws = ndraws,
            ndrawsNoThin = ndrawsNoThin,
            thin = thin,
            burnin = burnin,
            burninBeta = burninBeta, # May be null.
            burninGammas = burninGammas, # May be null.
            diagnostics = diagnostics, # May be null.
            samplingTime = samplingTime,
            burninTime = burninTime
        ),
        class = "mspm"
    )
}

data_spec.mspm <- function(object, ...) {
    object$data_spec
}

ntargets.mspm <- function(object, ...) {
    ntargets(object$data_spec)
}

nlevels.mspm <- function(object, ...) {
    nlevels(object$data_spec)
}

levelNames.mspm <- function(object, ...) {
    levelNames(object$data_spec)
}

predictorNames.mspm <- function(object, ...) {
    predictorNames(object$data_spec)
}

responseNames.mspm <- function(object, ...) {
    responseNames(object$data_spec)
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

ndraws.mspm <- function(object, ...) {
    object$ndraws
}

ndrawsNoThin.mspm <- function(object, ...) {
    object$ndrawsNoThin
}

burnin.mspm <- function(object, ...) {
    object$burnin
}

burninBeta.mspm <- function(object, ...) {
    object$burninBeta
}

burninGammas.mspm <- function(object, ...) {
    object$burninGammas
}

thin.mspm <- function(object, ...) {
    object$thin
}

samplingTime.mspm <- function(object, ...) {
    object$samplingTime
}

burninTime.mspm <- function(object, ...) {
    object$burninTime
}

diagnostics.mspm <- function(object, ...) {
    object$diagnostics
}

essBeta.mspm <- function(object, ...) {
    essBeta(diagnostics(object))
}

essGammas.mspm <- function(object, ...) {
    essGammas(diagnostics(object))
}

gewekeBeta.mspm <- function(object, ...) {
    gewekeBeta(diagnostics(object))
}

gewekeGammas.mspm <- function(object, ...) {
    gewekeGammas(diagnostics(object))
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
#' @param burnin The number of burn-in iterations used in MCMC sampling.
#' @param burninBeta An MCMC object containing the burn-in samples for the beta parameters.
#' @param burninGammas A list of MCMC objects containing the burn-in samples for the gamma 
#' parameters of each target dataset.
#' @param diagnostics An mspm_single_chain_diag object containing diagnostic information about the 
#' MCMC sampling, such as effective sample size (ess) and Geweke diagnostics for each parameter.
#' @param samplingTime The time taken for the sampling phase of MCMC sampling.
#' @param burninTime The time taken for the burn-in phase of MCMC sampling
#' @param call The original function call used to create the model.
new_mspm_pt <- function(
    data_spec,
    beta,
    gammas,
    meanPrior,
    precPrior,
    seed,
    ndraws,
    ndrawsNoThin,
    thin,
    ntemperatures,
    burnin,
    burninBeta,
    burninGammas,
    diagnostics,
    samplingTime,
    burninTime,
    initialTemperatureLadder,
    adjustedTemperatureLadder,
    targetSwapAcceptRatio,
    actualSwapAcceptRatio,
    swapRatios,
    samplingSwapAcceptRatio,
    nacceptedSwaps,
    nproposedSwaps,
    ladderLearningRate,
    initialWindowSize,
    windowGrowthFactor,
    completeSwapping,
    call
) {
    structure(
        list(
            data_spec = data_spec,
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
            burninBeta = burninBeta, # May be null.
            burninGammas = burninGammas, # May be null.

            ntemperatures = ntemperatures,
            initialTemperatureLadder = initialTemperatureLadder,
            adjustedTemperatureLadder = adjustedTemperatureLadder,
            targetSwapAcceptRatio = targetSwapAcceptRatio,
            actualSwapAcceptRatio = actualSwapAcceptRatio,
            swapRatios = swapRatios,
            samplingSwapAcceptRatio = samplingSwapAcceptRatio,
            nacceptedSwaps = nacceptedSwaps,
            nproposedSwaps = nproposedSwaps,
            ladderLearningRate = ladderLearningRate,
            initialWindowSize = initialWindowSize,
            windowGrowthFactor = windowGrowthFactor,
            
            completeSwapping = completeSwapping,

            samplingTime = samplingTime,
            burninTime = burninTime,
            diagnostics = diagnostics
        ),
        class = c("mspm_pt", "mspm")
    )
}

ntemperatures.mspm_pt <- function(object, ...) {
    object$ntemperatures
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

data_spec.mspm_labeled_prediction <- function(object, ...) {
    object$data_spec
}

ntargets.mspm_labeled_prediction <- function(object, ...) {
    ntargets(object$data_spec)
}

nlevels.mspm_labeled_prediction <- function(object, ...) {
    nlevels((object$data_spec))
}

predictorNames.mspm_labeled_prediction <- function(object, ...) {
    predictorNames(object$data_spec)
}

responseNames.mspm_labeled_prediction <- function(object, ...) {
    responseNames(object$data_spec)
}

levelNames.mspm_labeled_prediction <- function(object, ...) {
    levelNames(object$data_spec)
}

ndraws.mspm_labeled_prediction <- function(object, ...) {
    object$ndraws
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

data_spec.mspm_labeled_evaluation <- function(object, ...) {
    object$data_spec
}

ntargets.mspm_labeled_evaluation <- function(object, ...) {
    ntargets(object$data_spec)
}

nlevels.mspm_labeled_evaluation <- function(object, ...) {
    nlevels(object$data_spec)
}

predictorNames.mspm_labeled_evaluation <- function(object, ...) {
    predictorNames(object$data_spec)
}

responseNames.mspm_labeled_evaluation <- function(object, ...) {
    responseNames(object$data_spec)
}

levelNames.mspm_labeled_evaluation <- function(object, ...) {
    levelNames(object$data_spec)
}

ndraws.mspm_labeled_evaluation <- function(object, ...) {
    object$ndraws
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

#' Constructor for creating a new multi-scale probit model cross-validation result object.
#'
#' @param data_spec The data specification object containing information about the predictors, 
#' responses, levels, etc.
#' @param nsplits The number of cross-validation splits performed.
#' @param allEvaluations A list containing the evaluation results for each split. Each element 
#' is an mspm_labeled_evaluation object.
#' @param means A list containing the mean performance scores over the draws for each split. Each 
#' element is a matrix where each row is a split and each column is a metric. 
#' @param meansOnly A boolean indicating whether only the mean performance scores are included in 
#' the result (TRUE) or if the full evaluation results for each split are included (FALSE).
#' @param seed The random seed used for reproducibility.
#' @param gelmanRhatBeta A numeric vector containing the Gelman-Rubin R-hat diagnostic for the beta
#' parameters across all splits.
#' @param gelmanRhatGammas A list of numeric vectors containing the Gelman-Rubin R-hat diagnostic for
#' the gamma parameters for each target dataset across all splits.
#' @param allBetas A list of MCMC objects containing the posterior draws for the beta parameters for
#' each split. This may be NULL.
#' @param allGammas A list of lists of MCMC objects containing the posterior draws for the gamma
#' parameters for each target dataset for each split. This may be NULL.
#' @param allBurninBetas A list of MCMC objects containing the burn-in draws for the beta parameters 
#' for each split. This may be NULL.
#' @param allBurninGammas A list of lists of MCMC objects containing the burn-in draws for the gamma 
#' parameters for each target dataset for each split. This may be NULL.
#' @param allTuneParams A list containing the proposal variances used for each split.
#' @param allBurninTimes A vector containing the burn-in times for each split. Given in seconds.
#' @param allSamplingTimes A vector containing the sampling times for each split. Given in seconds.
new_mspm_cv_result <- function(
    data_spec,
    nsplits,
    metrics,
    allEvaluations,
    means,
    meansOnly,
    seed,
    gelmanRhatBeta,
    gelmanRhatGammas,
    allBetas,
    allGammas,
    allBurninBetas,
    allBurninGammas,
    allTuneParams,
    allBurninTimes = allBurninTimes,
    allSamplingTimes = allSamplingTimes,
    call
) {
    structure(
        list(
            data_spec = data_spec,
            nsplits = nsplits,
            metrics = metrics,
            allEvaluations = allEvaluations,
            means = means,
            meansOnly = meansOnly,
            seed = seed,
            gelmanRhatBeta = gelmanRhatBeta,
            gelmanRhatGammas = gelmanRhatGammas,
            allBetas = allBetas,
            allGammas = allGammas,
            allBurninBetas = allBurninBetas,
            allBurninGammas = allBurninGammas,
            allTuneParams = allTuneParams,
            allBurninTimes = allBurninTimes,
            samplingTime = allSamplingTimes,
            call = call
        ),
        class = "mspm_cv_result"
    )
}

data_spec.mspm_cv_result <- function(object, ...) {
    object$data_spec
}

ntargets.mspm_cv_result <- function(object, ...) {
    ntargets(object$data_spec)
}

evalMetrics.mspm_cv_result <- function(object, ...) {
    object$metrics
}

nsplits.mspm_cv_result <- function(object, ...) {
    object$nsplits
}

cvAllEvaluations.mspm_cv_result <- function(object, ...) {
    object$allEvaluations
}

cvMeans.mspm_cv_result <- function(object, ...) {
    object$means
}

cvAllDraws.mspm_cv_result <- function(object, ...) {

    allEvals <- cvAllEvaluations(object)
    if (is.null(allEvals)) {
        return(NULL)
    }

    # Get metrics and number of splits.
    metrics <- evalMetrics(object)
    nsplits <- nsplits(object)

    # Aggregate all evaluation draws for each metric across all splits.
    results <- list()
    for (metric in metrics) {
        results[[metric]] <- do.call(rbind, lapply(allEvals, function(eval) {
            evalDrawResults(eval)[[metric]]
        }))
    }

    results
}

gelmanRhatBeta.mspm_cv_result <- function(object, ...) {
    object$gelmanRhatBeta
}

gelmanRhatGammas.mspm_cv_result <- function(object, ...) {
    object$gelmanRhatGammas
}