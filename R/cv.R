source("R/data.R")
source("R/fit.R")
source("R/eval.R")
source("R/util.R")

# Libraries for parallelization.
library(doParallel)
library(foreach)

library(coda)

#' Cross-validate an mspm model using the specified data. This is a Monte Carlo cross-validation,
#' also known as repeated random sub-sampling validation. 
#'
#' Note that this function can be computationally intensive, especially with a large number of 
#' splits and posterior draws.
#'
#' @param data An mspm_data object containing the dataset to be used for cross-validation.
#' @param nsplits An integer specifying the number of random splits to perform.
#' @param prop A numeric value between 0 and 1 indicating the proportion of data to be used for 
#' training in each split. The rest will be used for testing.
#' @param sampler A function that fits an mspm model and returns an object containing posterior draws.
#' @param samplerArgs A list of additional arguments to pass to the sampler function, such
#' @param ndraws An integer specifying the number of posterior draws to use when fitting the model.
#' @param metrics A character vector specifying the evaluation metrics to compute. Default 
#' is c("f1", "kendall").
#' @param nworkers An integer specifying the number of worker processes to use for parallelization.
#' If nworkers = 1, cross-validation will run sequentially. If nworkers > 1, it will run in 
#' parallel across splits. Default is 1.
#' @param meansOnly A logical value indicating whether to return only the mean evaluations for 
#' each target (TRUE) or to return the full evaluation objects for each split (FALSE).
#' @param computeDiagnostics A logical value indicating whether to compute diagnostics for each 
# fitted model.
#'
#' @return An object of class 'mspm_cv_result' containing the results of cross-validation, 
#' including mean evaluation metrics for each target and, optionally, the full evaluation objects 
#' for each split.
cross_validate <- function(
    data,
    nsplits,
    prop,
    sampler,
    samplerArgs,
    ndraws,
    ...,
    pretune_function = NULL,
    pretune_args = NULL,
    metrics = c("f1", "kendall"),
    nworkers = 1,
    meansOnly = TRUE,
    computeDiagnostics = TRUE
) {
    # Use random seed. Note, however, that we cannot guarantee perfect reproducibility, even for
    # when nworkers = 0. This is most likely caused by underlying parallelization of the math
    # libraries used by the sampler. We will investigate this issue further and try to find a 
    # solution. For cross-validation, however, we usually only care about average performance
    # rather than specific performance, so lack of seeding should not be negatively impact the 
    # utility of this function.
    seed <- sample.int(1e6, 1)
    
    # 1 worker -> no parallelization, run sequentially.
    if (nworkers == 1) {
        res <- list()
        for (i in 1:nsplits) {
            res[[i]] <- .do_cv_trial(
                data = data,
                prop = prop,
                ndraws = ndraws,
                sampler = sampler,
                samplerArgs = samplerArgs,
                pretune_function = pretune_function,
                pretune_args = pretune_args,
                seed = seed + i,
                metrics = metrics,
                meansOnly = meansOnly
            )
        }
    }
    # >1 workers -> parallelize across splits.
    else if (nworkers > 1) {
        cl<-makeForkCluster(nworkers)
        on.exit(stopCluster(cl))
        registerDoParallel(cl)
        clusterSetRNGStream(cl, 12345)

        # Export necessary environment variables to the workers.
        exportVars <- c(
            ".do_cv_trial"
        )
        # Combine results as a list.
        combine <- function(x, y) {
            is_res_x <- "means" %in% names(x)
            is_res_y <- "means" %in% names(y)
            if (is_res_x && is_res_y) return(c(list(x), list(y)))
            if (is_res_x) return(c(list(x), y))
            if (is_res_y) return(c(x, list(y)))
            return(c(x, y))
        }

        # Run the cross-validation splits in parallel.
        res <- foreach(i = 1:nsplits, .combine = combine, .export = exportVars) %dopar% {
            .do_cv_trial(
                data = data,
                prop = prop,
                ndraws = ndraws,
                sampler = sampler,
                samplerArgs = samplerArgs,
                pretune_function = pretune_function,
                pretune_args = pretune_args,
                seed = seed + i, # Use different seed for each split.
                metrics = metrics,
                meansOnly = meansOnly
            )
        }
    }
    else {
        stop("nworkers must be a positive integer.")
    }


    # Aggregate means as a matrix, with rows for each trial and cols for each metric.
    allMeans <- list()
    for (i in 1:ntargets(data)) {
        targetMeans <- matrix(ncol = length(metrics) + 1, nrow = nsplits)
        colnames(targetMeans) <- if (length(metrics) > 1) c(metrics, "HarmonicMean") else metrics
        for (j in 1:nsplits) {
            targetMeans[j, ] <- res[[j]]$means[[i]]
        }
        allMeans[[i]] <- targetMeans
    }
    
    # Aggregate all eval objects as a list.
    allEvals <- NULL
    if (!meansOnly) {
        allEvals <- list()
        for (i in 1:nsplits) {
            allEvals[[i]] <- res[[i]]$results
        }
    }

    # Compute diagnostics if requested.
    rhatBeta <- NULL
    rhatGammas <- NULL
    if (computeDiagnostics) {
        betaList <- list()
        gammasList <- list()

        # Collect MCMC chains.
        for (i in 1:nsplits) {
            betaList[[i]] <- res[[i]]$fit$beta
            gammasList[[i]] <- res[[i]]$fit$gammas
        }

        # Compute Gelman-Rubin R-hat diagnostics.
        rhatBeta <- gelman.diag(betaList)
        rhatGammas <- list()
        for (i in 1:ntargets(data)) {
            gammas_i <- lapply(gammasList, function(g) g[[i]])
            rhatGammas[[i]] <- gelman.diag(mcmc.list(gammas_i))
        }
    }

    # Aggregate tune results as list.
    all_tune_results <- NULL
    if (!is.null(pretune_function)) {
        all_tune_results <- list()
        for (i in 1:nsplits) {
            all_tune_results[[i]] <- res[[i]]$tune_results
        }
    }

    # Aggregate nlikelihood_calls as vector.
    all_nlikelihood_calls <- sapply(res, function(r) {
        if (!is.null(r$fit)) {
            return(get_nlikelihood_calls(r$fit))
        }
        return(NA)
    })

    # Return result.
    new_mspm_cv_result(
        data_spec = data_spec(data),
        nsplits = nsplits,
        metrics = metrics,
        allEvaluations = allEvals,
        means = allMeans,
        meansOnly = meansOnly,
        seed = seed,
        gelmanRhatBeta = rhatBeta,
        gelmanRhatGammas = rhatGammas,
        all_tune_results = all_tune_results,
        all_nlikelihood_calls = all_nlikelihood_calls,
        call = match.call()
    )
}

#' Run a single trial for Monte Carlo cross-validation.
#'
#' @param data An mspm_data object containing the dataset to be used for this trial.
#' @param prop A numeric value between 0 and 1 indicating the proportion of data to be used for 
#' training in this trial. The rest will be used for testing.
#' @param ndraws An integer specifying the number of posterior draws to use when fitting the model.
#' @param sampler A function that fits an mspm model and returns an object containing posterior draws.
#' @param samplerArgs A list of additional arguments to pass to the sampler function.
#' @param seed An integer random seed for reproducibility.
#' @param metrics A character vector specifying the evaluation metrics to compute.
#' @param meansOnly A logical value indicating whether to return only the mean evaluations for 
#' each target (TRUE) or to return the full evaluation objects for this trial (FALSE).
#'
#' @return A list containing either:
#' \item{means}{A list of matrices containing mean evaluation metrics for each target.}
#' \item{results}{An mspm_labeled_evaluation object containing the full evaluation results for 
#' this trial.}
.do_cv_trial <- function(
    data,
    prop,
    ndraws,
    sampler,
    samplerArgs,
    pretune_function,
    pretune_args,
    seed,
    metrics,
    meansOnly
) {
    # Split data into train and test sets.
    splits <- split_data(data, prop = prop, seed = seed)
    train_data <- splits$train
    test_data <- splits$test

    # Pre-tune the sampler.
    tune_results <- NULL
    if (!is.null(pretune_function)) {
        pretune_args <- c(list(data = train_data, seed = seed), pretune_args)
        tune_results <- do.call(pretune_function, pretune_args)

        # Ugly unpacking of the params. We may fix this later. 
        if (class(tune_results) == "mspm_tune_results") {
            samplerArgs <- c(samplerArgs, tune_results$final_tune)
        }
        else if (class(tune_results) == "mspm_tune_results_pt") {
            samplerArgs <- c(samplerArgs, tune_results$final_tune)
            # Todo: extract the tuned temperature ladder as well.
            stop("Pre-tuning for parallel tempering is not yet implemented in the cross-validation function.")
        }
    }

    # Fit model on training data.
    args <- list(
        data = train_data,
        ndraws = ndraws,
        seed = seed
    )
    args <- c(args, samplerArgs)
    fit <- do.call(sampler, args)

    # Predict on test data.
    predictions <- predict_mspm(fit, newdata = test_data)

    # Evaluate predictions using specified metrics.
    results <- eval_mspm_prediction_draws(predictions, test_data, metrics = metrics)

    # Extract means into matrix.
    means <- list()
    nmetrics <- length(metrics)
    ncols <- if (nmetrics > 1) nmetrics + 1 else nmetrics
    colnames <- if (nmetrics > 1) c(metrics, "HarmonicMean") else metrics
    for (i in 1:ntargets(data)) {
        targetMeans <- evalTargetMeans(results)
        meansMatrix <- matrix(ncol = ncols, nrow = 1)
        colnames(meansMatrix) <- colnames
        
        for (j in 1:nmetrics) {
            meansMatrix[1, j] <- targetMeans[[metrics[j]]][i]
        }

        # Compute harmonic mean across metrics for this target.
        if (nmetrics > 1) {
            metricValues <- meansMatrix[1, 1:nmetrics]
            meansMatrix[1, nmetrics + 1] <- harmonic_mean(metricValues)
        }

        means[[i]] <- meansMatrix
    }

    # Return the results.
    if (meansOnly) {
        return (list(
            means = means,
            fit = fit
        ))
    }
    else {
        return(list(
            results = results,
            means = means,
            fit = fit,
            tune_results = tune_results
        ))
    }
}


