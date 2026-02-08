source("R/data.R")
source("R/fit.R")
source("R/eval.R")
source("R/util.R")

# Libraries for parallelization.
library(doParallel)
library(foreach)

#' Cross-validate an mspm model using the specified data. This is a Monte Carlo cross-validation,
#' also known as repeated random sub-sampling validation. 
#'
cross_validate <- function(
    data,
    nsplits,
    prop,
    ndraws,
    burnin,
    thin,
    ...,
    meanPrior = NULL,
    precPrior = NULL,
    seed = NULL,
    metrics = c("f1", "kendall"),
    nworkers = 1,
    meansOnly = TRUE
) {
    
    # 1 worker -> no parallelization, run sequentially.
    if (nworkers == 1) {
        res <- list()
        for (i in 1:nsplits) {
            res[[i]] <- .do_cv_trial(
                data = data,
                prop = prop,
                ndraws = ndraws,
                burnin = burnin,
                thin = thin,
                meanPrior = meanPrior,
                precPrior = precPrior,
                seed = seed + i,
                metrics = metrics,
                meansOnly = meansOnly
            )
        }
    }
    # >1 workers -> parallelize across splits.
    else if (nworkers > 1) {
        cl<-makeForkCluster(8)
        registerDoParallel(cl)
        clusterSetRNGStream(cl, 12345)

        # Export necessary environment variables to the workers.
        exportVars <- c(
            "data",
            "prop", 
            "ndraws", 
            "burnin", 
            "thin", 
            "meanPrior", 
            "precPrior", 
            "metrics", 
            "seed",
            "meansOnly",
            ".do_cv_trial"
        )

        # Run the cross-validation splits in parallel.
        res <- foreach(i = 1:nsplits, .combine = 'list', .export = exportVars) %dopar% {
            .do_cv_trial(
                data = data,
                prop = prop,
                ndraws = ndraws,
                burnin = burnin,
                thin = thin,
                meanPrior = meanPrior,
                precPrior = precPrior,
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

    # Return result.
    new_mspm_cv_result(
        data_spec = data_spec(data),
        nsplits = nsplits,
        allEvaluations = allEvals,
        means = allMeans,
        meansOnly = meansOnly,
        seed,
        call = match.call()
    )
}

#' Run a single trial for Monte Carlo cross-validation.
#'
.do_cv_trial <- function(
    data,
    prop,
    ndraws,
    burnin,
    thin,
    meanPrior,
    precPrior,
    seed,
    metrics,
    meansOnly
) {
    # Split data into train and test sets.
    splits <- split_data(data, prop = prop, seed = seed)
    train_data <- splits$train
    test_data <- splits$test

    # Fit model on training data.
    fit <- fit_mspm(
        data = train_data,
        ndraws = ndraws,
        burnin = burnin,
        thin = thin,
        meanPrior = meanPrior,
        precPrior = precPrior,
        seed = seed
    )

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
            # if (any(metricValues == 0, na.rm = TRUE)) {
            #     meansMatrix[1, nmetrics + 1] <- 0
            # } else {
            #     #meansMatrix[1, nmetrics + 1] <- length(metricValues) / sum(1/metricValues, na.rm = TRUE)
                
            # }
        }

        means[[i]] <- meansMatrix
    }

    # Return the results.
    if (meansOnly) {
        return (list(
            means = means
        ))
    }
    else {
        return(list(
            results = results,
            means = means
        ))
    }
}


