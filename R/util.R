source("R/internal.R")

#' Compute the harmonic mean of a numeric vector, handling zeros and NAs appropriately. If any
#' value is zero, the harmonic mean is defined to be zero. If any value is NA, the harmonic mean 
#' is NA.
#'
#' @param x A numeric vector for which to compute the harmonic mean.
#' @return The harmonic mean of the input vector, or 0 if any value is zero, or NA if any value 
#' is NA.
harmonic_mean <- function(x) {
    if (any(x == 0, na.rm = TRUE)) {
        return(0)
    } 
    else if (any(is.na(x))) {
        return(NA)
    }
    else {
        return(length(x) / sum(1/x, na.rm = TRUE))
    }
}

# Simulate data from a probit model with given X, beta, and gamma thresholds.
# Returns FALSE if unable to generate valid categories within the limit of attempts,  
# or a list containing the generated data otherwise.
# 
# Arguments:
# X     : Design matrix (n x p). Note: n >= c for c categories.
# beta  : Coefficient vector (p x 1)
# gamma : Thresholds for categories (vector of length c-1 for c categories)
# limit : Maximum number of attempts to generate valid categories.
roprobit <- function(X, beta, gamma = 0, limit = 10) {

    # Check if X has enough rows for the number of categories.
    if (nrow(X) < (length(gamma) + 1)) {
        stop("Error in roprobit: Number of observations (rows of X) must be at least equal to number of categories.")
    }

    # Generate latent variable and observed categories until valid or limit reached.
    Y <- rep(0, nrow(X))
    counter <- 0
    while (length(table(Y)) != length(gamma)+1) {
            y.star <- X %*% beta + rnorm(nrow(X))
            Y <- findInterval(y.star, gamma)
            counter = counter + 1
            if (counter > limit) {
                return(FALSE)
            }
    }
    return(list(X = X, beta = beta, gamma = gamma,
                y.star = y.star, Y = Y))
}


# Compute the minimum value for each row in a data frame or matrix.
# This is useful for quickly extracting row-wise minimums from tabular data.
# Returns a vector of minimum values.
#
# Arguments:
# df : A data frame or matrix.
rowMin <- function(df) {
    apply(df, 1, min)
}


# Evaluate the predictive performance of a simulated probit model. 
# Returns a list of performance metrics (F1 score, Kendall's tau, harmonic mean).
# 
# Arguments:
# sim      : A list containing simulation results with beta coefficients, gamma thresholds, 
#            design matrices (X), and observed outcomes (y) for each dataset. This is the object returned
#            by the hprobit_gibbs function.
# test_set : An optional list of test datasets with design matrices (X) and observed outcomes (Y).
eval_sim_small <- function(sim, test_set = NULL) {
    # Get the number of posterior draws and datasets
    ndraws <- nrow(sim$beta)
    ndatasets <- length(sim$gammas)
    # Determine the number of result columns (for multi-scale, add one)
    if(ndatasets>1) {
        nresults <- ndatasets+1
    } else {
        nresults <- ndatasets
    }

    # Initialize matrices to store metrics for each draw and dataset
    F1.train <- matrix(0, nrow = ndraws, ncol = nresults)
    F1.test <- matrix(0, nrow = ndraws, ncol = nresults)
    kendall.train <- matrix(0, nrow = ndraws, ncol = nresults)
    kendall.test <- matrix(0, nrow = ndraws, ncol = nresults)
    harm.train <- matrix(0, nrow = ndraws, ncol = nresults)
    harm.test <- matrix(0, nrow = ndraws, ncol = nresults)

    # Loop over each dataset (group)
    for (i in 1:ndatasets) {
        # Compute the linear predictor for all draws
        y <- tcrossprod(sim$X[[i]], sim$beta)
        # For each draw, discretize the latent variable into categories using the sampled gammas
        y.train <- sapply(1:ndraws, function(j) {
            return(findInterval(y[, j], sim$gammas[[i]][j, ]))
        })

        # Reference labels (shifted to zero-based)
        ref <- sim$y[[i]]-1
        ref.fac <- factor(ref)
        # Compute F1 score for each draw
        F1.train[, i] <- cpp_fmeasure_distribution(y.train, ref, length(levels(ref.fac)))
        # Compute Kendall's tau for each draw
        kendall.train[, i] <- apply(y.train, 2, function(draw) {
            # return(fastkendall(draw, y = ref)$kendall.tau)
            return(suppressWarnings(cor(draw, ref, method = "kendall")))
        })

        # If a test set is provided, evaluate on test data
        if (!is.null(test_set)) {
        # Compute the linear predictor for test data
        y <- tcrossprod(test_set[[i]]$X, sim$beta)
        # Discretize using sampled gammas for each draw
        y.test <- sapply(1:ndraws, function(j) {
            return(findInterval(y[, j], sim$gammas[[i]][j, ]))
        })
        # Reference test labels
        ref <- test_set[[i]]$Y
        ref.fac <- factor(ref)
        # Compute F1 score for test data
        F1.test[, i] <- cpp_fmeasure_distribution(y.test, ref, length(levels(ref.fac)))
        # Compute Kendall's tau for test data
        kendall.test[, i] <- apply(y.test, 2, function(draw) {
            return(suppressWarnings(cor(draw, ref, method = "kendall")))
        })

        # Identify valid draws (Kendall's tau >= 0 and not NA)
        ok <- kendall.test[, i] >= 0
        ok[which(is.na(ok))] <- FALSE
        # Compute harmonic mean of F1 and Kendall's tau for valid draws
        harm.test[ok, i] <- cpp_harmonic_rowmeans(as.matrix(data.frame(F1 = F1.test[ok, i],
                                                                        kendall = kendall.test[ok, i])))
        }
    }
    # (Optional) Code for computing harmonic means across datasets is commented out
    # Return all computed metrics and posterior samples
    return(list(
        F1.train = F1.train,
        F1.test = F1.test,
        kendall.train = kendall.train,
        kendall.test = kendall.test,
        harm.train = harm.train,
        harm.test = harm.test,
        beta = sim$beta,
        gammas = sim$gammas
    ))
}

# Evaluate the predictive performance of a simulated probit model.
# Returns a list of performance metrics (F1 score, Jaccard index, Youden's index, harmonic mean, etc.).
# 
# Arguments:
# sim      : A list containing simulation results with beta coefficients, gamma thresholds,
#            design matrices (X), and observed outcomes (y) for each dataset. This is the object returned
#            by the hprobit_gibbs function.
# test_set : An optional list of test datasets with design matrices (X) and observed outcomes (Y).
eval_sim <- function(sim, test_set = NULL) {
    # Note: no code is using this function. All calls are commented out. Only 
    # eval_sim_small is used.

    # Get the number of posterior draws and datasets
    ndraws <- nrow(sim$beta)
    ndatasets <- length(sim$gammas)
    # Determine the number of result columns (for multi-scale, add one)
    if(ndatasets>1) {
        nresults <- ndatasets+1
    } else {
        nresults <- ndatasets
    }

    # Initialize containers for all metrics and classification performance
    clperf.train <- vector("list", ndatasets)
    clperf.test <- vector("list", ndatasets)
    F1.train <- matrix(0, nrow = ndraws, ncol = nresults)
    F1.test <- matrix(0, nrow = ndraws, ncol = nresults)
    jaccard.train <- matrix(0, nrow = ndraws, ncol = nresults)
    jaccard.test <- matrix(0, nrow = ndraws, ncol = nresults)
    youden.train <- matrix(0, nrow = ndraws, ncol = nresults)
    youden.test <- matrix(0, nrow = ndraws, ncol = nresults)
    clharm.train <- matrix(0, nrow = ndraws, ncol = nresults)
    clharm.test <- matrix(0, nrow = ndraws, ncol = nresults)
    kendall.train <- matrix(0, nrow = ndraws, ncol = nresults)
    kendall.test <- matrix(0, nrow = ndraws, ncol = nresults)
    spearman.train <- matrix(0, nrow = ndraws, ncol = nresults)
    spearman.test <- matrix(0, nrow = ndraws, ncol = nresults)
    rharm.train <- matrix(0, nrow = ndraws, ncol = nresults)
    rharm.test <- matrix(0, nrow = ndraws, ncol = nresults)
    harm.train <- matrix(0, nrow = ndraws, ncol = nresults)
    harm.test <- matrix(0, nrow = ndraws, ncol = nresults)

    # Loop over each dataset (group)
    for (i in 1:ndatasets) {
        # Compute the linear predictor for all draws
        y <- tcrossprod(sim$X[[i]], sim$beta)
        # For each draw, discretize the latent variable into categories using the sampled gammas
        y.train <- sapply(1:ndraws, function(j) {
            return(findInterval(y[, j], sim$gammas[[i]][j, ]))
        })

        # Reference labels (shifted to zero-based)
        ref <- sim$y[[i]]-1
        ref.fac <- factor(ref)

        # Compute all classification metrics for each draw
        clperf.train[[i]] <- cpp_classification_metric_distributions(y.train, ref, length(levels(ref.fac)))
        F1.train[, i] <- clperf.train[[i]]$F1
        jaccard.train[, i] <- clperf.train[[i]]$Jaccard
        youden.train[, i] <- clperf.train[[i]]$Youden
        # (Optional) Compute harmonic means of metrics (commented out)
        # Compute correlation metrics if more than two classes (commented out)

        # If a test set is provided, evaluate on test data
        if (!is.null(test_set)) {
            # Compute the linear predictor for test data
            y <- tcrossprod(test_set[[i]]$X, sim$beta)

            # Discretize using sampled gammas for each draw
            y.test <- sapply(1:ndraws, function(j) {
                return(findInterval(y[, j], sim$gammas[[i]][j, ]))
            })

            # Reference test labels
            ref <- test_set[[i]]$Y
            ref.fac <- factor(ref)

            # Compute all classification metrics for test data
            clperf.test[[i]] <- cpp_classification_metric_distributions(y.test, ref, length(levels(ref.fac)))
            F1.test[, i] <- clperf.test[[i]]$F1
            jaccard.test[, i] <- clperf.test[[i]]$Jaccard
            youden.test[, i] <- clperf.test[[i]]$Youden

            # (Optional) Compute harmonic means of metrics (commented out)
            # Compute Kendall's tau for test data
            kendall.test[, i] <- apply(y.test, 2, function(draw) {
                # return(fastkendall(draw, y = ref)$kendall.tau)
                return(suppressWarnings(cor(draw, ref, method = "kendall")))
            })

            # Compute Spearman's rho for test data
            spearman.test[, i] <- apply(y, 2, function(draw) {
                return(cor(draw, ref, method = "spearman"))
            })

            # (Optional) Compute harmonic means of correlation metrics (commented out)
            # Identify valid draws (Kendall's tau and Spearman's rho >= 0 and not NA)
            ok <- kendall.test[, i] >= 0 & spearman.test[, i] >= 0
            ok[which(is.na(ok))] <- FALSE

            # Compute harmonic mean of all metrics for valid draws
            harm.test[ok, i] <- cpp_harmonic_rowmeans(as.matrix(data.frame(F1 = F1.test[ok, i],
                                                                    jaccard = jaccard.test[ok, i],
                                                                    youden = youden.test[ok, i],
                                                                    kendall = kendall.test[ok, i],
                                                                    spearman = spearman.test[ok, i])))
        }
    }
    # (Optional) Code for computing harmonic means across datasets is commented out
    # Return all computed metrics and posterior samples
    return(list(F1.train = F1.train,
                F1.test = F1.test,
                jaccard.train = jaccard.train,
                jaccard.test = jaccard.test,
                youden.train = youden.train,
                youden.test = youden.test,
                clharm.train = clharm.train,
                clharm.test = clharm.test,
                kendall.train = kendall.train,
                kendall.test = kendall.test,
                spearman.train = spearman.train,
                spearman.test = spearman.test,
                rharm.train = rharm.train,
                rharm.test = rharm.test,
                harm.train = harm.train,
                harm.test = harm.test))
}


# Compute the root square mean error (RMSE) between the true beta coefficients and the sampled 
# beta coefficients. For each simulation, it computes the RMSE for every posterior draw of beta, 
# returning a matrix or list of these errors for all simulations. This helps assess how close the 
# estimated betas are to the true values across all MCMC samples.
# Returns a matrix or list of RMSE values.
# 
# Arguments:
# sims      : A list of simulation results, each containing sampled beta coefficients.
# true_beta : The true beta coefficients used in the data generation.
compute_beta_errors <- function(sims, true_beta) {
  return(sapply(sims, function (sim) {
    return(apply(sim$beta, 1, function(beta) {
      return(sqrt(mean((true_beta-beta)^2)))
    }))
  }))
}


# Reorganize the results from multiple simulation trials (the output of run_n_trials) and restructures 
# them into a more convenient format for analysis. It extracts and organizes:
# - Summary statistics (means and totals) for metrics, beta coefficients, and gamma thresholds across all trials.
# - The true (known) beta and gamma values if the data is simulated.
# - The training and test data (X, Y, y.star) for each group and trial.
# 
# Returns a large list where each element (e.g., metric_means, beta_means, X.train, y.train, etc.) is 
# organized by group and trial.
# 
# Arguments:
# trial     : A list of simulation trial results, where each trial contains metrics, beta/gamma samples, and data.
# simulated : A boolean indicating whether the data is simulated (to extract known beta/gamma).
reorganize <- function(trial, simulated = TRUE) {
  # Get trial, covariate, draw, and gamma group sizes
  ntrials <- length(trial)
  ncov <- length(trial[[1]][[1]][[2]][[1]])
  ndraws <- nrow(trial[[1]][[2]][[2]][[1]])
  ngammas <- sapply(trial[[1]][[1]][[3]], length)

  # Initialize containers for summary statistics and data
  metric_means <- list()
  metric_totals <- list()
  beta_means <- list()
  beta_totals <- list()
  gamma_means <- list()
  gamma_totals <- list()

  # Initialize containers for known parameters if simulated data
  if (simulated) {
    known_beta <- matrix(0, nrow = ntrials, ncol = ncov)
    known_gamma <- list()
  } else {
    known_beta <- NA
    known_gamma <- NA
  }

  # Initialize containers for raw data (X, Y, y.star)
  X.train <- list()
  X.test <- list()
  y.train <- list()
  y.test <- list()
  y_star.train <- list()
  y_star.test <- list()

  # Pre-allocate matrices for each metric/gamma group
  for (j in 1:6) {
    metric_means[[j]] <- matrix(0, nrow = ntrials, ncol = 6,
                                dimnames = list(NULL,
                                                names(trial[[1]][[1]][[1]][[1]])))
    metric_totals[[j]] <- matrix(0, nrow = ntrials*ndraws, ncol = 6,
                                 dimnames = list(NULL,
                                                 names(trial[[1]][[1]][[1]][[1]])))
    gamma_means[[j]] <- matrix(0, nrow = ntrials, ncol = ngammas[j])
    gamma_totals[[j]] <- matrix(0, nrow = ntrials*ndraws, ncol = ngammas[j])
  }
  # Pre-allocate matrices for beta coefficients
  for (j in 1:4) {
    beta_means[[j]] <- matrix(0, nrow = ntrials, ncol = ncov)
    beta_totals[[j]] <- matrix(0, nrow = ntrials*ndraws, ncol = ncov)
  }
  # Pre-allocate containers for each group for known gamma and data
  for (j in 1:3) {
    if (simulated) {
      known_gamma[[j]] <- matrix(0, nrow = ntrials, ncol = length(trial[[1]]$data$train[[j]]$gamma))
    }
    X.train[[j]] <- list()
    X.test[[j]] <- list()
    y.train[[j]] <- list()
    y.test[[j]] <- list()
    y_star.train[[j]] <- list()
    y_star.test[[j]] <- list()
  }

  # Fill in all containers with data from each trial
  for (i in 1:ntrials) {
    # Store known beta for simulated data
    if (simulated) {
      known_beta[i, ] <- trial[[i]]$data$train[[1]]$beta
    }
    # Store metrics and gamma values for each group
    for (j in 1:6) {
      metric_means[[j]][i, ] <- trial[[i]][[1]][[1]][[j]]
      metric_totals[[j]][(i*ndraws - ndraws + 1):(i*ndraws), ]  <- trial[[i]][[2]][[1]][[j]]
      gamma_means[[j]][i, ] <- trial[[i]][[1]][[3]][[j]]
      gamma_totals[[j]][(i*ndraws - ndraws + 1):(i*ndraws), ] <- trial[[i]][[2]][[3]][[j]]
    }
    # Store beta values for each group
    for (j in 1:4) {
      beta_means[[j]][i, ] <- trial[[i]][[1]][[2]][[j]]
      beta_totals[[j]][(i*ndraws - ndraws + 1):(i*ndraws), ] <- trial[[i]][[2]][[2]][[j]]
    }
    # Store known gamma and all data for each group
    for (j in 1:3) {
      if (simulated) {
        known_gamma[[j]][i, ] <- trial[[i]]$data$train[[j]]$gamma
        y_star.train[[j]][[i]] <- trial[[i]]$data$train[[j]]$y.star
        y_star.test[[j]][[i]] <- trial[[i]]$data$test[[j]]$y.star
      }
      X.train[[j]][[i]] <- trial[[i]]$data$train[[j]]$X
      X.test[[j]][[i]] <- trial[[i]]$data$test[[j]]$X
      y.train[[j]][[i]] <- trial[[i]]$data$train[[j]]$Y
      y.test[[j]][[i]] <- trial[[i]]$data$test[[j]]$Y
    }
  }

  # Return all organized results as a list
  return(list(metric_means = metric_means,
              beta_means = beta_means,
              gamma_means = gamma_means,
              metric_totals = metric_totals,
              beta_totals = beta_totals,
              gamma_totals = gamma_totals,
              known_beta = known_beta,
              known_gamma = known_gamma,
              X.train = X.train,
              X.test = X.test,
              y.train = y.train,
              y.test = y.test,
              y_star.train = y_star.train,
              y_star.test = y_star.test))
}
