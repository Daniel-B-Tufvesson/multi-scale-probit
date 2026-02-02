# MOVE FILE TO scripts/experiments_sim.R

# library(coda)
# library(MASS)
# library(caret)
# library(Kendall)
# library(e1071)
# library(fastkendall)

source("R/hprobit_gibbs.R")
source("R/util.R")

# Libraries for parallelization.
library(parallel)
library(doParallel)
library(foreach)


# Function to generate synthetic experiment data. Returns a list with training and test sets.
#
# Arguments:
#   nobs    : Number of observations to generate for the dataset.
#   ncov    : Number of covariates (features) in the dataset.
#   train_id: Indices for training set rows.
#   test_id : Indices for test set rows.
#   data_sd : Standard deviation for generated covariate values (default: 1).
#   beta    : Numeric vector of regression coefficients (optional; generated if NULL).
#   gammas  : List of threshold values for each group (optional; generated if NULL).
#   ngamma  : Vector specifying number of gamma thresholds per group (default: c(1, 3, 3)).
generate_experiment <- function(nobs,
                                ncov,
                                train_id,
                                test_id,
                                data_sd = 1,
                                beta = NULL,
                                gammas = NULL,
                                ngamma = c(1, 3, 3)) {
  # If beta is not provided, generate random beta coefficients
  if (is.null(beta)) {
    beta <- rnorm(ncov, 0, 1)
  }
  # If gammas is not provided, generate random gamma thresholds for each group
  if (is.null(gammas)) {
    gammas <- list()
    for (i in 1:length(ngamma)) {
      gammas[[i]] <- sort(rnorm(ngamma[i], 0, 5))
    }
  }
  # Create covariate names
  X.names <- paste("X", 1:ncov, sep = "")
  train <- list()
  test <- list()
  # For each gamma group, generate data and apply roprobit
  for (i in 1:length(gammas)) {
    # Generate covariate matrix X
    X <- as.matrix(sapply(1:(ncov), function(i) {
      return(rnorm(nobs, sd = data_sd))
    }))
    colnames(X) <- X.names
    
    # Apply roprobit to training and test sets
    tmp_train <- roprobit(X[train_id, ], beta, gammas[[i]])
    tmp_test <- roprobit(X[test_id, ], beta, gammas[[i]])
    # If roprobit fails, regenerate gammas and try again
    while ((!is.list(tmp_train)) | (!is.list(tmp_test))) {
      gammas[[i]] <- sort(rnorm(ngamma[i], 0, 5))
      tmp_train <- roprobit(X[train_id, ], beta, gammas[[i]])
      tmp_test <- roprobit(X[test_id, ], beta, gammas[[i]])
    }
    # Store results
    train[[i]] <- tmp_train
    test[[i]] <- tmp_test
  }
  # Return training and test sets as a list
  return(list(train = train,
              test = test))
}

# -----------------------------------------------------------------------------------------------------
# Test generate data. This is not used for fitting.
set.seed(12345)
beta <- rnorm(48, 0, 1)
gamma_bin <- 3
gamma_rnk1 <- c(-5, 0, 5)
gamma_rnk2 <- c(-7, -3, 2)
gammas <- list(gamma_bin, gamma_rnk1, gamma_rnk2)
dat <- generate_experiment(nobs = 600,
                           ncov = 48,
                           train_id = 1:400,
                           test_id = 401:600,
                           data_sd = 1,
                           beta = beta,
                           gammas = gammas)

sort(dat$train[[1]]$Y)
sort(dat$train[[1]]$Y)
sort(dat$train[[2]]$Y)
sort(dat$train[[3]]$Y)


# ------------------------------------------------------------------------------------------------------
# Test fitting on synthetic data. This is not the final experiment.

# Initial parameter tuning code
burnin <- 50000
ndraws <- 50000
thin <- 100
meanPrior <- rep(0, 48)
precPrior = diag(rep(0.1, 48))
dfPrior <- 10 # nu0
variancePrior <- 1 # sigma20
etaPrior <- 1 # eta_0
lambdaPrior <- 1 #lambda_0
verbose <- 0

# Generate synthetic data for testing.
dat <- generate_experiment(nobs = 60,
                           ncov = 48,
                           train_id = 1:40,
                           test_id = 41:60,
                           data_sd = 1,
                           ngamma = c(1, 3, 3))
train <- dat$train
test <- dat$test
precPrior = NULL
tune_all = list(0.1, 0.1, 0.1, c(5.0, 1.9, 1.9)) # 48 covariates, 3 datasets, 2/3 training, 60 dataset, precPrior = diag(rep(0.1, 48)), no fixed gamma or beta, with estimated shrinkage
tune_lambda = 30

set = 2

# 
if (set == length(tune_all)) {
    # Fit using all training dataset (multi-scale).
    print(system.time(tst_chain <- hprobit_gibbs(train, 
                                                burnin = burnin,
                                                ndraws = ndraws, 
                                                thin = thin,
                                                tune = tune_all[[set]],
                                                meanPrior = meanPrior,
                                                precPrior = diag(rep(0.1, 48)),
                                                verbose = verbose)))
} else {
    # Fit using only one of the training datasets (single-scale).
    print(system.time(tst_chain <- hprobit_gibbs(train[set], 
                                                burnin = burnin,
                                                ndraws = ndraws, 
                                                thin = thin,
                                                tune = tune_all[[set]],
                                                meanPrior = meanPrior,
                                                precPrior = diag(rep(0.1, 48)),
                                                verbose = verbose)))
}
beep(4)


# -------------------------------------------------------------------------------------------------------
# Do full experiment on synthetic data.


# Function to run one complete trial. Returns a list with training set, test set, tunings, and simulation results.
#
# Arguments:
#   testset  : List of test datasets (one per group).
#   trainset : List of training datasets (one per group).
#   tunings  : List of tuning parameters (one per group, plus one for multi-scale).
#   burnin   : Number of burn-in iterations for MCMC.
#   ndraws   : Number of MCMC samples to draw.
#   thin     : Thinning interval for MCMC samples.
#   verbose  : Verbosity level for model fitting output.
run_experiment <- function(testset, trainset, tunings,
                           burnin,
                           ndraws,
                           thin,
                           verbose) {
  # Initialize result list
  result <- list()
  # Fit model for each test set individually
  for (i in 1:length(testset)) {
    sim <- hprobit_gibbs(trainset[i],
                         burnin = burnin,
                         ndraws = ndraws,
                         thin = thin,
                         tune = tunings[[i]],
                         meanPrior = meanPrior,
                         precPrior = precPrior,
                         verbose = verbose)
    print(paste("Set", i))
    result[[i]] <- sim
  }
  # Fit model using all training sets combined (multi-scale model)
  sim <- hprobit_gibbs(trainset,
                       burnin = burnin,
                       ndraws = ndraws,
                       thin = thin,
                       tune = tunings[[length(tunings)]],
                       meanPrior = meanPrior,
                       precPrior = precPrior,
                       verbose = verbose)
  print("All")
  result[[length(tunings)]] <- sim

  # Return results and input parameters as a list
  return(list(trainset = trainset,
              simulations = result,
              testset = testset,
              tunings = tunings,
              ndraws = ndraws,
              thin = thin))
}





# Function to run multiple trials. Returns a list of results for each trial. The trials are run 
# independently in parallel. For each trial it: 
# 1. Generates synthetic training and test datasets using generate_experiment.
# 2. Fits models to the data using run_experiment (which calls hprobit_gibbs).
# 3. Evaluates the results for each group and for the combined model.
# 4. Collects summary statistics (means, totals) for each trial.
# 5. Optionally saves the results to an RDS file if a filename is provided.
# 
# Arguments:
# n         : Number of trials to run.
# tune_all  : List of tuning parameters (one per group, plus one for multi-scale).
# burnin    : Number of burn-in iterations for MCMC.
# ndraws    : Number of MCMC samples to draw.
# meanPrior : Prior mean vector for regression coefficients.
# precPrior : Prior precision matrix for regression coefficients.
# thin      : Thinning interval for MCMC samples (default: 100).
# ncov      : Number of covariates (features) in the dataset (default: 48).
# nobs      : Number of observations to generate for the dataset.
# train_id  : Indices for training set rows.
# test_id   : Indices for test set rows.
# filename  : Optional filename to save results as RDS.
# logfile   : Log file to record progress and results (default: "log.txt").
run_n_trials <- function(n = 1,
                         tune_all,
                         burnin,
                         ndraws,
                         meanPrior,
                         precPrior,
                         thin = 100,
                         ncov = 48,
                         nobs,
                         train_id,
                         test_id,
                         filename = NULL,
                         logfile = "log.txt") {
  
  # Set up gamma group structure
  ngamma <- c(1, 3, 3)
  gamma_lengths <- ngamma
  gamma_lengths <- c(gamma_lengths, gamma_lengths)
  print(gamma_lengths)

  # Calculate and print number of model parameters
  nestimates <- ncov + sum(gamma_lengths)
  print(paste("Estimating", nestimates, "model parameters"))

  # Initialize storage for summary statistics
  means <- list(metrics = list(), beta = list(), gamma = list())
  totals <- list(metrics = list(), beta = list(), gamma = list())

  # Prepare logfile
  writeLines(c(""), logfile)

  # Run n trials in parallel using foreach and %dopar%
  res <- foreach(i=1:n, .export = c("generate_experiment",
                                    "run_experiment",
                                    "roprobit",
                                    "hprobit_gibbs",
                                    "meanPrior",
                                    "precPrior",
                                    "means",
                                    "totals")) %dopar% {
    # Redirect output to logfile
    sink(logfile, append=TRUE)

    # Generate synthetic data for this trial
    dat <- generate_experiment(nobs = nobs,
                               ncov = ncov,
                               train_id = train_id,
                               test_id = test_id,
                               data_sd = 1,
                               ngamma = ngamma)
    train <- dat$train
    test <- dat$test
    print(paste("Running trial", i))
    # Fit models and collect results
    results <- run_experiment(test, train, tune_all,
                              burnin,
                              ndraws,
                              thin,
                              0)
    print("--Evaluating results")
    # Evaluate single scale models for each group
    print("--Evaluating single scale models")
    for (j in 1:length(ngamma)) {
        print(paste("----Evaluating set", j))
        res <- eval_sim_small(results$simulations[[j]], results$testset[j])
        means$metrics[[j]] <- sapply(res[1:6], mean, na.rm = TRUE)
        means$beta[[j]] <- colMeans(as.matrix(res[[7]]), na.rm = TRUE)
        means$gamma[[j]] <- colMeans(as.matrix(res[[8]][[1]]), na.rm = TRUE)
        totals$metrics[[j]] <- matrix(unlist(res[1:6]), ncol = 6, byrow = FALSE, dimnames = list(NULL, names(res[1:6])))
        totals$beta[[j]] <- as.matrix(res[[7]])
        totals$gamma[[j]] <- as.matrix(res[[8]][[1]])
    }
    # Evaluate multi-scale model (all groups combined)
    print("----Evaluating multi scale models")
    res <- eval_sim_small(results$simulations[[length(ngamma)+1]], results$testset)
    metric_means <- sapply(res[1:6], colMeans)
    for (j in 1:length(ngamma)) {
      means$metrics[[length(ngamma)+j]] <- metric_means[j, ]
      totals$metrics[[3+j]] <- sapply(res[1:6], function(metric) {
        return(metric[, j])
      })
      means$gamma[[3+j]] <- colMeans(as.matrix(res[[8]][[j]]), na.rm = TRUE)
      totals$gamma[[3+j]] <- as.matrix(res[[8]][[j]])
    }
    means$beta[[length(ngamma)+1]] <- colMeans(as.matrix(res[[7]]), na.rm = TRUE)
    totals$beta[[length(ngamma)+1]] <- as.matrix(res[[7]])
    # Restore output and return results for this trial
    sink()
    list(means = means, totals = totals, data = dat)
  }

  # Optionally save results to file
  if (!is.null(filename)) {
    saveRDS(res, filename)
  }

  # Return all trial results
  return(res)
}

# Comment out this code since duplicate code is below.
#Finns ej definierade
# cl<-makeForkCluster(8)
# registerDoParallel(cl)
# clusterSetRNGStream(cl, 12345)

# worker.init <- function() {
#   source("hprobit_gibbs.R")
#   source("util.R")
# }
# clusterCall(cl, worker.init)


# Comment out this code since the "parallel" package is not loaded.
# Set up workers for parallelization.
cl<-makeForkCluster(8)
registerDoParallel(cl)
clusterSetRNGStream(cl, 12345)


print(system.time(trial_results60 <- run_n_trials(n = 500,
                                                  tune_all = list(5.0, 1.9, 1.9, c(5.0, 1.9, 1.9)),
                                                  burnin = 50000,
                                                  ndraws = 50000,
                                                  meanPrior = rep(0, 48),
                                                  precPrior = diag(rep(0.1, 48)),
                                                  nobs = 60,
                                                  train_id = 1:40,
                                                  test_id = 41:60,
                                                  filename = "sim_60_dataset_all500.RDS",
                                                  logfile = "sim_60_log.txt")))

# print(system.time(trial_results600 <- run_n_trials(n = 500,
#                               tune_all = list(0.57, 0.2, 0.22, c(0.41, 0.19, 0.19)),
#                               burnin = 50000,
#                               ndraws = 50000,
#                               meanPrior = rep(0, 48),
#                               precPrior = diag(rep(0.1, 48)),
#                               nobs = 600,
#                               train_id = 1:400,
#                               test_id = 401:600,
#                               filename = "sim_600_dataset_all500.RDS",
#                               logfile = "sim_600_log.txt")))

print(system.time(trial_results600 <- run_n_trials(n = 500,
                              tune_all = list(1.0, 0.3, 0.3, c(1.0, 0.3, 0.3)),
                              burnin = 50000,
                              ndraws = 50000,
                              meanPrior = rep(0, 48),
                              precPrior = diag(rep(1, 48)),
                              nobs = 600,
                              train_id = 1:400,
                              test_id = 401:600,
                              filename = "sim_600_dataset_all500_noshrinkage.RDS",
                              logfile = "sim_600_log.txt")))



str(trial_results60[[1]])

#sapply(trial_results[[1]][[1]][[3]], length)

stopCluster(cl)

trial_results60 <- readRDS("sim_60_dataset_all500.RDS")
reorganized60 <- reorganize(trial_results60)
saveRDS(reorganized60, "sim_60_dataset_all500_finished2_reorganized.RDS")

trial_results600 <- readRDS("sim_600_dataset_all500_noshrinkage_finished2.RDS")
reorganized600 <- reorganize(trial_results600)
saveRDS(reorganized600, "sim_600_dataset_all500_noshrinkage_finished2_reorganized.RDS")


