# MOVE FILE TO scripts/experiments_test.R

source("hprobit_gibbs.R")
source("util.R")

library(foreach)
library(doParallel)


# features <- readRDS("featureset47_names.RDS")
# 
# fh2 <- readRDS("nyponvilja_data/nypon_1.RDS")
# fh3 <- readRDS("nyponvilja_data/nypon_2.RDS")
# fh4 <- readRDS("nyponvilja_data/nypon_3.RDS")
# 
# fleg3_9 <- readRDS("nyponvilja_data/vilja_medium.RDS")
# fleg10_12 <- readRDS("nyponvilja_data/vilja_large.RDS")
# fleg13_19 <- readRDS("nyponvilja_data/vilja_x-large.RDS")
# 
# fll <- readRDS("nyponvilja_data/vilja_small.RDS")
# 
# norstedts <- readRDS("nyponvilja_data/nypon_5.RDS")
# bonnier77_78 <- readRDS("nyponvilja_data/nypon_4.RDS")
# bonnier80_81 <- readRDS("nyponvilja_data/nypon_6.RDS")

fh2 <- readRDS("subsampled_data/fh2.RDS")
fh3 <- readRDS("subsampled_data/fh3.RDS")
fh4 <- readRDS("subsampled_data/fh4.RDS")

fleg3_9 <- readRDS("subsampled_data/fleg3_9.RDS")
fleg10_12 <- readRDS("subsampled_data/fleg10_12.RDS")
fleg13_19 <- readRDS("subsampled_data/fleg13_19.RDS")

fll <- readRDS("subsampled_data/fll.RDS")

norstedts <- readRDS("subsampled_data/norstedts.RDS")
bonnier77_78 <- readRDS("subsampled_data/bonnier77_78.RDS")
bonnier80_81 <- readRDS("subsampled_data/bonnier80_81.RDS")

fh.length <- nrow(fh2) + nrow(fh3) + nrow(fh4)
fleg.length <- nrow(fleg3_9) + nrow(fleg10_12) + nrow(fleg13_19)
fll.length <- nrow(fll)
f.length <- fh.length + fleg.length + fll.length
fa.length <- nrow(norstedts) + nrow(bonnier77_78) + nrow(bonnier80_81)

#master <- rbind(fh2, fh3, fh4, fleg3_9, fleg10_12, fleg13_19, fll, norstedts, bonnier77_78, bonnier80_81)
#features = names(master)[1:118]

initData <- function(newdata, label, trainprop = NA, trainsize = NA, scale = FALSE) {
  dataset <- list(trainset = list(),
                  testset = list(),
                  trainprop = trainprop,
                  trainsize = trainsize)
  if (scale) {
    tmpmat <- scale(na.omit(data.matrix(newdata[, features])), center = sc.center, scale = sc.scale)
  } else {
    tmpmat <- na.omit(data.matrix(newdata[, features]))
  }
  
  n <- dim(tmpmat)[1]
  if (!is.na(trainprop)) {
    id <- sample(1:n, floor(n*trainprop))  
  } else {
    if (!is.na(trainsize)) {
      id <- sample(1:n, trainsize)
    } else {
      stop("Either trainprop or trainsize has to be provided")
    }
  }
  
  #print(id)
  dataset$trainset[['X']] <- tmpmat[id, ]
  #print(id)
  #print(nrow(tmpmat[id, ]))
  dataset$trainset[['Y']] <- rep(0, length(id))
  dataset$testset[['X']] <- tmpmat[-id, ]
  if (is.null(nrow(tmpmat[-id, ]))) {
    dataset$testset[['Y']] <- c(label)
  } else {
    dataset$testset[['Y']] <- rep(label, nrow(tmpmat[-id, ]))
  }
  return(dataset)
}

addData <- function(dataset, newdata, label, ratio = 1, trainsize = NA, scale = FALSE) {
  if (scale) {
    tmpmat <- scale(na.omit(data.matrix(newdata[, features])), center = sc.center, scale = sc.scale)
  } else {
    tmpmat <- na.omit(data.matrix(newdata[, features]))
  }
  
  n <- dim(tmpmat)[1]
  #print(dim(tmpmat))
  if (is.na(trainsize)) {
    trainprop <- dataset$trainprop
    trainsize <- dataset$trainsize
    if (!is.na(trainprop)) {
      id <- sample(1:n, floor(n*trainprop))  
    } else {
      if (!is.na(trainsize)) {
        id <- sample(1:n, floor(trainsize*ratio))
      } else {
        stop("Either trainprop or trainsize has to be provided")
      }
    }
  } else {
    id <- sample(1:n, floor(trainsize*ratio))
  }
  
  dataset$trainset[['X']] <- rbind(dataset$trainset[['X']], tmpmat[id, ])
  #print(id)
  #print(nrow(tmpmat[id, ]))
  dataset$trainset[['Y']] <- c(dataset$trainset[['Y']], rep(label, length(id)))
  dataset$testset[['X']] <- rbind(dataset$testset[['X']], tmpmat[-id, ])
  if (is.null(nrow(tmpmat[-id, ]))) {
    dataset$testset[['Y']] <- c(dataset$testset[['Y']], label)
  } else {
    dataset$testset[['Y']] <- c(dataset$testset[['Y']], rep(label, nrow(tmpmat[-id, ])))
  }
  return(dataset)
}



generate_text_experiment <- function(global_trainprop = 2/3) {
  #global_trainprop <- 1/3
  global_trainsize <- NA#1
  
  etr_fiction_ratios <- c(fh.length, fleg.length, fll.length)/f.length
  print(etr_fiction_ratios)
  
  fiction_adult <- rbind(rbind(norstedts, bonnier77_78), bonnier80_81)
  samp <- rmultinom(nrow(fiction_adult), size = 1, prob = etr_fiction_ratios)
  
  fll_adult = fiction_adult[samp[1, ], ]
  fleg_adult = fiction_adult[samp[2, ], ]
  fh_adult = fiction_adult[samp[3, ], ]
  
  fiction_ll <- initData(fll, 0, trainprop = global_trainprop, trainsize = global_trainsize)
  fiction_ll <- addData(fiction_ll, fll_adult, 1)
  #str(fiction_ll)
  #dim(fiction_ll$trainset$X)
  
  fiction_leg <- initData(fleg3_9, 0, trainprop = global_trainprop, trainsize = global_trainsize)
  fiction_leg <- addData(fiction_leg, fleg10_12, 1)
  fiction_leg <- addData(fiction_leg, fleg13_19, 2)
  fiction_leg <- addData(fiction_leg, fleg_adult, 3)
  #str(fiction_leg)
  
  fiction_h <- initData(fh2, 0, trainprop = global_trainprop, trainsize = global_trainsize)
  fiction_h <- addData(fiction_h, fh3, 1)
  fiction_h <- addData(fiction_h, fh4, 2)
  fiction_h <- addData(fiction_h, fh_adult, 3)
  #str(fiction_h)
  
  fiction_train <- list()
  fiction_train[[1]] <- fiction_ll$trainset
  fiction_train[[2]] <- fiction_leg$trainset
  fiction_train[[3]] <- fiction_h$trainset
  
  fiction_test <- list()
  fiction_test[[1]] <- fiction_ll$testset
  fiction_test[[2]] <- fiction_leg$testset
  fiction_test[[3]] <- fiction_h$testset
  
  train <- list()
  train[[1]] <- fiction_ll$trainset
  train[[2]] <- fiction_leg$trainset
  train[[3]] <- fiction_h$trainset
  
  test <- list()
  test[[1]] <- fiction_ll$testset
  test[[2]] <- fiction_leg$testset
  test[[3]] <- fiction_h$testset
  
  scale_master <- do.call(rbind, list(train[[1]]$X,
                                      train[[2]]$X,
                                      train[[3]]$X
  ))
  scaled <- scale(scale_master)
  sc.center <- attr(scaled, "scaled:center")
  sc.scale <- attr(scaled, "scaled:scale")
  print("sc.center:", sc.center)
  print("sc.scale:", sc.scale)
  print("Scaled data:")
  print(scaled)
  
  for (i in 1:3) {
    train[[i]]$X <- scale(train[[i]]$X, center = sc.center, scale = sc.scale)
    test[[i]]$X <- scale(test[[i]]$X, center = sc.center, scale = sc.scale)
  }
  
  saveRDS(train, paste("subsampled_data/train_1third_47_", Sys.Date(), ".RDS", sep=""))
  saveRDS(test, paste("subsampled_data/test_1third_47_", Sys.Date(), ".RDS", sep=""))
  return(list(train = train,
              test = test))
}

set.seed(12345)
dat <- generate_text_experiment(1/3)
train <- dat$train
test <- dat$test

train[[1]]$Y
train[[2]]$Y
train[[3]]$Y

print(train[[3]])



# Initial parameter tuning code
burnin <- 50000
ndraws <- 50000
thin <- 100
meanPrior <- rep(0, 47)
precPrior <- diag(rep(0.1, 47))
verbose <- 0
# tune_all = list(3.6, 1.5, 1.5, c(3.4, 1.6, 1.6)) # 48 covariates, 3 datasets

# tune_all = list(16.0, 3.1, 3.4, c(22.5, 3.6, 4.4)) # 48 covariates, 3 datasets

# tune_all = list(55.0, 10.0, 13.0, c(75.0, 10.0, 15.0)) # 48 covariates, 3 datasets
# 
# set = 4
# if (set == length(tune_all)) {
#   print(system.time(tst_chain <- hprobit_gibbs(train, burnin = burnin, ndraws = ndraws, thin = thin,
#                                                tune = tune_all[[length(tune_all)]], meanPrior = meanPrior,
#                                                precPrior = precPrior, verbose = verbose)))
#   print(system.time(tst_chain <- hprobit_gibbs(train, burnin = burnin, ndraws = ndraws, thin = thin,
#                                                tune = tune_all[[length(tune_all)]], meanPrior = meanPrior,
#                                                precPrior = precPrior, verbose = verbose)))
# } else {
#   print(system.time(tst_chain <- hprobit_gibbs(train[set], burnin = burnin, ndraws = ndraws, thin = thin,
#                                                tune = tune_all[[set]], meanPrior = meanPrior, precPrior = precPrior,
#                                                verbose = verbose)))
#   print(system.time(tst_chain <- hprobit_gibbs(train[set], burnin = burnin, ndraws = ndraws, thin = thin,
#                                                tune = tune_all[[set]], meanPrior = meanPrior, precPrior = precPrior,
#                                                verbose = verbose)))
# }


run_experiment <- function(testset, trainset, tunings,
                           burnin,
                           ndraws,
                           thin,
                           verbose,
                           all = TRUE) {
  result <- list()
  if (all) {
    for (i in 1:length(testset)) {
      print(paste("Set", i))
      sim <- hprobit_gibbs(trainset[i],
                           burnin = burnin,
                           ndraws = ndraws,
                           thin = thin,
                           tune = tunings[[i]],
                           meanPrior = meanPrior,
                           precPrior = precPrior,
                           verbose = verbose)
      result[[i]] <- sim
    }
  }
  print("All")
  sim <- hprobit_gibbs(trainset,
                       burnin = burnin,
                       ndraws = ndraws,
                       thin = thin,
                       tune = tunings[[length(tunings)]],
                       meanPrior = meanPrior,
                       precPrior = precPrior,
                       verbose = verbose)
  result[[length(tunings)]] <- sim
  return(list(trainset = trainset,
              simulations = result,
              testset = testset,
              tunings = tunings,
              ndraws = ndraws,
              thin = thin))
}

run_n_text_trials <- function(n = 1,
                              training.ratio = 2/3,
                              tune_all,
                              burnin,
                              ndraws,
                              meanPrior,
                              precPrior,
                              thin = 100,
                              ncov = 47,
                              filename = NULL,
                              logfile = "text_trial_log.txt") {
  gammas <- list(c(NA), c(NA, NA, NA), c(NA, NA, NA))
  gamma_lengths <- c(1, 3, 3)
  nestimates <- ncov + sum(gamma_lengths)
  gamma_lengths <- c(gamma_lengths, gamma_lengths)
  print(gamma_lengths)
  print(paste("Estimating", nestimates, "model parameters"))
  means <- list(metrics = list(),
                beta = list(),
                gamma = list())
  totals <- list(metrics = list(),
                 beta = list(),
                 gamma = list())
  # for (j in 1:(2*3)) {
  #   means$metrics[[j]] <- matrix(nrow = n, ncol = 6) #based on number of results from eval_sim_small
  #   means$gamma[[j]] <- matrix(nrow = n, ncol = gamma_lengths[j])
  #   totals$metrics[[j]] <- matrix(nrow = n*ndraws/thin, ncol = 6)
  #   totals$gamma[[j]] <- matrix(nrow = n*ndraws/thin, ncol = gamma_lengths[j])
  # }
  # for (j in 1:(3+1)) {
  #   means$beta[[j]] <- matrix(nrow = n, ncol = ncov) #based on number of results from eval_sim_small
  #   totals$beta[[j]] <- matrix(nrow = n*ndraws/thin, ncol = ncov)
  # }
  writeLines(c(""), logfile)
  res <- foreach(i=1:n, .export = c("generate_text_experiment",
                                    "initData",
                                    "addData",
                                    "fh.length",
                                    "fleg.length",
                                    "fll.length",
                                    "f.length",
                                    "fa.length",
                                    "norstedts",
                                    "bonnier77_78",
                                    "bonnier80_81",
                                    "fh2",
                                    "fh3",
                                    "fh4",
                                    "fleg3_9",
                                    "fleg10_12",
                                    "fleg13_19",
                                    "fll",
                                    "features",
                                    "run_experiment",
                                    "hprobit_gibbs",
                                    "meanPrior",
                                    "precPrior",
                                    "means",
                                    "totals")) %dopar% {
                                      sink(logfile, append=TRUE)
                                      dat <- generate_text_experiment(training.ratio)
                                      train <- dat$train
                                      test <- dat$test
                                      print(paste("Running trial", i))
                                      results <- run_experiment(test, train, tune_all,
                                                                burnin,
                                                                ndraws,
                                                                thin,
                                                                0)
                                      #print("--Evaluating results")
                                      # res <- eval_sim(results$simulations[[4]], results$testset)
                                      # print(str(results$simulations[[1]]))
                                      # print(str(results$simulations[[2]]))
                                      # print(str(results$simulations[[3]]))
                                      # print(str(results$simulations[[4]]))
                                      #print("--Evaluating single scale models")
                                      for (j in 1:length(gammas)) {
                                        #print(paste("----Evaluating set", j))
                                        res <- eval_sim_small(results$simulations[[j]], results$testset[j])
                                        means$metrics[[j]] <- sapply(res[1:6], mean, na.rm = TRUE)
                                        means$beta[[j]] <- colMeans(as.matrix(res[[7]]), na.rm = TRUE)
                                        means$gamma[[j]] <- colMeans(as.matrix(res[[8]][[1]]), na.rm = TRUE)
                                        totals$metrics[[j]] <- matrix(unlist(res[1:6]), ncol = 6, byrow = FALSE, dimnames = list(NULL, names(res[1:6])))
                                        totals$beta[[j]] <- as.matrix(res[[7]])
                                        totals$gamma[[j]] <- as.matrix(res[[8]][[1]])
                                      }
                                      #print("----Evaluating multi scale models")
                                      res <- eval_sim_small(results$simulations[[length(gammas)+1]], results$testset)
                                      metric_means <- sapply(res[1:6], colMeans)
                                      for (j in 1:length(gammas)) {
                                        means$metrics[[length(gammas)+j]] <- metric_means[j, ]
                                        totals$metrics[[3+j]] <- sapply(res[1:6], function(metric) {
                                          return(metric[, j])
                                        })
                                        means$gamma[[3+j]] <- colMeans(as.matrix(res[[8]][[j]]), na.rm = TRUE)
                                        totals$gamma[[3+j]] <- as.matrix(res[[8]][[j]])
                                      }
                                      means$beta[[length(gammas)+1]] <- colMeans(as.matrix(res[[7]]), na.rm = TRUE)
                                      totals$beta[[length(gammas)+1]] <- as.matrix(res[[7]])
                                      sink()
                                      list(means = means, totals = totals, data = dat)
                                    }
  
  if (!is.null(filename)) {
    saveRDS(res, filename)
  }
  return(res)
}


cl<-makeForkCluster(20)
registerDoParallel(cl)
clusterSetRNGStream(cl, 12345)
# worker.init <- function() {
#   source("hprobit_gibbs.R")
#   source("util.R")
# }
# clusterCall(cl, worker.init)
print(system.time(trial_results <- run_n_text_trials(n = 2,
                                                     training.ratio = 2/3,
                                                     #tune_all = list(16.0, 3.1, 3.4, c(22.5, 3.6, 4.4)),
                                                     tune_all = list(1, 1, 1, c(1, 1, 1)),
                                                     burnin = 50000,
                                                     ndraws = 25000,
                                                     meanPrior = meanPrior,
                                                     precPrior = precPrior,
                                                     filename = "text_2thirds_dataset_all2_2.RDS")))

print(system.time(trial_results <- run_n_text_trials(n = 2,
                                                     training.ratio = 1/3,
                                                     #tune_all = list(55.0, 10.0, 13.0, c(75.0, 10.0, 15.0)),
                                                     #tune_all = list(16.0, 3.1, 3.4, c(22.5, 3.6, 4.4)),
                                                     tune_all = list(1, 1, 1, c(1, 1, 1)),
                                                     burnin = 50000,
                                                     ndraws = 25000,
                                                     meanPrior = meanPrior,
                                                     precPrior = precPrior,
                                                     filename = "text_1thirds_dataset_all2_2.RDS")))
stopCluster(cl)


trial_results_text <- readRDS("text_2thirds_dataset_all2.RDS")
trial_results_text2 <- readRDS("text_2thirds_dataset_all2_2.RDS")
trial_results_text22 <- readRDS("text_2thirds_dataset_all2_2_2.RDS")
reorganized_text <- reorganize(trial_results_text22, simulated = FALSE)
print(trial_results_text2[[1]][[1]][[3]])
saveRDS(reorganized_text, "text_2thirds_dataset_all2_2_2_reorganized.RDS")

trial_results_text <- readRDS("text_1thirds_dataset_all500.RDS")
trial_results_text2 <- readRDS("text_1thirds_dataset_all2_2.RDS")
trial_results_text22 <- readRDS("text_1thirds_dataset_all2_2_2.RDS")
reorganized_text <- reorganize(trial_results_text22, simulated = FALSE)
saveRDS(reorganized_text, "text_1thirds_dataset_all2_2_reorganized.RDS")



# run_n_text_trials_old <- function(n = 1,
#                          tune_all = list(16.0, 3.1, 3.4, c(22.5, 3.6, 4.4)),
#                          burnin = 50000,
#                          ndraws = 500000,
#                          thin = 100,
#                          filename = NULL,
#                          filename2 = NULL) {
#   model_means <- list()
#   model_totals <- list()
#   for (j in 1:(2*3)) {
#     model_means[[j]] <- matrix(nrow = n, ncol = 6) #based on number of results from eval_sim
#     model_totals[[j]] <- matrix(nrow = n*ndraws/thin, ncol = 6)
#   }
#   for (i in 1:n) {
#     dat <- generate_text_experiment(1/3)
#     train <- dat$train
#     test <- dat$test
#     print(paste("Running trial", i))
#     results <- run_experiment(test, train, tune_all,
#                               burnin,
#                               ndraws,
#                               thin,
#                               0)
#     print("--Evaluating results")
#     # res <- eval_sim(results$simulations[[4]], results$testset)
#     # print(str(results$simulations[[1]]))
#     # print(str(results$simulations[[2]]))
#     # print(str(results$simulations[[3]]))
#     # print(str(results$simulations[[4]]))
#     print("--Evaluating single scale models")
#     for (j in 1:3) {
#       print(paste("----Evaluating set", j))
#       res <- eval_sim_small(results$simulations[[j]], results$testset[j])
#       means <- sapply(res, mean, na.rm = TRUE)
#       # print(dim(model_totals[[j]][((i-1)*ndraws/thin+1):(i*ndraws/thin), ]))
#       # print(dim(matrix(unlist(res), ncol = 16, byrow = FALSE)))
#       model_means[[j]][i, ] <- means
#       model_totals[[j]][((i-1)*ndraws/thin+1):(i*ndraws/thin), ] <- matrix(unlist(res), ncol = 6, byrow = FALSE)
#     }
#     print("----Evaluating multi scale models")
#     res <- eval_sim_small(results$simulations[[3+1]], results$testset)
#     means <- sapply(res, colMeans)
#     for (j in 1:3) {
#       model_means[[3+j]][i, ] <- means[j, ]
#       for (k in 1:6) {
#         # print((i-1)*ndraws/thin+1:i*ndraws/thin)
#         # print(length(res[[k]][, j]))
#         model_totals[[3+j]][((i-1)*ndraws/thin+1):(i*ndraws/thin), k] <- res[[k]][, j]
#       }
#     }
#     if (!is.null(filename)) {
#       saveRDS(model_means, filename)
#     }
#     if (!is.null(filename2)) {
#       saveRDS(model_totals, filename2)
#     }
#   }
#   return(list(model_means, model_totals))
# }






set.seed(12345)
dat <- generate_text_experiment(1)
train <- dat$train
test <- dat$test

train[[1]]$Y
train[[2]]$Y
train[[3]]$Y


# Initial parameter tuning code
burnin <- 200000
ndraws <- 50000
thin <- 100
meanPrior <- rep(0, 47)
precPrior <- diag(rep(0.1, 47))
verbose <- 0


# tune_all = list(12.3, 1.8, 2.1, c(15.8, 1.7, 2.3)) # 48 covariates, 3 datasets, diag(rep(1, 48))

# tune_all = list(4.09, 0.88, 1.06, c(5.2, 0.9, 1.3)) # 48 covariates, 3 datasets, diag(rep(10, 48))

tune_all = list(20.0, 1.8, 2.1, c(15.8, 1.7, 2.3)) # 48 covariates, 3 datasets, diag(rep(0.1, 48))

set = 1
if (set == length(tune_all)) {
  print(system.time(tst_chain <- hprobit_gibbs(train, burnin = burnin, ndraws = ndraws, thin = thin,
                                               tune = tune_all[[length(tune_all)]], meanPrior = meanPrior,
                                               precPrior = precPrior, verbose = verbose)))
  print(system.time(tst_chain <- hprobit_gibbs(train, burnin = burnin, ndraws = ndraws, thin = thin,
                                               tune = tune_all[[length(tune_all)]], meanPrior = meanPrior,
                                               precPrior = precPrior, verbose = verbose)))
  beep(4)
} else {
  print(system.time(tst_chain <- hprobit_gibbs(train[set], burnin = burnin, ndraws = ndraws, thin = thin,
                                               tune = tune_all[[set]], meanPrior = meanPrior, precPrior = precPrior,
                                               verbose = verbose)))
  print(system.time(tst_chain <- hprobit_gibbs(train[set], burnin = burnin, ndraws = ndraws, thin = thin,
                                               tune = tune_all[[set]], meanPrior = meanPrior, precPrior = precPrior,
                                               verbose = verbose)))
  beep(4)
}



run_n_chains <- function(n = 1,
                         training.ratio = 2/3,
                         tune_all,
                         burnin,
                         ndraws,
                         meanPrior,
                         precPrior,
                         thin = 100,
                         ncov = 47,
                         filename = NULL,
                         logfile = "text_trial_log.txt") {
  gammas <- list(c(NA), c(NA, NA, NA), c(NA, NA, NA))
  gamma_lengths <- c(1, 3, 3)
  nestimates <- ncov + sum(gamma_lengths)
  gamma_lengths <- c(gamma_lengths, gamma_lengths)
  print(gamma_lengths)
  print(paste("Estimating", nestimates, "model parameters"))
  writeLines(c(""), logfile)
  res <- foreach(i=1:n, .export = c("generate_text_experiment",
                                    "initData",
                                    "addData",
                                    "fh.length",
                                    "fleg.length",
                                    "fll.length",
                                    "f.length",
                                    "fa.length",
                                    "norstedts",
                                    "bonnier77_78",
                                    "bonnier80_81",
                                    "fh2",
                                    "fh3",
                                    "fh4",
                                    "fleg3_9",
                                    "fleg10_12",
                                    "fleg13_19",
                                    "fll",
                                    "features",
                                    "run_experiment",
                                    "hprobit_gibbs",
                                    "meanPrior",
                                    "precPrior")) %dopar% {
                                      sink(logfile, append=TRUE)
                                      dat <- generate_text_experiment(training.ratio)
                                      train <- dat$train
                                      test <- dat$test
                                      print(paste("Running trial", i))
                                      results <- run_experiment(test, train, tune_all,
                                                                burnin,
                                                                ndraws,
                                                                thin,
                                                                0,
                                                                all = FALSE)
                                      sink()
                                      results
                                    }
  
  if (!is.null(filename)) {
    saveRDS(res, filename)
  }
  return(res)
}




print(system.time(trial_results <- run_n_chains(n = 8,
                                                training.ratio = 1,
                                                tune_all = list(4.09, 0.88, 1.06, c(5.2, 0.9, 1.3)),
                                                # tune_all = list(12.3, 1.8, 2.1, c(15.8, 1.7, 2.3)),
                                                burnin = 200000,
                                                ndraws = 200000,
                                                meanPrior = meanPrior,
                                                precPrior = precPrior,
                                                filename = "text_full_dataset_all8.RDS")))

saveRDS(trial_results, "8_chains_fitted_on_all_text_data_prec10.RDS")

str(trial_results[1])


str(trial_results[[1]]$simulations[[4]])


