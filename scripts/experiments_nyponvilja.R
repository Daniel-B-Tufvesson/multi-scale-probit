# MOVE FILE TO scripts/experiments_nyponvilja.R

source("hprobit_gibbs.R")
source("util.R")

library(foreach)
library(doParallel)

# Load features
features <- readRDS("featureset48_names.RDS")

# Load data
nypon1 <- readRDS("nyponvilja_data/nypon_1.RDS")
nypon2 <- readRDS("nyponvilja_data/nypon_2.RDS")
nypon3 <- readRDS("nyponvilja_data/nypon_3.RDS")
nypon4 <- readRDS("nyponvilja_data/nypon_4.RDS")
nypon5 <- readRDS("nyponvilja_data/nypon_5.RDS")

viljaxs <- readRDS("nyponvilja_data/vilja_x-small.RDS")
viljas <- readRDS("nyponvilja_data/vilja_small.RDS")
viljam <- readRDS("nyponvilja_data/vilja_medium.RDS")
viljal <- readRDS("nyponvilja_data/vilja_large.RDS")
viljaxl <- readRDS("nyponvilja_data/vilja_x-large.RDS")
viljaxxl <- readRDS("nyponvilja_data/vilja_xx-large.RDS")

suc0_63 <- readRDS("suc_data/suc0-63.RDS")
suc64_126 <- readRDS("suc_data/suc64-126.RDS")

nypon.length <- nrow(nypon1) + nrow(nypon2) + nrow(nypon3) + nrow(nypon4) + nrow(nypon5)
vilja.length <- nrow(viljaxs) + nrow(viljas) + nrow(viljam) + nrow(viljal) + nrow(viljaxl) + nrow(viljaxxl)

ll.length <- nypon.length + vilja.length
suc.length <- nrow(suc0_63) + nrow(suc64_126)

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
  
  etr_ratios <- c(nypon.length, vilja.length)/ll.length
  print(etr_ratios)
  
  fiction_suc <- rbind(suc0_63, suc64_126)
  samp <- rmultinom(nrow(fiction_suc), size = 1, prob = etr_ratios)
  
  nypon_adult = fiction_suc[samp[1, ], ]
  vilja_adult = fiction_suc[samp[2, ], ]
  
  nypon_data <- initData(nypon1, 0, trainprop = global_trainprop, trainsize = global_trainsize)
  nypon_data <- addData(nypon_data, nypon2, 1)
  nypon_data <- addData(nypon_data, nypon3, 2)
  nypon_data <- addData(nypon_data, nypon4, 3)
  nypon_data <- addData(nypon_data, nypon5, 4)
  nypon_data <- addData(nypon_data, nypon_adult, 5)
  #str(nypon_data)
  #dim(nypon_data$trainset$X)
  
  vilja_data <- initData(viljaxs, 0, trainprop = global_trainprop, trainsize = global_trainsize)
  vilja_data <- addData(vilja_data, viljas, 1)
  vilja_data <- addData(vilja_data, viljam, 2)
  vilja_data <- addData(vilja_data, viljal, 3)
  vilja_data <- addData(vilja_data, viljaxl, 4)
  vilja_data <- addData(vilja_data, viljaxxl, 5)
  vilja_data <- addData(vilja_data, vilja_adult, 6)
  #str(vilja_data)
  #dim(vilja_data$trainset$X)
  
  fiction_train <- list()
  fiction_train[[1]] <- nypon_data$trainset
  fiction_train[[2]] <- vilja_data$trainset
  
  fiction_test <- list()
  fiction_test[[1]] <- nypon_data$testset
  fiction_test[[2]] <- vilja_data$testset
  
  train <- list()
  train[[1]] <- nypon_data$trainset
  train[[2]] <- vilja_data$trainset
  
  test <- list()
  test[[1]] <- nypon_data$testset
  test[[2]] <- vilja_data$testset
  
  scale_master <- do.call(rbind, list(train[[1]]$X,
                                      train[[2]]$X
  ))
  scaled <- scale(scale_master)
  sc.center <- attr(scaled, "scaled:center")
  sc.scale <- attr(scaled, "scaled:scale")
  
  for (i in 1:2) {
    train[[i]]$X <- scale(train[[i]]$X, center = sc.center, scale = sc.scale)
    test[[i]]$X <- scale(test[[i]]$X, center = sc.center, scale = sc.scale)
  }
  
  saveRDS(train, paste("nyponvilja_data/train_1third_48_", Sys.Date(), ".RDS", sep=""))
  saveRDS(test, paste("nyponvilja_data/test_1third_48_", Sys.Date(), ".RDS", sep=""))
  return(list(train = train,
              test = test))
}

set.seed(12345)
dat <- generate_text_experiment(1/3)
train <- dat$train
test <- dat$test


train[[1]]$Y
train[[2]]$Y


# Initial parameter tuning code
burnin <- 50000
ndraws <- 50000
thin <- 100
meanPrior <- rep(0, 48)
precPrior <- diag(rep(0.1, 48))
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
                              ncov = 48,
                              filename = NULL,
                              logfile = "text_trial_log.txt") {
  gammas <- list(c(NA, NA, NA, NA, NA), c(NA, NA, NA, NA, NA, NA)) # class labels per scale
  gamma_lengths <- c(5, 6)
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
                                    "nypon.length",
                                    "vilja.length",
                                    "suc0_63",
                                    "suc64_126",
                                    "nypon1",
                                    "nypon2",
                                    "nypon3",
                                    "nypon4",
                                    "nypon5",
                                    "viljaxs",
                                    "viljas",
                                    "viljam",
                                    "viljal",
                                    "viljaxl",
                                    "viljaxxl",
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
                                        means$metrics[[j]] <- sapply(res[1:4], mean, na.rm = TRUE)
                                        means$beta[[j]] <- colMeans(as.matrix(res[[7]]), na.rm = TRUE)
                                        means$gamma[[j]] <- colMeans(as.matrix(res[[8]][[1]]), na.rm = TRUE)
                                        totals$metrics[[j]] <- matrix(unlist(res[1:4]), ncol = 4, byrow = FALSE, dimnames = list(NULL, names(res[1:4])))
                                        totals$beta[[j]] <- as.matrix(res[[7]])
                                        totals$gamma[[j]] <- as.matrix(res[[8]][[1]])
                                      }
                                      #print("----Evaluating multi scale models")
                                      res <- eval_sim_small(results$simulations[[length(gammas)+1]], results$testset)
                                      metric_means <- sapply(res[1:4], colMeans)
                                      for (j in 1:length(gammas)) {
                                        means$metrics[[length(gammas)+j]] <- metric_means[j, ]
                                        totals$metrics[[2+j]] <- sapply(res[1:4], function(metric) {
                                          return(metric[, j])
                                        })
                                        means$gamma[[2+j]] <- colMeans(as.matrix(res[[8]][[j]]), na.rm = TRUE)
                                        totals$gamma[[2+j]] <- as.matrix(res[[8]][[j]])
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
print(system.time(trial_results <- run_n_text_trials(n = 500,
                                                     training.ratio = 2/3,
                                                     tune_all = list(1, 1, c(1, 1)), #1 for now to initialize
                                                     burnin = 50000,
                                                     ndraws = 25000, #50000
                                                     meanPrior = meanPrior,
                                                     precPrior = precPrior,
                                                     filename = "text_2thirds_dataset_all500.RDS")))

print(system.time(trial_results <- run_n_text_trials(n = 500,
                                                     training.ratio = 1/3,
                                                     tune_all = list(1, 1, c(1, 1)), #1 for now to initialize
                                                     burnin = 50000,
                                                     ndraws = 50000, #50000
                                                     meanPrior = meanPrior,
                                                     precPrior = precPrior,
                                                     filename = "text_1thirds_dataset_all500.RDS")))
stopCluster(cl)


trial_results_text <- readRDS("text_2thirds_dataset_all500.RDS")
reorganized_text <- reorganize(trial_results_text, simulated = FALSE)
saveRDS(reorganized_text, "text_2thirds_dataset_all500_reorganized.RDS")

trial_results_text <- readRDS("text_1thirds_dataset_all500.RDS")
reorganized_text <- reorganize(trial_results, simulated = FALSE)
saveRDS(reorganized_text, "text_1thirds_dataset_all500_reorganized.RDS")



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


# Initial parameter tuning code
burnin <- 200000
ndraws <- 50000
thin <- 100
meanPrior <- rep(0, 48)
precPrior <- diag(rep(0.1, 48))
verbose <- 0


# tune_all = list(12.3, 1.8, 2.1, c(15.8, 1.7, 2.3)) # 48 covariates, 3 datasets, diag(rep(1, 48))

# tune_all = list(4.09, 0.88, 1.06, c(5.2, 0.9, 1.3)) # 48 covariates, 3 datasets, diag(rep(10, 48))

##tune_all = list(20.0, 1.8, 2.1, c(15.8, 1.7, 2.3)) # 48 covariates, 3 datasets, diag(rep(0.1, 48))
tune_all = list(1, 1, c(1, 1))

set = 1
if (set == length(tune_all)) {
  print(system.time(tst_chain <- hprobit_gibbs(train, burnin = burnin, ndraws = ndraws, thin = thin,
                                               tune = tune_all[[length(tune_all)]], 
                                               meanPrior = meanPrior,
                                               precPrior = precPrior, verbose = verbose)))
  print(system.time(tst_chain <- hprobit_gibbs(train, burnin = burnin, ndraws = ndraws, thin = thin,
                                               tune = tune_all[[length(tune_all)]], 
                                               meanPrior = meanPrior,
                                               precPrior = precPrior, verbose = verbose)))
  beep(4)
} else {
  print(system.time(tst_chain <- hprobit_gibbs(train[set], burnin = burnin, ndraws = ndraws, thin = thin,
                                               tune = tune_all[[set]], 
                                               meanPrior = meanPrior, precPrior = precPrior,
                                               verbose = verbose)))
  print(system.time(tst_chain <- hprobit_gibbs(train[set], burnin = burnin, ndraws = ndraws, thin = thin,
                                               tune = tune_all[[set]], 
                                               meanPrior = meanPrior, precPrior = precPrior,
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
                         ncov = 48,
                         filename = NULL,
                         logfile = "text_trial_log.txt") {
  gammas <- list(c(NA, NA, NA, NA, NA), c(NA, NA, NA, NA, NA, NA)) #class labels per scale
  gamma_lengths <- c(5, 6)
  nestimates <- ncov + sum(gamma_lengths)
  gamma_lengths <- c(gamma_lengths, gamma_lengths)
  print(gamma_lengths)
  print(paste("Estimating", nestimates, "model parameters"))
  writeLines(c(""), logfile)
  res <- foreach(i=1:n, .export = c("generate_text_experiment",
                                    "initData",
                                    "addData",
                                    "nypon.length",
                                    "vilja.length",
                                    "suc0_63",
                                    "suc64_126",
                                    "nypon1",
                                    "nypon2",
                                    "nypon3",
                                    "nypon4",
                                    "nypon5",
                                    "viljaxs",
                                    "viljas",
                                    "viljam",
                                    "viljal",
                                    "viljaxl",
                                    "viljaxxl",
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
                                      results <- run_experiment(test, train, #tune_all,
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
                                                #tune_all = list(4.09, 0.88, 1.06, c(5.2, 0.9, 1.3)),
                                                # tune_all = list(12.3, 1.8, 2.1, c(15.8, 1.7, 2.3)),
                                                tune_all = list(1, 1, c(1, 1)),
                                                burnin = 200000,
                                                ndraws = 200000,
                                                meanPrior = meanPrior,
                                                precPrior = precPrior,
                                                filename = "text_full_dataset_all8.RDS")))

saveRDS(trial_results, "8_chains_fitted_on_all_text_data_prec10.RDS")

str(trial_results[1])


str(trial_results[[1]]$simulations[[4]])


