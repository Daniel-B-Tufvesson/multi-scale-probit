library(Rcpp)
sourceCpp("rutil.cpp")

roprobit <- function(X, beta, gamma = 0, limit = 10) {
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

rowMin <- function(df) {
  apply(df, 1, min)
}

eval_sim_small <- function(sim, test_set = NULL) {
  ndraws <- nrow(sim$beta)
  ndatasets <- length(sim$gammas)
  if(ndatasets>1) {
    nresults <- ndatasets+1
  } else {
    nresults <- ndatasets
  }
  # clperf.train <- vector("list", ndatasets)
  # clperf.test <- vector("list", ndatasets)
  F1.train <- matrix(0, nrow = ndraws, ncol = nresults)
  F1.test <- matrix(0, nrow = ndraws, ncol = nresults)
  # jaccard.train <- matrix(0, nrow = ndraws, ncol = nresults)
  # jaccard.test <- matrix(0, nrow = ndraws, ncol = nresults)
  # youden.train <- matrix(0, nrow = ndraws, ncol = nresults)
  # youden.test <- matrix(0, nrow = ndraws, ncol = nresults)
  # clharm.train <- matrix(0, nrow = ndraws, ncol = nresults)
  # clharm.test <- matrix(0, nrow = ndraws, ncol = nresults)
  kendall.train <- matrix(0, nrow = ndraws, ncol = nresults)
  kendall.test <- matrix(0, nrow = ndraws, ncol = nresults)
  # spearman.train <- matrix(0, nrow = ndraws, ncol = nresults)
  # spearman.test <- matrix(0, nrow = ndraws, ncol = nresults)
  # rharm.train <- matrix(0, nrow = ndraws, ncol = nresults)
  # rharm.test <- matrix(0, nrow = ndraws, ncol = nresults)
  harm.train <- matrix(0, nrow = ndraws, ncol = nresults)
  harm.test <- matrix(0, nrow = ndraws, ncol = nresults)
  for (i in 1:ndatasets) {
    y <- tcrossprod(sim$X[[i]], sim$beta)
    y.train <- sapply(1:ndraws, function(j) {
      return(findInterval(y[, j], sim$gammas[[i]][j, ]))
    })
    ref <- sim$y[[i]]-1
    ref.fac <- factor(ref)
    F1.train[, i] <- cpp_fmeasure_distribution(y.train, ref, length(levels(ref.fac)))
    kendall.train[, i] <- apply(y.train, 2, function(draw) {
      # return(fastkendall(draw, y = ref)$kendall.tau)
      return(suppressWarnings(cor(draw, ref, method = "kendall")))
    })

    if (!is.null(test_set)) {
      y <- tcrossprod(test_set[[i]]$X, sim$beta)
      y.test <- sapply(1:ndraws, function(j) {
        return(findInterval(y[, j], sim$gammas[[i]][j, ]))
      })
      ref <- test_set[[i]]$Y
      ref.fac <- factor(ref)
      F1.test[, i] <- cpp_fmeasure_distribution(y.test, ref, length(levels(ref.fac)))
      kendall.test[, i] <- apply(y.test, 2, function(draw) {
        return(suppressWarnings(cor(draw, ref, method = "kendall")))
      })
      ok <- kendall.test[, i] >= 0
      ok[which(is.na(ok))] <- FALSE
      harm.test[ok, i] <- cpp_harmonic_rowmeans(as.matrix(data.frame(F1 = F1.test[ok, i],
                                                                     kendall = kendall.test[ok, i])))
    }
  }
  # if(ndatasets>1) {
  #   F1.train[, nresults] <- cpp_harmonic_rowmeans(F1.train[, 1:nresults-1])
  #   F1.test[, nresults] <- cpp_harmonic_rowmeans(F1.test[, 1:nresults-1])
  #   jaccard.train[, nresults] <- cpp_harmonic_rowmeans(jaccard.train[, 1:nresults-1])
  #   jaccard.test[, nresults] <- cpp_harmonic_rowmeans(jaccard.test[, 1:nresults-1])
  #   youden.train[, nresults] <- cpp_harmonic_rowmeans(youden.train[, 1:nresults-1])
  #   youden.test[, nresults] <- cpp_harmonic_rowmeans(youden.test[, 1:nresults-1])
  #   clharm.train[, nresults] <- cpp_harmonic_rowmeans(clharm.train[, 1:nresults-1])
  #   clharm.test[, nresults] <- cpp_harmonic_rowmeans(clharm.test[, 1:nresults-1])
  #   kendall.train[, nresults] <- cpp_harmonic_rowmeans(kendall.train[, 1:nresults-1])
  #   kendall.test[, nresults] <- cpp_harmonic_rowmeans(kendall.test[, 1:nresults-1])
  #   spearman.train[, nresults] <- cpp_harmonic_rowmeans(spearman.train[, 1:nresults-1])
  #   spearman.test[, nresults] <- cpp_harmonic_rowmeans(spearman.test[, 1:nresults-1])
  #   rharm.train[, nresults] <- cpp_harmonic_rowmeans(rharm.train[, 1:nresults-1])
  #   rharm.test[, nresults] <- cpp_harmonic_rowmeans(rharm.test[, 1:nresults-1])
  #   harm.train[, nresults] <- cpp_harmonic_rowmeans(harm.train[, 1:nresults-1])
  #   harm.test[, nresults] <- cpp_harmonic_rowmeans(harm.test[, 1:nresults-1])
  # }
  return(list(F1.train = F1.train,
                F1.test = F1.test,
                kendall.train = kendall.train,
                kendall.test = kendall.test,
                harm.train = harm.train,
                harm.test = harm.test,
                beta = sim$beta,
                gammas = sim$gammas))
}

eval_sim <- function(sim, test_set = NULL) {
  ndraws <- nrow(sim$beta)
  ndatasets <- length(sim$gammas)
  if(ndatasets>1) {
    nresults <- ndatasets+1
  } else {
    nresults <- ndatasets
  }
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
  for (i in 1:ndatasets) {
    y <- tcrossprod(sim$X[[i]], sim$beta)
    y.train <- sapply(1:ndraws, function(j) {
      return(findInterval(y[, j], sim$gammas[[i]][j, ]))
    })
    ref <- sim$y[[i]]-1
    ref.fac <- factor(ref)
    clperf.train[[i]] <- cpp_classification_metric_distributions(y.train, ref, length(levels(ref.fac)))
    F1.train[, i] <- clperf.train[[i]]$F1
    jaccard.train[, i] <- clperf.train[[i]]$Jaccard
    youden.train[, i] <- clperf.train[[i]]$Youden
    # clharm.train[, i] <- cpp_harmonic_rowmeans(as.matrix(data.frame(F1 = F1.train[, i],
    #                                                                 jaccard = jaccard.train[, i],
    #                                                                 youden = youden.train[, i])))
    # if (length(ref.fac) < 3) {
    #   cor.train[, i] <- apply(y, 2, function(draw) {
    #     if (all(draw == ref)) {
    #       return(1)
    #     } else {
    #       return(0)
    #     }
    #   })
    # } else {
    # kendall.train[, i] <- apply(y.train, 2, function(draw) {
    #   # return(fastkendall(draw, y = ref)$kendall.tau)
    #   return(suppressWarnings(cor(draw, ref, method = "kendall")))
    # })
    # spearman.train[, i] <- apply(y, 2, function(draw) {
    #   return(cor(draw, ref, method = "spearman"))
    # })
    # rharm.train[, i] <- cpp_harmonic_rowmeans(as.matrix(data.frame(kendall = kendall.train[, i],
    #                                                                spearman = spearman.train[, i])))
    # df <- data.frame(F1 = F1.train[, i],
    #                  jaccard = jaccard.train[, i],
    #                  youden = youden.train[, i],
    #                  kendall = kendall.train[, i],
    #                  spearman = spearman.train[, i])
    # df <- df[!(rowMin(df) < 0), ]
    # harm.train[, i] <- cpp_harmonic_rowmeans(as.matrix(na.omit(df)))
    # }
    
    if (!is.null(test_set)) {
      y <- tcrossprod(test_set[[i]]$X, sim$beta)
      y.test <- sapply(1:ndraws, function(j) {
        return(findInterval(y[, j], sim$gammas[[i]][j, ]))
      })
      ref <- test_set[[i]]$Y
      ref.fac <- factor(ref)
      clperf.test[[i]] <- cpp_classification_metric_distributions(y.test, ref, length(levels(ref.fac)))
      F1.test[, i] <- clperf.test[[i]]$F1
      jaccard.test[, i] <- clperf.test[[i]]$Jaccard
      youden.test[, i] <- clperf.test[[i]]$Youden
      # clharm.test[, i] <- cpp_harmonic_rowmeans(as.matrix(data.frame(F1 = F1.test[, i],
      #                                                                 jaccard = jaccard.test[, i],
      #                                                                 youden = youden.test[, i])))
      kendall.test[, i] <- apply(y.test, 2, function(draw) {
        # return(fastkendall(draw, y = ref)$kendall.tau)
        return(suppressWarnings(cor(draw, ref, method = "kendall")))
      })
      spearman.test[, i] <- apply(y, 2, function(draw) {
        return(cor(draw, ref, method = "spearman"))
      })
      # rharm.test[, i] <- cpp_harmonic_rowmeans(as.matrix(data.frame(kendall = kendall.test[, i],
      #                                                                spearman = spearman.test[, i])))
      ok <- kendall.test[, i] >= 0 & spearman.test[, i] >= 0
      ok[which(is.na(ok))] <- FALSE
      harm.test[ok, i] <- cpp_harmonic_rowmeans(as.matrix(data.frame(F1 = F1.test[ok, i],
                                                                   jaccard = jaccard.test[ok, i],
                                                                   youden = youden.test[ok, i],
                                                                   kendall = kendall.test[ok, i],
                                                                   spearman = spearman.test[ok, i])))
    }
  }
  # if(ndatasets>1) {
  #   F1.train[, nresults] <- cpp_harmonic_rowmeans(F1.train[, 1:nresults-1])
  #   F1.test[, nresults] <- cpp_harmonic_rowmeans(F1.test[, 1:nresults-1])
  #   jaccard.train[, nresults] <- cpp_harmonic_rowmeans(jaccard.train[, 1:nresults-1])
  #   jaccard.test[, nresults] <- cpp_harmonic_rowmeans(jaccard.test[, 1:nresults-1])
  #   youden.train[, nresults] <- cpp_harmonic_rowmeans(youden.train[, 1:nresults-1])
  #   youden.test[, nresults] <- cpp_harmonic_rowmeans(youden.test[, 1:nresults-1])
  #   clharm.train[, nresults] <- cpp_harmonic_rowmeans(clharm.train[, 1:nresults-1])
  #   clharm.test[, nresults] <- cpp_harmonic_rowmeans(clharm.test[, 1:nresults-1])
  #   kendall.train[, nresults] <- cpp_harmonic_rowmeans(kendall.train[, 1:nresults-1])
  #   kendall.test[, nresults] <- cpp_harmonic_rowmeans(kendall.test[, 1:nresults-1])
  #   spearman.train[, nresults] <- cpp_harmonic_rowmeans(spearman.train[, 1:nresults-1])
  #   spearman.test[, nresults] <- cpp_harmonic_rowmeans(spearman.test[, 1:nresults-1])
  #   rharm.train[, nresults] <- cpp_harmonic_rowmeans(rharm.train[, 1:nresults-1])
  #   rharm.test[, nresults] <- cpp_harmonic_rowmeans(rharm.test[, 1:nresults-1])
  #   harm.train[, nresults] <- cpp_harmonic_rowmeans(harm.train[, 1:nresults-1])
  #   harm.test[, nresults] <- cpp_harmonic_rowmeans(harm.test[, 1:nresults-1])
  # }
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


compute_beta_errors <- function(sims, true_beta) {
  return(sapply(sims, function (sim) {
    return(apply(sim$beta, 1, function(beta) {
      return(sqrt(mean((true_beta-beta)^2)))
    }))
  }))
}

# compute_beta_errors <- function(sims, true_beta) {
#   return(sapply(sims, function (sim) {
#     return(apply(sim$beta, 1, function(draw) {
#       sc = sqrt(mean((true_beta-draw)^2))
#       sc_draw = scale(true_beta, sc)
#       return(abs(mean(true_beta-sc_draw)))
#     }))
#   }))
# }
# 

# compute_beta_errors <- function(sims, true_beta) {
#   return(sapply(sims, function (sim) {
#     # print(str(sim$beta))
#     means <- colMeans(sim$beta)
#     return(apply(sim$beta, 1, function(beta) {
#       return(sqrt(mean((means-beta)^2)))
#     }))
#   }))
# }


reorganize <- function(trial, simulated = TRUE) {
  ntrials <- length(trial)
  ncov <- length(trial[[1]][[1]][[2]][[1]])
  ndraws <- nrow(trial[[1]][[2]][[2]][[1]])
  ngammas <- sapply(trial[[1]][[1]][[3]], length)
  metric_means <- list()
  metric_totals <- list()
  beta_means <- list()
  beta_totals <- list()
  gamma_means <- list()
  gamma_totals <- list()
  if (simulated) {
    known_beta <- matrix(0, nrow = ntrials, ncol = ncov)
    known_gamma <- list()
  } else {
    known_beta <- NA
    known_gamma <- NA
  }
  X.train <- list()
  X.test <- list()
  y.train <- list()
  y.test <- list()
  y_star.train <- list()
  y_star.test <- list()
  for (j in 1:4) {
    metric_means[[j]] <- matrix(0, nrow = ntrials, ncol = 6,
                                dimnames = list(NULL,
                                                names(trial[[1]][[1]][[1]][[1]])))
    metric_totals[[j]] <- matrix(0, nrow = ntrials*ndraws, ncol = 6,
                                 dimnames = list(NULL,
                                                 names(trial[[1]][[1]][[1]][[1]])))
    gamma_means[[j]] <- matrix(0, nrow = ntrials, ncol = ngammas[j])
    gamma_totals[[j]] <- matrix(0, nrow = ntrials*ndraws, ncol = ngammas[j])
  }
  for (j in 1:3) {
    beta_means[[j]] <- matrix(0, nrow = ntrials, ncol = ncov)
    beta_totals[[j]] <- matrix(0, nrow = ntrials*ndraws, ncol = ncov)
  }
  for (j in 1:2) {
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
  for (i in 1:ntrials) {
    if (simulated) {
      known_beta[i, ] <- trial[[i]]$data$train[[1]]$beta
    }
    for (j in 1:4) {
      metric_means[[j]][i, ] <- trial[[i]][[1]][[1]][[j]]
      #print(paste("i:", i, ",j:", j))
      #print(dim(metric_totals[[j]]))
      #print(dim(trial[[i]][[2]][[1]][[j]]))
      metric_totals[[j]][(i*ndraws - ndraws + 1):(i*ndraws), ]  <- trial[[i]][[2]][[1]][[j]]
      gamma_means[[j]][i, ] <- trial[[i]][[1]][[3]][[j]]
      gamma_totals[[j]][(i*ndraws - ndraws + 1):(i*ndraws), ] <- trial[[i]][[2]][[3]][[j]]
      
    }
    for (j in 1:3) {
      beta_means[[j]][i, ] <- trial[[i]][[1]][[2]][[j]]
      beta_totals[[j]][(i*ndraws - ndraws + 1):(i*ndraws), ] <- trial[[i]][[2]][[2]][[j]]
    }
    for (j in 1:2) {
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
