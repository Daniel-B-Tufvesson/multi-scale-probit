library(MASS)
library(coda)
#library(Rcpp)
library(beepr)
#.libPaths("/home/elsan830/R/x86_64-pc-linux-gnu-library/3.6")
#sourceCpp("src/hprobit_gibbs.cpp")

# Compute the log-likelihood for an ordered probit model. 
# Returns the log-likelihood value.
# 
# Arguments:
# gamma: A vector of threshold parameters.
# y: A vector of observed responses.
# Xbeta: A matrix of linear predictors (X * beta).
# ncat: The number of categories in the response variable.
llgamma <- function(gamma, y, Xbeta, ncat) {
    # Note: this function appears to be unused.
    print("---")
    print(gamma)
    print(length(y))
    print(dim(Xbeta))
    print(ncat)
    ll = 0.0
    for (i in nrow(Xbeta)){
        if (y[i] == ncat){
        ll = ll + log(1.0  - pnorm(gamma[y[i]-1] - Xbeta[i]))
        }
        else if (y[i] == 1){
        ll = ll + log(pnorm(gamma[y[i]] - Xbeta[i]))
        }
        else{
        ll = ll + log(pnorm(gamma[y[i]] - Xbeta[i]) -
                        pnorm(gamma[y[i]-1] - Xbeta[i]))
        }
    }
    return(ll)
}

# Fit a multi-scale probit model using Gibbs sampling.
# Returns a list containing posterior draws and model information.
# 
# Arguments:
# dat: A list of data frames, each containing 'X' (predictors) and 'Y' (responses).
# burnin: Number of burn-in iterations.
# ndraws: Number of posterior draws to collect.
# thin: Thinning interval for MCMC sampling.
# meanPrior: Prior mean for regression coefficients.
# precPrior: Prior precision for regression coefficients.
# fix.zero: Index of the threshold to fix at zero.
# tune: Tuning parameter for the sampler.
# seed: Random seed for reproducibility.
# beta.start: Initial values for regression coefficients.
# gamma.start: Initial values for threshold parameters.
# verbose: Verbosity level for output.
hprobit_gibbs <- function(dat,
                          burnin,
                          ndraws,
                          thin,
                          meanPrior,
                          precPrior,
                          fix.zero = 1,
                          tune = NULL,
                          seed = NA,
                          beta.start = NULL,
                          gamma.start = NULL,
                          verbose = 0) {
  # Set the random seed for reproducibility
  if (is.na(seed)) {
    seed <- round(runif(1, 0, .Machine$integer.max))
    set.seed(round(runif(1, 0, .Machine$integer.max)))
  } else {
    set.seed(seed)
  }

  # Initialize lists and extract data dimensions
  X.names = list()
  ncat <- c()
  N <- c()
  K <- c()
  ntargets <- length(dat)
  Xlist <- list()
  ylist <- list()
  for (i in 1:ntargets) {
    X.names[[i]] <- dimnames(dat[[i]]$X)[[2]]
    ncat[i] <- length(unique(dat[[i]]$Y))
    N[i] <- nrow(dat[[i]]$X)
    K[i] <- ncol(dat[[i]]$X)
    Xlist[[i]] = dat[[i]]$X
    ylist[[i]] = dat[[i]]$Y+1
  }

  # Check that all covariate matrices have the same columns
  if (length(unique(X.names)) != 1) {
    stop("All covariate matrices needs to have the same columns. \n")
  } else {
    K <- K[1]
  }

  # Set starting values for gamma and beta if not provided
  if (is.null(gamma.start)) {
    gamma.start <- list()
    for (i in 1:ntargets) {
      gamma.start[[i]] <- rep(0, ncat[i]+1)
      gamma.start[[i]][1] <- -300
      gamma.start[[i]][2:ncat[i]] <- 0:(ncat[i]-2)
      gamma.start[[i]][ncat[i] + 1] <- 300
    }
  }
  if (is.null(beta.start)) {
    beta.start <- rep(0, K)
  }

  # Set tuning parameter for the sampler
  if (is.null(tune)) {
    tune <- 0.05/ncat
  } else if (length(tune) != ntargets) {
    tune <- rep(tune, ntargets)
  }

  # Run the main Gibbs sampler (C++ implementation)
  sim = cpp_hprobit(Xlist,
                    ylist,
                    meanPrior,
                    precPrior,
                    fix.zero,
                    ncat,
                    gamma.start,
                    beta.start,
                    tune,
                    ndraws,
                    burnin,
                    thin,
                    seed,
                    verbose)

  # Format and store MCMC results
  colnames(sim$storebeta) <- c(X.names[[1]])
  beta <- mcmc(sim$storebeta, start = burnin + 1, thin = thin)
  gammas <- list()
  for (i in 1:ntargets) {
    gammas[[i]] <- mcmc(sim$storegamma[[i]], start = burnin + 1,
                        thin = thin)
  }

  # Print acceptance rate and return results
  print(paste('Acceptance rate:', paste(sim$accepts/ndraws, collapse = " "), sep = " "))
  #beep(4)
  return(list(beta = beta,
              gammas = gammas,
              X = Xlist,
              y = ylist
        ))
}


  