# An initial cross validation experiment with synthetic data and the standard Gibbs sampler.

source("R/data.R")
source("R/fit.R")

data <- generate_synthetic_data(
    nobs = 500,
    ncov = 48,
    ngamma = c(1, 5, 4),
    seed = 1234
)

start_time <- Sys.time()
cv_res <- cross_validate(
    data = data,
    nsplits = 500,
    prop = 0.66,
    ndraws = 50000,
    burnin = 50000,
    thin = 100,
    seed = 1234,
    nworkers = 7, # 8 cores on M2 macbook air
    meansOnly = FALSE
)
end_time <- Sys.time()
print(end_time - start_time)
# Time difference of 1.597387 hours

#saveRDS(cv_res, file = "synth_gibbs_cv")
