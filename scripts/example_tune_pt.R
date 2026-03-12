source("R/fit.R")
source("R/data.R")


# Generate the data.
data <- generate_synthetic_data(
    nobs = 400 * 3,
    ncov = 48,
    ngamma = c(1, 3, 3),
    seed = 1234
)

# Find the optimal tuning parameters for the sampler.
tune_results <- tune_mspm_pt(
    data = data,
    iterations = 50000,
    ntemperatures = 5,
    verbose = 1000,
    target_epsilon = 0.01,
    stop_early = TRUE
)

# Fit using all the data and the tuned parameters.
fit <- fit_mspm_pt(
    data = data,
    ndraws = 2000,
    burnin = 2000,
    thin = 10,
    ntemperatures = 5,
    inv_temperature_ladder = tune_results$inv_temperature_ladder,
    proposal_variance = tune_results$proposal_variance,
    seed = 1234,
    verbose = 100
)
#fit
