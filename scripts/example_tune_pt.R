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
    iterations = 10000,
    ntemperatures = 5,
    verbose = 1000,
    target_epsilon = 0.01,
    stop_early = TRUE,
    tune_proposal_variance = TRUE,
    temperature_window_size = 100,
    temperature_ladder_learning_rate = 0.1,
    temperature_window_growth_factor = 1.3
)

# Fit using all the data and the tuned parameters.
fit <- fit_mspm_pt(
    data = data,
    ndraws = 2000,
    burnin = 2000,
    thin = 10,
    ntemperatures = 5,
    inv_temperature_ladder = get_inv_temperatures(tune_results),
    proposal_variance = get_proposal_variance(tune_results),
    seed = 1234,
    verbose = 100
)
#fit
