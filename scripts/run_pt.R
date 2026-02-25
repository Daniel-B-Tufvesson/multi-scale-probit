
source("R/fit.R")
source("R/plot.R")
source("R/data.R")


# Generate the data.
data <- generate_synthetic_data(
    nobs = 60,
    ncov = 6,
    ngamma = c(1, 3, 3),
    seed = 1234
)

# Fit using all the data.
fit <- fit_mspm_pt(
    data = data,
    ndraws = 500,
    burnin = 10000,
    thin = 2,
    tune = 0.1,
    seed = 1234,
    ntemperatures = 5,
    verbose = 100,
    temperature_ladder_learning_rate = 0.01,
    tempperature_window_growth_factor = 2,
    complete_param_swapping = TRUE
)
#fit