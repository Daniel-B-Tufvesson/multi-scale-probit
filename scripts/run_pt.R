
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
    ndraws = 2000,
    burnin = 2000,
    thin = 10,
    tune = 0.1,
    seed = 1234,
    ntemperatures = 3
)
fit