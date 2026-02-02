# This script generates a test fit of the multi-scale probit model using synthetic data.
# The fitted model's posterior distributions for beta coefficients and gamma thresholds
# are then plotted.

source("R/fit.R")
source("R/util.R")
source("R/plot.R")
source("R/data.R")


# Initial parameter tuning code
burnin <- 2000
ndraws <- 2000
thin <- 10
verbose <- 0

# Generate the data.
mspm.data <- generate_synthetic_data(
    nobs = 60,
    ncov = 48,
    ngamma = c(1, 3, 3),
    seed = 1234
)

# Fit using all the data.
tst_fit = mspm.fit <- fit_mspm(
    data = mspm.data,
    ndraws = ndraws,
    burnin = burnin,
    thin = thin,
    tune = 0.1,
    seed = 1234,
    verbose = verbose
)
tst_fit

# Plot fit.
plot_posteriors_beta(tst_fit)
plot_posterior_gammas(tst_fit)
