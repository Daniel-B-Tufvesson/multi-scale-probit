# This script generates a test fit of the multi-scale probit model using synthetic data.
# The fitted model's posterior distributions for beta coefficients and gamma thresholds
# are then plotted.

source("R/fit.R")
source("R/util.R")
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
fit <- fit_mspm(
    data = data,
    ndraws = 2000,
    burnin = 2000,
    thin = 10,
    tune = 0.1,
    seed = 1234,
)
fit

# Plot fit.
#plot_posteriors_beta(fit)
#plot_posterior_gammas(fit)

# Plot chains.
plot_beta_chains(fit)
plot_gamma_chains(fit, seperateGraphs = TRUE)
plot_gamma_chains(fit)