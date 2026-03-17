
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
fit_pt <- fit_mspm_pt(
    data = data,
    ndraws = 10000,
    burnin = 10000,
    thin = 10,
    proposal_variance = 0.1,
    seed = 1234,
    ntemperatures = 10,
    complete_param_swapping = FALSE,
    verbose = 1000
)

# Plot fit.
plot_posteriors_beta(fit)
plot_posterior_gammas(fit)

# Plot chains.
plot_beta_chains(fit)
plot_gamma_chains(fit, seperateGraphs = TRUE)
plot_gamma_chains(fit)