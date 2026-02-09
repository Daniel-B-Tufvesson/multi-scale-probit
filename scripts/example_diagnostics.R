source("R/data.R")
source("R/fit.R")

library(coda)

# Generate the data.
data <- generate_synthetic_data(
    nobs = 60,
    ncov = 6,
    ngamma = c(1, 3, 3),
    seed = 1234
)

# Fit the model.
fit1 <- fit_mspm(
    data = data,
    ndraws = 2000,
    burnin = 2000,
    thin = 10,
    tune = 0.1,
    seed = 123,
)
fit2 <- fit_mspm(
    data = data,
    ndraws = 2000,
    burnin = 2000,
    thin = 10,
    tune = 0.1,
    seed = 42,
)

# Compute Geweke diagnostic for beta coefficients.
geweke_result <- geweke.diag(beta(fit1))
geweke_result

# Compute the Gelman-Rubin R-hat for beta coefficients.
rhat_beta <- gelman.diag(mcmc.list(beta(fit1), beta(fit2)))
rhat_beta

#rhat_gammas <- gelman.diag(mcmc.list(gammas(fit1), gammas(fit2)))
#rhat_gammas