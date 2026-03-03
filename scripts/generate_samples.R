# Script that generates a bunch of chains for later analysis. Saved to json.
library(jsonlite)

source("R/data.R")
source("R/fit.R")

# Generate synthetic data.
data <- generate_synthetic_data(
    nobs = 100,
    ncov = 10,
    ngamma = c(1, 3, 3),
    seed = 1234
)

# Fit several models to the data with different seeds.
nchains <- 4
results <- list(nchains = nchains, chains = list(), type = "gibbsmh")
for (i in 1:nchains) {
    seed_i <- 1234 + i
    fit <- fit_mspm(
        data = data,
        ndraws = 100000,
        burnin = 100000,
        thin = 100,
        seed = seed_i,
        saveBurninSamples = TRUE
    )
    results$chains[[i]] <- list(
        seed = seed_i,
        samplingTime = fit$samplingTime,
        burninTime = fit$burninTime,
        beta = as.matrix(fit$beta),
        gammas = lapply(fit$gammas, as.matrix),
        burninBeta = as.matrix(fit$burninBeta),
        burninGammas = lapply(fit$burninGammas, as.matrix)

    )
}

# Save samples to JSON
json_data <- toJSON(results, pretty = TRUE, auto_unbox = TRUE)
write(json_data, file = "gibbs_march01.json")