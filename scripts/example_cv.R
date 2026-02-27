# This script demonstrates how to perform cross-validation on an mspm model using synthetic data.

source("R/data.R")
source("R/fit.R")
source("R/cv.R")
source("R/plot.R")


# Generate the data.
data <- generate_synthetic_data(
    nobs = 60,
    ncov = 6,
    ngamma = c(1, 3, 3),
    seed = 1234
)

# Perform cross-validation on the standard mh-within-gibbs sampler.
cv_res1 <- cross_validate(
    data = data,
    nsplits = 50,
    prop = 0.66,
    ndraws = 50000,
    sampler = fit_mspm,
    samplerArgs = list(
        burnin = 50000,
        thin = 10
    ),
    seed = 42,
    nworkers = 5,
    meansOnly = FALSE
)

# Perform cross-validation on the parallel tempering sampler.
cv_res_pt <- cross_validate(
    data = data,
    nsplits = 50,
    prop = 0.66,
    ndraws = 50000,
    sampler = fit_mspm_pt,
    samplerArgs = list(
        burnin = 50000,
        thin = 10,
        ntemperatures = 5
    ),
    seed = 42,
    nworkers = 5,
    meansOnly = FALSE
)

print("done")


# Plot difference.
plot_cv_diff(cv_res1, cv_res_pt, title = "Means per split")
# cvAllDraws(cv_res1)
plot_cv_diff(cv_res1, cv_res_pt, plotData = "allDraws", title = "Total draws")

