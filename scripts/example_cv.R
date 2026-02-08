# This script demonstrates how to perform cross-validation on an mspm model using synthetic data.

source("R/data.R")
source("R/cv.R")


# Generate the data.
data <- generate_synthetic_data(
    nobs = 60,
    ncov = 6,
    ngamma = c(1, 3, 3),
    seed = 1234
)

# Perform cross-validation.
cv_res <- cross_validate(
    data = data,
    nsplits = 10,
    prop = 0.7,
    ndraws = 100,
    burnin = 100,
    thin = 1,
    seed = 42,
    nworkers = 10
)
cv_res
