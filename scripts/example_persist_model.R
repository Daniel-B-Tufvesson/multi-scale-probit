source("R/persistence.R")
source("R/data.R")
source("R/fit.R")

# Generate synthetic data and fit model.
data <- generate_synthetic_data(
    nobs = 100,
    ncov = 5,
    ngamma = c(2, 3),
    seed = 42
)

# Fit model.
fit = fit_mspm(
    data = mspm.data,
    ndraws = 100,
    burnin = 100,
    thin = 2,
    tune = 0.1,
    seed = 1234,
    verbose = 0
)
fit

# Save model to temporary file.
temp_file <- tempfile(fileext = ".rds")
save_model(fit, temp_file)

# Load model from file.
loaded_fit <- load_model(temp_file)
loaded_fit