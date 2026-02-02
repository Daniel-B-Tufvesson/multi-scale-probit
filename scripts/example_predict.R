source("R/fit.R")
source("R/util.R")
source("R/plot.R")
source("R/data.R")
source("R/predict.R")
source("R/eval.R")

# Generate the data.
mspm.data <- generate_synthetic_data(
    nobs = 100,
    ncov = 48,
    ngamma = c(1, 3, 3),
    seed = 1234
)

# Split data.
split_data = split_data.mspm_data(
    data = mspm.data,
    prop = 0.8,
    seed = 1234
)
mspm.train = split_data$train
mspm.test = split_data$test

# Initial parameter tuning code
burnin <- 2000
ndraws <- 2000
thin <- 10
verbose <- 0

# Fit using the training data.
mspm.fit <- fit_mspm(
    data = mspm.train,
    ndraws = ndraws,
    burnin = burnin,
    thin = thin,
    tune = 0.1,
    seed = 1234,
    verbose = verbose
)

# Predict latent ystar values for training data.
mspm.latent.train <- predict_mspm(
    fit = mspm.fit,
    latentOnly = TRUE
)

# Predict labels for training data.
mspm.labels.train <- predict_mspm(
    fit = mspm.fit
)

# Predict labels for test data.
mspm.labels.test <- predict_mspm(
  fit = mspm.fit,
  newdata = mspm.test
)


# Evaluate predictions.
train.eval <- eval_mspm_prediction_draws(
    predictions = mspm.labels.train
)

test.eval <- eval_mspm_prediction_draws(
    predictions = mspm.labels.test
)

# Plot performance metrics for training data.
plot_eval_draws(
    eval = list(train.eval, test.eval),
    plotMean = TRUE
)