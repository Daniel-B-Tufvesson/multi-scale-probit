source("R/fit.R")
source("R/util.R")
source("R/plot.R")
source("R/data.R")
source("R/predict.R")
source("R/eval.R")

# Generate the data.
data <- generate_synthetic_data(
    nobs = 100,
    ncov = 48,
    ngamma = c(1, 3, 3),
    seed = 1234
)

# Split data.
split_data = split_data(
    data = data,
    prop = 0.8,
    seed = 1234
)
train = split_data$train
test = split_data$test

# Fit using the training data.
fit <- fit_mspm(
    data = train,
    ndraws = 2000,
    burnin = 2000,
    thin = 10,
    tune = 0.1,
    seed = 1234
)

# Predict latent ystar values for training data.
latent.train <- predict_mspm(fit, newdata = train, latentOnly = TRUE)

# Predict labels for training data.
labels.train <- predict_mspm(fit, newdata = train)

# Predict labels for test data.
labels.test <- predict_mspm(fit, newdata = test)

# Evaluate predictions.
train.eval <- eval_mspm_prediction_draws(predictions = labels.train, test_data = train)
test.eval <- eval_mspm_prediction_draws(predictions = labels.test, test_data = test)

# Plot performance metrics for training data.
plot_eval_draws(
    eval = list(train.eval, test.eval),
    plotMean = TRUE
)