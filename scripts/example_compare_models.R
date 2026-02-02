# Compare the predictive performance of two models.

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
split_data = split_data.mspm_data(
    data = data,
    prop = 0.8,
    seed = 1234
)
train = split_data$train
test = split_data$test

# Fit two models using the training data.
fit1 <- fit_mspm(
    data = train,
    ndraws = 2000,
    burnin = 2000,
    thin = 10,
    tune = 0.1,
    seed = 1234,
    verbose = 0
)
fit2 <- fit_mspm(
    data = train,
    ndraws = 200, # Fewer draws
    burnin = burnin,
    thin = 1, # No thinning
    tune = 0.1,
    seed = 1234,
    verbose = 0
)

# Predict labels for test data.
pred1 <- predict_mspm(
  fit = fit1,
  newdata = test
)
pred2 <- predict_mspm(
  fit = fit2,
  newdata = test
)


# Evaluate predictions.
eval1 <- eval_mspm_prediction_draws(
    predictions = pred1
)
eval2 <- eval_mspm_prediction_draws(
    predictions = pred2
)


# Plot performance metrics.
# plot_eval_draws(
#     eval = list(eval1, eval2),
#     plotMean = TRUE
# )

# Plot diff in metrics.
plot_eval_draws_diff(
    eval1 = eval1,
    eval2 = eval2,
    label1 = "Thinning",
    label2 = "No thinning"
)