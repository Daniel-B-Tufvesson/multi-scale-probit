# Implementation of the Multi-Scale Probit

The Multi-Scale Probit (MSP) is an extension of the Bayesian Ordered Probit. It is suitable for cross-source fitting on different ordinal scales. 

## Mathematical Definition
MSP models the relationship between the continuous input features $\bf{x}$ and an ordinal response variable $y$ via a latent variable $`y^*`$. For a given observation $i$, the value $`y^*_i`$ is defined as 

$$y^*_i=\bf{x}^T_i \bf{\beta}+\epsilon$$

where $\bf{x}_i$ is the observed vector of the features, $\bf{\beta}$ is a vector of coefficients corresponding to $\bf{x}$ and $\epsilon \sim N(0,1)$ which is independent noise. 

What makes MSP different from the standard Bayesian Ordered Probit is that it can model $`y*`$ over multiple datasets with different ordinal scales $s$ for $y$. Let $s_i$ be the scale for observation $i$. Then predicting $y_i$ for a given $`y_i^*`$ is defined as

$$
y_i = 
\begin{cases}
    1 & \text{if } y_i^* \le \gamma_1^{(s_i)}, \\
    2 & \text{if } \gamma_1^{(s_i)} \le y_i^* \le \gamma_2^{(s_i)}, \\
    \vdots \\
    C^{(s_i)} & \text{if } \gamma_{C^(s_i)-1}^{(s_i)} <  y_i^*,
\end{cases}
$$

where $s_i$ indicates the ordinal scale of $y_i$, $C^{(s_i)}$ is the number of ordinal values for the scale $s_i$, and $\gamma^{(s_i)}$ are the thresholds between each value on the scale $s_i$.

(Footnote: this section is taken directly from my thesis.) 

## Usage
```R
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

# Predict labels for training data.
labels.train <- predict_mspm(fit, newdata = train)

# Predict labels for test data.
labels.test <- predict_mspm(fit, newdata = test)

# Evaluate predictions.
train.eval <- eval_mspm_prediction_draws(predictions = labels.train, test_data = train)
test.eval <- eval_mspm_prediction_draws(predictions = labels.test, test_data = test)

# Plot performance metrics for training data as distributions.
plot_eval_draws(
    eval = list(train.eval, test.eval),
    plotMean = TRUE
)
```

