# Implementation of the Multi-Scale Probit

The Multi-Scale Probit (MSP) is an extension of the Bayesian Ordered Probit. It is suitable for cross-source fitting on different ordinal scales. 

## Use Case: Modeling Text Complexity

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

## Description
Let people know what your project can do specifically. Provide context and add a link to any reference visitors might be unfamiliar with. A list of Features or a Background subsection can also be added here. If there are alternatives to your project, this is a good place to list differentiating factors.

## Badges
On some READMEs, you may see small images that convey metadata, such as whether or not all the tests are passing for the project. You can use Shields to add some to your README. Many services also have instructions for adding a badge.

## Visuals
Depending on what you are making, it can be a good idea to include screenshots or even a video (you'll frequently see GIFs rather than actual videos). Tools like ttygif can help, but check out Asciinema for a more sophisticated method.

## Installation
Within a particular ecosystem, there may be a common way of installing things, such as using Yarn, NuGet, or Homebrew. However, consider the possibility that whoever is reading your README is a novice and would like more guidance. Listing specific steps helps remove ambiguity and gets people to using your project as quickly as possible. If it only runs in a specific context like a particular programming language version or operating system or has dependencies that have to be installed manually, also add a Requirements subsection.

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

## Support
Tell people where they can go to for help. It can be any combination of an issue tracker, a chat room, an email address, etc.

## Roadmap
If you have ideas for releases in the future, it is a good idea to list them in the README.

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
For open source projects, say how it is licensed.

## Project status
If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.
