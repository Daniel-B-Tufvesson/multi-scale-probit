# Baseline experiment.
# 1. Run 500 chains for 50000 iterations + 50000 burn-in iterations. Adaptive burn-in is done for each run.
# 2. Compute R-hat over all chains for each iteration and record when it reaches < 1.05. This is done on the burn-in chain.
# 3. Record ESS per second and iteration for post-burn-in samples.
# 4. Verify and count each chain that is converged during post-burn-in. This is only for reporting.

source("R/data.R")
source("R/fit.R")
source("R/cv.R")
source("R/plot.R")

# Generate the data.
data <- generate_synthetic_data(
    nobs = 400 * 3, # 400 per target.
    ncov = 48,
    ngamma = c(1, 3, 3),
    seed = 1234
)

cv_res = cross_validate(
    data = data,
    nsplits = 5,
    ndraws = 5000,
    prop = 0.66,
    sampler = fit_mspm,
    samplerArgs = list(
        burnin = 5000,
        thin = 1,
        saveBurninSamples = TRUE,
        adapt_tune = TRUE
    ),
    seed = 42,
    nworkers = 5,
    meansOnly = FALSE,
    saveSamples = TRUE
)

saveRDS(cv_res, file = "results/cv_gibbs_5000.rds")
print("Saved results")

cv_loaded = readRDS("results/cv_gibbs_5000.rds")
print("Loaded results")
