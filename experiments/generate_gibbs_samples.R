source("R/data.R")
source("R/fit.R")

library(doParallel)
library(foreach)

# Generate samples for the experiment.
generate_experiment_samples <- function(data) {

    # Set up parallelization workers.
    cl<-makeForkCluster(7) # 7 workers.
    on.exit(stopCluster(cl))
    registerDoParallel(cl)
    clusterSetRNGStream(cl, 42)

    proposal_variance <- tune_gibbs(data, 10000)

    # Generate samples for 500 trials.
    res <- foreach(i = 1:500) %dopar% {
        samples <- generate_gibbs_samples(
            data = data,
            ndraws = 100000,
            thin = 10,
            proposal_variance = proposal_variance
        )
        return(samples)
    }

    return(res)
}

tune_gibbs <- function(data, iterations) {
    tune_results <- tune_mspm(
        data = data,
        iterations = iterations
    )
    return (get_proposal_variance(tune_results))
}

generate_gibbs_samples <- function(
    data,
    ndraws,
    thin,
    proposal_variance
) {
    sampled_fit <- fit_mspm(
        data = data,
        ndraws = ndraws,
        burnin = 0,
        thin = thin,
        proposal_variance = proposal_variance,
        compute_diagnostics = FALSE
    )
    return(sampled_fit)
}

# Execute script --------------------------------------------------------------------------

# Generate the data.
data <- generate_synthetic_data(
    nobs = 400,
    ncov = 48,
    ngamma = c(1, 3, 3),
    seed = 42
)

# Generate samples.
gibbs_runs <- generate_experiment_samples(data)
saveRDS(gibbs_runs, "experiments/results/gibbs-samples-big.rds")