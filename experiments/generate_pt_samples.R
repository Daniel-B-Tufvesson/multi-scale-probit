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

    ntemperatures <- 10

    print("Tune sampler.")
    tune_results <- tune_pt(data, 10000, ntemperatures)

    # Generate samples for 500 trials.
    print("Start parallel samplers")
    res <- foreach(i = 1:500) %dopar% {
        samples <- generate_pt_samples(
            data = data,
            ndraws = 100000,
            thin = 10,
            ntemperatures = ntemperatures,
            inv_temperature_ladder = get_inv_temperatures(tune_results),
            proposal_variance = get_proposal_variance(tune_results),
            complete_param_swapping = TRUE
        )
        return(samples)
    }

    return(res)
}

tune_pt <- function(data, iterations, ntemperatures) {
    tune_results <- tune_mspm_pt(
        data = data,
        iterations = iterations,
        ntemperatures = ntemperatures
    )
    return (tune_results)
}

generate_pt_samples <- function(
    data,
    ndraws,
    thin,
    ntemperatures,
    inv_temperature_ladder,
    proposal_variance,
    complete_param_swapping
) {
    # Start params at random state.
    beta_start <- rnorm(48, sd = 4)
    gamma_start <- list()
    for (i in 1:ntargets(data)) {
        gamma_start[[i]] <- sort(rnorm(nlevels(data)[i], sd = 4))
    } 

    sampled_fit <- fit_mspm_pt(
        data = data,
        ndraws = ndraws,
        burnin = 0,
        thin = thin,
        ntemperatures = ntemperatures,
        inv_temperature_ladder = inv_temperature_ladder,
        proposal_variance = proposal_variance,
        compute_diagnostics = FALSE,
        beta_start = beta_start,
        gamma_start = gamma_start,
        complete_param_swapping = complete_param_swapping
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
pt_runs <- generate_experiment_samples(data)
saveRDS(pt_runs, "experiments/results/pt-samples-10-chains-complete-big.rds")
print("Done")