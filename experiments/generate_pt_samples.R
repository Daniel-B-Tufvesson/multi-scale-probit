source("R/data.R")
source("R/fit.R")

library(doParallel)
library(foreach)

# Generate samples for the experiment.
generate_experiment_samples <- function(
    data, 
    ntemperatures, 
    ndraws, 
    nsplits, 
    tune_iterations,
    complete_param_swapping
) {

    # Set up parallelization workers.
    ncores <- parallel::detectCores() - 1
    cat("Setting up parallelization with ", ncores, " cores.\n", sep = "")

    cl<-makeForkCluster(ncores) # 7 workers.
    on.exit(stopCluster(cl))
    registerDoParallel(cl)
    clusterSetRNGStream(cl, 42)

    ntemperatures <- 3

    print("Tune sampler.")
    tune_results <- tune_pt(data, tune_iterations, ntemperatures)

    # Generate samples for 500 trials.
    print("Start parallel samplers")
    res <- foreach(i = 1:nsplits) %dopar% {
        samples <- generate_pt_samples(
            data = data,
            ndraws = ndraws,
            thin = 10,
            ntemperatures = ntemperatures,
            inv_temperature_ladder = get_inv_temperatures(tune_results),
            proposal_variance = get_proposal_variance(tune_results),
            complete_param_swapping = complete_param_swapping
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

execute_experiment <- function(
    ntemperatures,
    output_directory = "experiments/results",
    nsplits = 500,
    ndraws = 100000,
    tune_iterations = 10000,
    complete_param_swapping = TRUE
) {
    # Generate the data.
    data <- generate_synthetic_data(
        nobs = 400,
        ncov = 48,
        ngamma = c(1, 3, 3),
        seed = 42
    )

    # Generate samples.
    pt_runs <- generate_experiment_samples(
        data,
        ntemperatures = ntemperatures,
        ndraws = ndraws,
        nsplits = nsplits,
        tune_iterations = tune_iterations,
        complete_param_swapping = complete_param_swapping
    )
    complete_swapping <- ifelse(complete_param_swapping, "complete", "partial")
    path <- file.path(output_directory, paste0("pt-samples-", ntemperatures, "-temperatures-", nsplits, "-splits-", complete_swapping, "-", ndraws,"-draws.rds"))
    saveRDS(pt_runs, path)
    print("Done")
}


# Run dummy experiment.
# execute_experiment(
#     ntemperatures = 3,
#     output_directory = "experiments/results",
#     nsplits = 3,
#     ndraws = 1000,
#     tune_iterations = 100
# )