source("R/data.R")
source("R/fit.R")
source("R/plot.R")

library(doParallel)
library(foreach)
library(coda)

# Generate samples for the experiment.
generate_experiment_samples <- function(data) {

    # Set up parallelization workers.
    cl<-makeForkCluster(5) # 5 workers.
    on.exit(stopCluster(cl))
    registerDoParallel(cl)
    clusterSetRNGStream(cl, 42)

    # Generate samples for 500 trials.
    res <- foreach(i = 1:500) %dopar% {
        trial_res <- generate_trial_samples(
            data = data,
            train_split_prop = 0.666, # 2/3
            burnin = 50000,
            ndraws = 50000,
            thin = 10,
            tune_iterations = 10000
        )
        return(trial_res)
    }

    return(res)
}

# Generate samples for a single trial of the experiment.
generate_trial_samples <- function(
    data,
    train_split_prop,
    burnin,
    ndraws,
    thin,
    tune_iterations
) {
    # Split data into train and test sets.
    splits <- split_data(data, prop = train_split_prop)
    train_data <- splits$train
    test_data <- splits$test

    # Tune the sampler.
    tune_results <- tune_mspm(
        data = train_data,
        iterations = tune_iterations
    )

    # Run burn-in.
    burnin_fit <- fit_mspm(
        data = train_data,
        ndraws = burnin,
        burnin = 0, # Set 0 burnin because we don't want to discard any samples.
        thin = thin,
        proposal_variance = get_proposal_variance(tune_results),
        compute_diagnostics = FALSE
    )

    # Set the starting values for the proper sampling to the last values of the burn-in chain.
    beta_start = as.numeric(tail(as.matrix(burnin_fit$beta), 1))
    gamma_start <- lapply(burnin_fit$gammas, function(gamma_chain) {
        as.numeric(tail(as.matrix(gamma_chain), 1))
    })

    # Run proper sampling.
    sampled_fit <- fit_mspm(
        data = train_data,
        ndraws = ndraws,
        burnin = 0,
        thin = thin,
        proposal_variance = get_proposal_variance(tune_results),
        beta_start = beta_start,
        gamma_start = gamma_start,
        compute_diagnostics = FALSE
    )

    return(list(
        burnin_fit = burnin_fit,
        sampled_fit = sampled_fit,
        estimated_proposal_variance = get_proposal_variance(tune_results),

        burnin_time = burnin_fit$samplingTime,
        sampled_time = sampled_fit$samplingTime,

        burnin_nlikelihood_calls = burnin_fit$nlikelihood_calls,
        sampled_nlikelihood_calls = sampled_fit$nlikelihood_calls

    ))
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
experiment_samples <- generate_experiment_samples(data)
saveRDS(experiment_samples, "experiments/results/gibbs_samples.rds")

print("Done")