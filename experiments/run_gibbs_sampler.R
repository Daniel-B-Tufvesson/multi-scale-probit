source("R/data.R")
source("R/fit.R")
source("R/plot.R")
source("R/predict.R")
source("R/eval.R")

library(doParallel)
library(foreach)
library(coda)

# Generate samples for the experiment.
run_gibbs <- function(data) {

    # Set up parallelization workers.
    cl<-makeForkCluster(5) # 5 workers.
    on.exit(stopCluster(cl))
    registerDoParallel(cl)
    clusterSetRNGStream(cl, 42)

    # Generate samples for 500 trials.
    res <- foreach(i = 1:500) %dopar% {
        trial_res <- run_gibbs_trial(
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
run_gibbs_trial <- function(
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

    # Predict labels for test data.
    predicted_labels <- predict_mspm(sampled_fit, newdata = test_data)

    # Evaluate the predictions.
    evaluation <- eval_mspm_prediction_draws(predictions = predicted_labels, 
        test_data = test_data)

    return(list(
        burnin_fit = burnin_fit,
        sampled_fit = sampled_fit,
        estimated_proposal_variance = get_proposal_variance(tune_results),
        evaluation = evaluation
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
gibbs_runs <- run_gibbs(data)
saveRDS(gibbs_runs, "experiments/results/gibbs_runs2.rds")

loaded_runs = readRDS("experiments/results/gibbs_runs2.rds")

print("Done")