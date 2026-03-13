/**
 * @file hprobit_gibbs.cpp
 * @author Johan Falkenjack
 * @author Daniel Tufvesson
 * @version 4.0
 * @brief A Metropolis-Hastings-withing-Gibbs sampler for the Multi-Scale Probit model.
 */

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "sampling.hpp"
#include <chrono>


/**
 * Do the burn-in phase of the MCMC sampler, which runs the sampler for a specified number of 
 * iterations without storing the samples.
 * 
 * @param chain The MspmChain object representing the current state of the MCMC chain.
 * @param data The Data object containing the feature matrices and response vectors for each target.
 * @param burnin The number of burn-in iterations to perform.
 * @param gen The GSL random number generator to use for sampling.
 * @param verbose The verbosity level for printing progress during burn-in. A value of 0 means
 * no progress will be printed, while higher values will print progress every `verbose` iterations.
 * 
 * @return The time taken for the burn-in phase in seconds.
 */
double do_burnin(
    MspmChain& chain,
    const Data& data,
    int burnin,
    gsl_rng* gen,
    int verbose
) {
    // Measure burnin sampling time.
    auto start_time_burnin = std::chrono::high_resolution_clock::now();

    // Burnin loop.
    for (int iter = 0; iter < burnin; iter++) {
        chain.simulate_step(data, gen);

        // Print progress.
        if(verbose > 0 && iter % verbose == 0){
            Rcpp::Rcout << "probit_gibbs burnin iteration " << (iter+1) 
                << " of " << burnin << std::endl;

            // Print acceptance rates for each target.
            for (int target = 0; target < data.ntargets; ++target) {
                double acceptance_rate = chain.cumulative_acceptance_probabilities(target) / static_cast<double>(iter+1);
                Rcpp::Rcout << "Metropolis acceptance rate for gamma (target " << target << ") = " 
                    << acceptance_rate << std::endl;
            }
        }

    }
    // Measure burnin time in seconds.
    auto end_time_burnin = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_burnin = end_time_burnin - start_time_burnin;
    return elapsed_burnin.count();
}

/**
 * Perform the sampling phase of the MCMC sampler, after burn-in is complete. This runs the MCMC
 * sampler for a specified number of iterations, and stores the sampled beta and gamma values.
 * 
 * @param iterations The number of MCMC iterations to perform during the sampling phase.
 * @param chain The MspmChain object representing the current state of the MCMC chain.
 * @param data The Data object containing the feature matrices and response vectors for each target.
 * @param sample_storage The SampleStorage object where the sampled beta and gamma values will be 
 * stored.
 * @param thin The thinning interval for storing samples.
 * @param gen The GSL random number generator to use for sampling.
 * @param verbose The verbosity level for printing progress during sampling. A value of 0 means
 * no progress will be printed, while higher values will print progress every `verbose` iterations.
 * 
 * @return The time taken for the sampling phase in seconds.
 */
double do_sampling(
    int iterations,
    MspmChain& chain,
    const Data& data,
    SampleStorage& sample_storage,
    const int thin,
    gsl_rng* gen,
    int verbose
) {
    // Measure sampling time.
    auto start_time = std::chrono::high_resolution_clock::now();

    // Sampling loop.
    for (int iter = 0; iter < iterations; iter++) {
        chain.simulate_step(data, gen);
        
        // Print progress.
        if(verbose > 0 && iter % verbose == 0){
            Rcpp::Rcout << "probit_gibbs iteration " << (iter+1) 
                << " of " << iterations << std::endl;

            // Print acceptance rates for each target.
            for (unsigned int target = 0; target < data.ntargets; ++target) {
                double acceptance_rate = chain.cumulative_acceptance_probabilities(target) / static_cast<double>(iter+1);
                Rcpp::Rcout << "Metropolis acceptance rate for gamma (target " << target << ") = " 
                    << acceptance_rate << std::endl;
            }
        }

        // Store sample.
        if ((iter % thin) == 0) {
            sample_storage.store_sample(chain);
        }
    }

    // Measure sampling time in seconds.
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}

/**
 * Run the sampler for a specified number of iterations, and store the sampled beta and gamma values.
 * This is the main function that is called from R, and it unpacks the input data, initializes the 
 * MCMC chain and storage, runs the burn-in and sampling phases, and then packs the results into a 
 * list to return to R.
 * 
 * @param Xlist A list of design matrices for each target.
 * @param Ylist A list of response vectors for each target.
 * @param mean_prior The prior mean for the regression coefficients.
 * @param prec_prior The prior precision matrix for the regression coefficients.
 * @param ncategories A vector containing the number of categories for each target.
 * @param gammas_start A list of initial values for the gamma parameters for each target.
 * @param beta_start The initial values for the regression coefficients.
 * @param proposal_variance The proposal variance for the Metropolis-Hastings updates of the
 * gamma parameters.
 * @param iterations The number of MCMC iterations to perform during the sampling phase.
 * @param burnin The number of burn-in iterations to perform before the sampling phase.
 * @param thin The thinning interval for storing samples. Only every `thin`-th sample
 * will be stored in the storage matrices. To store all samples without thinning, set `thin` to 1.
 * @param seed The random seed to use for the GSL random number generator.
 * @param verbose The verbosity level for printing progress during burn-in and sampling. A value of
 * 0 means no progress will be printed, while higher values will print progress every `verbose` 
 * iterations.
 * 
 * @return An R list containing the results, including the stored beta and gamma samples, acceptance
 * rates, total iterations, and time taken for burn-in and sampling.
 */
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List cpp_hprobit(
    const Rcpp::List& Xlist,
    const Rcpp::List& Ylist,
    const arma::colvec& mean_prior,
    const arma::mat& prec_prior,
    const arma::ivec& ncategories,
    const Rcpp::List& gammas_start,
    const arma::colvec& beta_start,
    const arma::vec& proposal_variance,
    const int iterations,
    const int burnin,
    const int thin,
    const int seed,
    const int verbose
) {
    if (verbose > 0){
        Rcpp::Rcout << "Starting MH sampler for multi-scale probit model..." << std::endl;
    }
  
    //--- GSL random init ---
    // Used to take efficient samples from truncated normal.
    gsl_rng_env_setup();                          // Read variable environnement
    const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
    gsl_rng *gen = gsl_rng_alloc (type);          // Rand generator allocation
    gsl_rng_set(gen, seed);
  
    // Unpack data and define constants.
    const Data data = unpack_data(Xlist, Ylist);
    const int ntargets = data.ntargets;
    const int nstore = iterations / thin;

    // Sample storage.
    SampleStorage sampling_storage(nstore, data.npredictors, ncategories);

    // Create MCMC chain.
    MspmChain chain(
        1.0, // inv_temperature is 1.0 for the Gibbs sampler.
        beta_start,
        unpack_gamma(gammas_start, ncategories),
        proposal_variance,
        mean_prior,
        prec_prior,
        ncategories,
        data.nobs
    );

    if (verbose > 0 && burnin > 0) {
        Rcpp::Rcout << "Starting burn-in phase..." << std::endl;
    }

    // Do burnin.
    double burnin_time = 0;
    if (burnin > 0) {
        burnin_time = do_burnin(chain, data, burnin, gen, verbose);
    }

    // Reset proposal acceptance probabilities before sampling.
    chain.cumulative_acceptance_probabilities.zeros();

    // Do sampling.
    double sampling_time = do_sampling(iterations, chain, data, sampling_storage, thin, gen, verbose);

    // Compute acceptance rates for sampling phase.
    arma::vec acceptance_rate(ntargets, arma::fill::zeros);
    compute_acceptance_rate(chain.cumulative_acceptance_probabilities, iterations, acceptance_rate);

    // Free pointer to GSL random generator.
    gsl_rng_free(gen);
  
    return Rcpp::List::create(
        _["storebeta"] = sampling_storage.store_beta,
        _["storegamma"] = sampling_storage.gamma_to_r_list(),
        _["acceptance_rate"] = acceptance_rate,
        _["total_iter"] = iterations + burnin,
        _["sampling_time"] = sampling_time,
        _["burnin_time"] = burnin_time,
        _["nlikelihood_calls"] = chain.nlikelihood_calls
    );
}

/**
 * Tune the gibbs sampler by finding the optimal proposal variance for the gamma parameters. 
 * 
 * @param x_list A list of feature matrices for each target.
 * @param y_list A list of response vectors for each target.
 * @param mean_prior The prior mean for the regression coefficients.
 * @param prec_prior The prior precision matrix for the regression coefficients.
 * @param ncategories A vector containing the number of categories for each target.
 * @param gamma_start A list of initial values for the gamma parameters for each target.
 * @param beta_start The initial values for the regression coefficients.
 * @param tune_start The initial values for the tuning parameters for the proposal distribution
 * for the gammas.
 * @param target_acceptance_rate The target acceptance rate for the proposed gammas during tuning.
 * @param target_epsilon The epsilon threshold for early stopping during tuning. If the acceptance 
 * rates for all targets are within `target_epsilon` of the `target_acceptance_rate`, tuning will 
 * stop early.
 * @param stop_early A boolean indicating whether to stop tuning early if the acceptance rates for 
 * all targets are within `target_epsilon` of the `target_acceptance_rate`.
 * @param max_iterations The maximum number of tuning iterations to perform. The tuning process 
 * will stop after this many iterations even if the acceptance rates have not yet converged to the 
 * target.
 * @param window_size The window size for computing acceptance rates and adjusting the proposal 
 * variance during tuning.
 * @param seed The random seed to use for the GSL random number generator during tuning.
 * @param verbose The verbosity level for printing progress during tuning. A value of 0 means no
 * progress will be printed, while higher values will print progress every `verbose` iterations.
 * 
 * @return An R list containing the results of tuning, including the final proposal variance, 
 * acceptance rates, and the number of iterations performed during tuning.
 */
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List cpp_hprobit_tune(
    const Rcpp::List& x_list,
    const Rcpp::List& y_list,
    const arma::colvec& mean_prior,
    const arma::mat& prec_prior,
    const arma::ivec& ncategories,
    const Rcpp::List& gamma_start,
    const arma::colvec& beta_start,
    const arma::vec& tune_start,
    double target_acceptance_rate,
    double target_epsilon,
    bool stop_early,
    int max_iterations,
    int window_size,
    const int seed,
    int verbose
) {

    //--- GSL random init ---
    // Used to take efficient samples from truncated normal.
    gsl_rng_env_setup();                          // Read variable environnement
    const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
    gsl_rng *gen = gsl_rng_alloc (type);          // Rand generator allocation
    gsl_rng_set(gen, seed);
  
    // Unpack data and define constants.
    const Data data = unpack_data(x_list, y_list);
    const int ntargets = data.ntargets;

    // Create chain.
    MspmChain chain(
        1.0, // inv_temperature is 1.0 for the Gibbs sampler.
        beta_start,
        unpack_gamma(gamma_start, ncategories),
        tune_start,
        mean_prior,
        prec_prior,
        ncategories,
        data.nobs
    );

    // Create tuner.
    ProposalTuner tuner(target_acceptance_rate, target_epsilon, window_size, ntargets);

    // Tuning loop.
    int iter = 0;
    for (; iter < max_iterations; iter++) {
        chain.simulate_step(data, gen);

        tuner.tune_step(chain, data, gen);

        // Print progress.
        if (verbose > 0 && iter % verbose == 0) {
            Rcpp::Rcout << "Tuning iteration " << (iter+1) << " of " << max_iterations << " ";

            double mean_accept_rate = arma::mean(tuner.acceptance_rates);
            Rcpp::Rcout << "Mean acceptance rate = " << mean_accept_rate << std::endl;
        }

        // Check for early stopping.
        if (stop_early && tuner.has_reached_target()) {
            if (verbose > 0) {
                double mean_acceptance_rate = arma::mean(tuner.acceptance_rates);
                Rcpp::Rcout << "Early stopping at iteration " << (iter+1) 
                    << " as acceptance rates are within epsilon of target. " 
                    << "Mean acceptance rate = " << mean_acceptance_rate << std::endl;
            }
            break;
        }
    }

    // Free pointer to GSL random generator.
    gsl_rng_free(gen);

    return Rcpp::List::create(
        _["target_acceptance_rate"] = target_acceptance_rate,
        _["final_iteration"] = iter,
        _["proposal_variance"] = chain.proposal_variance,
        _["acceptance_rates"] = tuner.acceptance_rates
    );
}