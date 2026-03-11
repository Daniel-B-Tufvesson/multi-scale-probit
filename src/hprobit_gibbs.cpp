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
 * Do the burn-in phase.
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
 * sampler for a specified number of iterations, and stores the sampled beta and gamma values in
 * the provided storage matrices.
 * 
 * @param iterations The number of MCMC iterations to perform during the sampling phase.
 * @param beta The current regression coefficients, which will be updated in place with the new
 * sampled values during the sampling phase.
 * @param gamma The current threshold parameters for each target, which will be updated in place
 * with the new sampled values during the sampling phase.
 * @param storebeta The matrix to store the sampled beta values. This will be updated in place with 
 * the sampled beta values for each iteration during the sampling phase.
 * @param storegamma The vector of matrices to store the sampled gamma values for each target. This
 * will be updated in place with the sampled gamma values for each iteration during the sampling
 * phase.
 * @param thin The thinning interval for storing samples. Only every `thin`-th sample will be stored
 * in the storage matrices. To store all samples without thinning, set `thin` to 1.
 * @param acceptance_rate The vector to store the acceptance rates for the proposed gammas for each target.
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
    int nstored = 0;

    // Measure sampling time.
    auto start_time = std::chrono::high_resolution_clock::now();

    // Sampling loop.
    for (int iter = 0; iter < iterations; iter++) {
        chain.simulate_step(data, gen);
        
        // Print progress.
        if(verbose > 0 && iter % verbose == 0){
            double avg_acceptance_rate = 0.0;
            Rcpp::Rcout << "probit_gibbs iteration " << (iter+1) 
                << " of " << iterations << std::endl;

            // Print acceptance rates for each target.
            for (unsigned int target = 0; target < data.ntargets; ++target) {
                double acceptance_rate = chain.cumulative_acceptance_probabilities(target) / static_cast<double>(iter+1);
                avg_acceptance_rate += acceptance_rate;
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
        _["burnin_time"] = burnin_time
    );
}



/**
 * The main function for the MCMC sampler, which is called from R. This function unpacks the input data,
 * initializes the parameters and storage matrices, runs the burn-in and sampling phases of the MCMC
 * sampler, and then packs the results into a list to return to R.
 * 
 * @param Xlist A list of design matrices for each target.
 * @param Ylist A list of response vectors for each target.
 * @param meanPrior The prior mean for the regression coefficients.
 * @param precPrior The prior precision matrix for the regression coefficients.
 * @param ncat A vector containing the number of categories for each target.
 * @param gammaStart A list of initial values for the gamma parameters for each target.
 * @param betaStart The initial values for the regression coefficients.
 * @param tune_start The initial values for the tuning parameters for the proposal distribution 
 * for the gammas.
 * @param adapt_tune A boolean indicating whether to perform adaptive tuning of the proposal variance
 * during burn-in.
 * @param tune_window_size The window size for computing acceptance rates and adjusting the proposal
 * variance during adaptive burn-in.
 * @param target_acceptance_rate The target acceptance rate for the proposed gammas during adaptive
 * burn-in.
 * @param iterations The number of MCMC iterations to perform during the sampling phase.
 * @param burnin The number of burn-in iterations to perform before the sampling phase.
 * @param thin The thinning interval for storing samples. Only every `thin`-th sample will be stored
 * in the storage matrices. To store all samples without thinning, set `thin` to 1.
 * @param save_burnin_samples A boolean indicating whether to store the samples from the burn-in phase
 * in the storage matrices.
 * @param seed The random seed to use for the GSL random number generator.
 * @param verbose The verbosity level for printing progress during burn-in and sampling. A value of
 * 0 means no progress will be printed, while higher values will print progress every `verbose` 
 * iterations.
 * 
 * @return An R list containing the results.
 */

// Rcpp::List cpp_hprobit_old(
//     const Rcpp::List& Xlist,
//     const Rcpp::List& Ylist,
//     const arma::colvec& meanPrior,
//     const arma::mat& precPrior,
//     const arma::ivec& ncat,
//     const Rcpp::List& gammaStart,
//     const arma::colvec& betaStart,
//     const arma::vec& tune_start,
//     const bool adapt_tune,
//     int tune_window_size,
//     double target_acceptance_rate,
//     const int iterations,
//     const int burnin,
//     const int thin,
//     const bool save_burnin_samples,
//     const int seed,
//     const int verbose
// ) {
  
//     //--- GSL random init ---
//     // Used to take efficient samples from truncated normal.
//     gsl_rng_env_setup();                          // Read variable environnement
//     const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
//     gsl_rng *gen = gsl_rng_alloc (type);          // Rand generator allocation
//     gsl_rng_set(gen, seed);
  
//     // Unpack data and define constants.
//     const Data data = unpack_data(Xlist, Ylist);
//     const int ntargets = data.ntargets;
//     const int nstore = iterations/thin;
//     const int nstore_burnin = save_burnin_samples ? burnin : 0;

//     // Unpack gamma.
//     std::vector<colvec> gamma = unpack_gamma(gammaStart, ncat);
  
//     // Storage matrices
//     mat storebeta(nstore, betaStart.n_rows, arma::fill::zeros);
//     std::vector<mat> storegamma(ntargets);
//     for (unsigned int target = 0; target < ntargets; ++target) {
//         storegamma[target] = mat(nstore, ncat(target)-1, arma::fill::zeros);
//     }

//     // Storage matrices for burnin.
//     mat storebeta_burnin(nstore_burnin, betaStart.n_rows, arma::fill::zeros);
//     std::vector<mat> storegamma_burnin(ntargets);
//     if (save_burnin_samples) {
//         for (unsigned int target = 0; target < ntargets; ++target) {
//             storegamma_burnin[target] = mat(burnin, ncat(target)-1, arma::fill::zeros);
//         }
//     }
  
//     // Set starting points
//     arma::colvec beta = betaStart;
//     arma::vec tune = tune_start;
  
//     // Do burnin.
//     double burnin_duration = 0;
//     arma::vec burnin_acceptance_rate (data.ntargets, arma::fill::zeros);
//     if (adapt_tune) {
//         burnin_duration = do_adaptive_burnin(burnin, beta, gamma, data, ncat, tune, 
//             burnin_acceptance_rate, target_acceptance_rate, tune_window_size, meanPrior, 
//             precPrior, save_burnin_samples, storebeta_burnin, storegamma_burnin, 
//             thin, gen, verbose);
//     }
//     else {
//         burnin_duration = do_burnin(burnin, beta, gamma, data, ncat, tune, burnin_acceptance_rate, 
//             meanPrior, precPrior, save_burnin_samples, storebeta_burnin, storegamma_burnin,
//             thin, gen, verbose);
//     }

//     // Do sampling.
//     arma::vec acceptance_rate (data.ntargets, arma::fill::zeros);
//     double sampling_time = do_sampling(iterations, beta, gamma, data, ncat, tune,
//         meanPrior, precPrior, storebeta, storegamma, thin, acceptance_rate, 
//         gen, verbose);
  
//     // Pack stored gammas.
//     Rcpp::List gammas = Rcpp::List::create();
//     for (unsigned int target = 0; target < ntargets; ++target) {
//         gammas.push_back(storegamma[target]);
//     }

//     // Pack stored burnin gammas.
//     Rcpp::List storegamma_burnin_list = Rcpp::List::create();
//     if (save_burnin_samples) {
//         for (unsigned int target = 0; target < ntargets; ++target) {
//             storegamma_burnin_list.push_back(storegamma_burnin[target]);
//         }
//     }

//     // Free pointer to GSL random generator.
//     gsl_rng_free(gen);
  
//     return Rcpp::List::create(
//         _["storebeta"] = storebeta,
//         _["storegamma"] = gammas,
//         _["tune"] = tune,
//         _["acceptance_rate"] = acceptance_rate,
//         _["burnin_acceptance_rate"] = burnin_acceptance_rate,
//         _["total_iter"] = iterations + burnin,
//         _["storebeta_burnin"] = storebeta_burnin,
//         _["storegamma_burnin"] = storegamma_burnin_list,
//         _["sampling_time"] = sampling_time,
//         _["burnin_time"] = burnin_duration
//     );
// }

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
 * @return An R list containing the tuned proposal variances and acceptance rates for each target, 
 * as well as information about the tuning process such as the number of iterations performed and 
 * whether early stopping was triggered.
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

    arma::vec acceptance_rates (data.ntargets, arma::fill::zeros);
    int window_step = 0;

    // Tuning loop.
    int iter = 0;
    int last_print = 0;
    for (; iter < max_iterations; iter++) {
        chain.simulate_step(data, gen);

        // Tune the proposal variance.
        window_step++;
        if (window_step == window_size) {
            window_step = 0;

            // Do adjustment.
            double learning_rate = 1.0 / std::sqrt(iter);
            compute_acceptance_rate(chain.cumulative_acceptance_probabilities, window_size, 
                acceptance_rates);
            chain.adjust_proposal_variance(acceptance_rates, target_acceptance_rate, learning_rate);
            chain.cumulative_acceptance_probabilities.zeros();

            // Print progress.
            if (verbose > 0 && (iter - last_print) >= verbose) {
                last_print = iter;
                Rcpp::Rcout << "Tuning iteration " << (iter+1) << " of " << max_iterations << " ";

                double mean_accept_rate = arma::mean(acceptance_rates);
                Rcpp::Rcout << "Mean acceptance rate = " << mean_accept_rate << std::endl;
            }

            // Check for early stopping.
            if (stop_early) {
                bool all_close = true;
                for (int target = 0; target < ntargets; target++) {
                    if (std::abs(acceptance_rates(target) - target_acceptance_rate) > target_epsilon) {
                        all_close = false;
                        break;
                    }
                }
                if (all_close) {
                    if (verbose > 0) {
                        double mean_accept_rate = arma::mean(acceptance_rates);
                        Rcpp::Rcout << "Early stopping at iteration " << (iter+1) 
                            << " as acceptance rates are within epsilon of target. " 
                            << "Mean acceptance rate = " << mean_accept_rate << std::endl;
                    }
                    break;
                }
            }
        }
    }

    // Free pointer to GSL random generator.
    gsl_rng_free(gen);

    return Rcpp::List::create(
        _["target_acceptance_rate"] = target_acceptance_rate,
        _["target_epsilon"] = target_epsilon,
        _["max_iterations"] = max_iterations,
        _["final_iteration"] = iter,
        _["final_tune"] = chain.proposal_variance,
        _["final_acceptance_rates"] = acceptance_rates
    );
}