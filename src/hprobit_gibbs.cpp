/**
 * @file hprobit_gibbs.cpp
 * @author Johan Falkenjack
 * @author Daniel Tufvesson
 * @version 3.0
 * @brief A Metropolis-Hastings-withing-Gibbs sampler for the Multi-Scale Probit model.
 */

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "sampling.hpp"
#include <chrono>

/**
 * Do a single MCMC step and update the beta and gammas.
 * 
 * @param iter The current iteration number of the MCMC sampler.
 * @param ncategories The number of categories for each target.
 * @param beta The current regression coefficients, which will be updated in place with the new 
 * sampled values.
 * @param gamma The current threshold parameters for each target, which will be updated in place
 * with the new sampled values.
 * @param tune The tuning parameters for the proposal distribution for the gammas. This controls
 * the standard deviation of the truncated normal distribution used for proposing new threshold
 * values.
 * @param data The data object.
 * @param acceptance_probabilities The vector to store the acceptance probabilities for the 
 * proposed gammas for each target. This will be updated in place with the acceptance probabilities
 * for the proposed gammas for each target in this iteration.
 * @param meanPrior The prior mean for the regression coefficients.
 * @param precPrior The prior precision matrix for the regression coefficients.
 * @param gen The GSL random number generator to use for sampling.
 */
void do_step(
    int iter,
    const arma::ivec& ncategories,
    arma::colvec& beta,
    std::vector<colvec>& gamma,
    const arma::vec& tune,
    const Data& data,
    arma::vec& acceptance_probabilities,
    const arma::colvec& meanPrior,
    const arma::mat& precPrior,
    gsl_rng* gen
) {
    // Step 1: Update gammas with Metropolis-Hastings.
    for (unsigned int t = 0; t < data.ntargets; ++t) {
        int target = (iter+t) % data.ntargets;

        // Make proposal for gamma.
        double acceptance_prob = 0;
        bool accepted = mh_update_gamma(
            gamma[target],
            beta,
            data.X[target],
            data.Y[target],
            ncategories(target),
            tune(target),
            1.0, // inv_temperature is 1.0 for the Gibbs sampler.
            acceptance_prob,
            gen
        );
        acceptance_probabilities(target) += acceptance_prob;
    }

    // Step 2: Update beta with Gibbs.
    gibbs_update_beta(beta, gamma, data, meanPrior, precPrior, gen);
}

void store_sample(
    const arma::colvec& beta,
    const std::vector<colvec>& gamma,
    arma::mat& storebeta,
    std::vector<mat>& storegamma,
    int& nstored
) {
    // Store beta.
    for (unsigned int j = 0; j < beta.n_rows; j++) {
        storebeta(nstored, j) = beta[j];
    }
    // Store gamma.
    for (unsigned int target = 0; target < gamma.size(); target++) {
        // Note: the first and last gammas should not be saved.
        for (unsigned int j = 1; j < gamma[target].n_rows - 1; j++) {
            storegamma[target](nstored, j-1) = gamma[target](j);
        }
    }
    nstored++;
}

/**
 * Do a single burnin step, which involves doing a MCMC step and then storing the sample if needed.
 */
void do_burnin_step(
    int iter,
    int burnin,
    const arma::ivec& ncategories,
    arma::colvec& beta,
    std::vector<colvec>& gamma,
    const arma::vec& tune,
    const Data& data,
    arma::vec& acceptance_probabilites,
    const arma::colvec& meanPrior,
    const arma::mat& precPrior,
    bool save_burnin_samples,
    arma::mat& storebeta_burnin,
    std::vector<arma::mat>& storegamma_burnin,
    int& nstored,
    int thin,
    gsl_rng* gen,
    int verbose
) {
    do_step(iter, ncategories, beta, gamma, tune, data, acceptance_probabilites, meanPrior,
        precPrior, gen);

    // Print progress.
    if(verbose > 0 && (iter % verbose) == 0){
        Rcpp::Rcout << "probit_gibbs burnin iteration " << (iter+1) 
            << " of " << burnin << std::endl;

        // Print acceptance rates for each target.
        for (unsigned int target = 0; target < data.ntargets; ++target) {
            double acceptance_rate = acceptance_probabilites(target) / static_cast<double>(iter+1);
            Rcpp::Rcout << "Metropolis acceptance rate for gamma (target " << target << ") = " 
                << acceptance_rate << std::endl;
        }
    }

    // Store burnin sample.
    if (save_burnin_samples && (iter % thin) == 0) {
        store_sample(beta, gamma, storebeta_burnin, storegamma_burnin, nstored);
    }
}

void compute_acceptance_rate(
    const arma::vec& acceptance_probabilities,
    int nsamples,
    arma::vec& acceptance_rates
) {
    for (unsigned int target = 0; target < acceptance_probabilities.n_rows; ++target) {
        acceptance_rates(target) = acceptance_probabilities(target) / static_cast<double>(nsamples);
    }
}

/**
 * Adjust the proposal variance for the gamma parameters based on the observed acceptance 
 * probabilities for the proposed gamma values. This is done by comparing the observed acceptance 
 * rates to a target acceptance rate, and adjusting the standard deviation of the truncated normal 
 * proposal distribution for the gammas accordingly.
 * 
 * @param tune The current tuning parameters for the proposal distribution for the gammas, which 
 * will be updated in place with the new tuned values. Each sigma element corresponds to a target.
 * @param acceptance_probabilites The sum of the acceptance probabilities for the proposed gammas 
 * for each target over the current window of burnin iterations. Each sum corresponed to a target.
 * @param target_acceptance_rate The target acceptance rate for the proposed gammas.
 * @param learning_rate The learning rate for adjusting the proposal variance.
 */
void adjust_proposal_variance(
    arma::vec& tune,
    const arma::vec& burnin_acceptance_rate,
    double target_acceptance_rate,
    double learning_rate
) {
    for (int i = 0; i < tune.n_rows; i++) {
        double acceptance_rate = burnin_acceptance_rate(i);
        double log_sigma = std::log(tune(i));
        log_sigma += learning_rate * (acceptance_rate - target_acceptance_rate);
        tune(i) = std::exp(log_sigma);
    }
}

/**
 * Do adaptive burnin by running the MCMC sampler for a specified number of burnin iterations, and 
 * automatically tuning the proposal distribution for the gammas based on the acceptance rates for 
 * the proposed gammas. The tuning is done by adjusting the standard deviation of the truncated 
 * Gaussian proposal distribution for the gammas based on the observed acceptance rates in windows 
 * of burnin iterations.
 * 
 * @param burnin The number of burn-in iterations to perform.
 * @param beta The current regression coefficients, which will be updated in place with the new
 * sampled values during burn-in.
 * @param gamma The current threshold parameters for each target, which will be updated in place
 * with the new sampled values during burn-in.
 */
double do_adaptive_burnin(
    int burnin,
    arma::colvec& beta,
    std::vector<arma::colvec>& gamma,
    const Data& data,
    const arma::ivec& ncategories,
    arma::vec& tune,
    arma::vec& burnin_acceptance_rate,
    double target_acceptance_rate,
    int window_size,
    const arma::colvec& meanPrior,
    const arma::mat& precPrior,
    bool save_burnin_samples,
    mat& storebeta_burnin,
    std::vector<mat>& storegamma_burnin,
    int thin,
    gsl_rng* gen,
    int verbose
) {
    arma::vec acceptance_probabilites (data.ntargets, arma::fill::zeros);
    int nstored = 0;
    int window_step = 0;

    // Measure burnin sampling time.
    auto start_time_burnin = std::chrono::high_resolution_clock::now();

    // Burnin loop.
    for (int iter = 0; iter < burnin; iter++) {
        do_burnin_step(iter, burnin, ncategories, beta, gamma, tune, data, acceptance_probabilites,
            meanPrior, precPrior, save_burnin_samples, storebeta_burnin, storegamma_burnin, nstored,
            thin, gen, verbose);

        // Tune the proposal variance.
        window_step++;
        if (window_step == window_size) {
            window_step = 0;
            double learning_rate = 1.0 / std::sqrt(iter);
            compute_acceptance_rate(acceptance_probabilites, window_size, burnin_acceptance_rate);
            adjust_proposal_variance(
                tune,
                burnin_acceptance_rate,
                target_acceptance_rate,
                learning_rate
            );
            acceptance_probabilites.zeros();
        }

    }
    // Measure burnin time in seconds.
    auto end_time_burnin = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_burnin = end_time_burnin - start_time_burnin;
    return elapsed_burnin.count();
}

/**
 * Do regular burn-in without any adaptation.
 */
double do_burnin(
    int burnin,
    arma::colvec& beta,
    std::vector<arma::colvec>& gamma,
    const Data& data,
    const arma::ivec& ncategories,
    arma::vec& tune,
    arma::vec& burnin_acceptance_rate,
    const arma::colvec& meanPrior,
    const arma::mat& precPrior,
    bool save_burnin_samples,
    arma::mat& storebeta_burnin,
    std::vector<arma::mat>& storegamma_burnin,
    int thin,
    gsl_rng* gen,
    int verbose
) {
    arma::vec acceptance_probabilites (data.ntargets, arma::fill::zeros);
    int nstored = 0;

    // Measure burnin sampling time.
    auto start_time_burnin = std::chrono::high_resolution_clock::now();

    // Burnin loop.
    for (int iter = 0; iter < burnin; iter++) {
        do_burnin_step(iter, burnin, ncategories, beta, gamma, tune, data, acceptance_probabilites,
            meanPrior, precPrior, save_burnin_samples, storebeta_burnin, storegamma_burnin,
            nstored, thin, gen, verbose);
    }
    // Measure burnin time in seconds.
    auto end_time_burnin = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_burnin = end_time_burnin - start_time_burnin;

    compute_acceptance_rate(acceptance_probabilites, burnin, burnin_acceptance_rate);

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
    arma::colvec& beta,
    std::vector<arma::colvec>& gamma,
    const Data& data,
    const arma::ivec& ncategories,
    arma::vec& tune,
    const arma::colvec& meanPrior,
    const arma::mat& precPrior,
    mat& storebeta,
    std::vector<mat>& storegamma,
    int thin,
    arma::vec& acceptance_rate,
    gsl_rng* gen,
    int verbose
) {
    int nstored = 0;
    arma::vec acceptance_probabilites (data.ntargets, arma::fill::zeros);

    // Measure sampling time.
    auto start_time = std::chrono::high_resolution_clock::now();

    // Sampling loop.
    for (int iter = 0; iter < iterations; iter++) {

        do_step(iter, ncategories, beta, gamma, tune, data, acceptance_probabilites,
            meanPrior, precPrior, gen);
        
        // Print progress.
        if(verbose > 0 && iter % verbose == 0){
            double avg_acceptance_rate = 0.0;
            Rcpp::Rcout << "probit_gibbs iteration " << (iter+1) 
                << " of " << iterations << std::endl;

            // Print acceptance rates for each target.
            for (unsigned int target = 0; target < data.ntargets; ++target) {
                double acceptance_rate = acceptance_probabilites(target) / static_cast<double>(iter+1);
                avg_acceptance_rate += acceptance_rate;
                Rcpp::Rcout << "Metropolis acceptance rate for gamma (target " << target << ") = " 
                    << acceptance_rate << std::endl;
            }
        }

        if ((iter % thin) == 0) {
            store_sample(beta, gamma, storebeta, storegamma, nstored);
        }
    }

    // Measure sampling time in seconds.
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    compute_acceptance_rate(acceptance_probabilites, iterations, acceptance_rate);

    return elapsed.count();
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
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List cpp_hprobit(
    const Rcpp::List& Xlist,
    const Rcpp::List& Ylist,
    const arma::colvec& meanPrior,
    const arma::mat& precPrior,
    const arma::ivec& ncat,
    const Rcpp::List& gammaStart,
    const arma::colvec& betaStart,
    const arma::vec& tune_start,
    const bool adapt_tune,
    int tune_window_size,
    double target_acceptance_rate,
    const int iterations,
    const int burnin,
    const int thin,
    const bool save_burnin_samples,
    const int seed,
    const int verbose
) {
  
    //--- GSL random init ---
    // Used to take efficient samples from truncated normal.
    gsl_rng_env_setup();                          // Read variable environnement
    const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
    gsl_rng *gen = gsl_rng_alloc (type);          // Rand generator allocation
    gsl_rng_set(gen, seed);
  
    // Unpack data and define constants.
    const Data data = unpack_data(Xlist, Ylist);
    const int ntargets = data.ntargets;
    const int nstore = iterations/thin;
    const int nstore_burnin = save_burnin_samples ? burnin : 0;

    // Unpack gamma.
    std::vector<colvec> gamma(ntargets);
    for (unsigned int target = 0; target < ntargets; ++target) {
        gamma[target] = Rcpp::as<colvec>(gammaStart[target]);
        gamma[target](0) = -std::numeric_limits<double>::max();
        gamma[target](ncat[target]) = std::numeric_limits<double>::max();
    }
  
    // Storage matrices
    mat storebeta(nstore, betaStart.n_rows, arma::fill::zeros);
    std::vector<mat> storegamma(ntargets);
    for (unsigned int target = 0; target < ntargets; ++target) {
        storegamma[target] = mat(nstore, ncat(target)-1, arma::fill::zeros);
    }

    // Storage matrices for burnin.
    mat storebeta_burnin(nstore_burnin, betaStart.n_rows, arma::fill::zeros);
    std::vector<mat> storegamma_burnin(ntargets);
    if (save_burnin_samples) {
        for (unsigned int target = 0; target < ntargets; ++target) {
            storegamma_burnin[target] = mat(burnin, ncat(target)-1, arma::fill::zeros);
        }
    }
  
    // Set starting points
    arma::colvec beta = betaStart;
    arma::vec tune = tune_start;
  
    // Do burnin.
    double burnin_duration = 0;
    arma::vec burnin_acceptance_rate (data.ntargets, arma::fill::zeros);
    if (adapt_tune) {
        burnin_duration = do_adaptive_burnin(burnin, beta, gamma, data, ncat, tune, 
            burnin_acceptance_rate, target_acceptance_rate, tune_window_size, meanPrior, 
            precPrior, save_burnin_samples, storebeta_burnin, storegamma_burnin, 
            thin, gen, verbose);
    }
    else {
        burnin_duration = do_burnin(burnin, beta, gamma, data, ncat, tune, burnin_acceptance_rate, 
            meanPrior, precPrior, save_burnin_samples, storebeta_burnin, storegamma_burnin,
            thin, gen, verbose);
    }

    // Do sampling.
    arma::vec acceptance_rate (data.ntargets, arma::fill::zeros);
    double sampling_time = do_sampling(iterations, beta, gamma, data, ncat, tune,
        meanPrior, precPrior, storebeta, storegamma, thin, acceptance_rate, 
        gen, verbose);
  
    // Pack stored gammas.
    Rcpp::List gammas = Rcpp::List::create();
    for (unsigned int target = 0; target < ntargets; ++target) {
        gammas.push_back(storegamma[target]);
    }

    // Pack stored burnin gammas.
    Rcpp::List storegamma_burnin_list = Rcpp::List::create();
    if (save_burnin_samples) {
        for (unsigned int target = 0; target < ntargets; ++target) {
            storegamma_burnin_list.push_back(storegamma_burnin[target]);
        }
    }

    // Free pointer to GSL random generator.
    gsl_rng_free(gen);
  
    return Rcpp::List::create(
        _["storebeta"] = storebeta,
        _["storegamma"] = gammas,
        _["tune"] = tune,
        _["acceptance_rate"] = acceptance_rate,
        _["burnin_acceptance_rate"] = burnin_acceptance_rate,
        _["total_iter"] = iterations + burnin,
        _["storebeta_burnin"] = storebeta_burnin,
        _["storegamma_burnin"] = storegamma_burnin_list,
        _["sampling_time"] = sampling_time,
        _["burnin_time"] = burnin_duration
    );
}
