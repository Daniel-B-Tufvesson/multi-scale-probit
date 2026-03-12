/**
 * @file hprobit_pt.cpp
 * @author Daniel Tufvesson
 * @version 3.0
 * @brief A Metropolis-Hastings-withing-Gibbs sampler for the Multi-Scale Probit model, 
 * that uses parallel tempering to improve mixing.
 */

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "sampling.hpp"

/**
 * A class managing several chains at different temperatures and performing parallel tempering 
 * state swaps between them.
 */
class PTChains {
public:
    /** The parallel chains. First chain (canonical) is the lowest temperature chain. */
    std::vector<MspmChain> chains;

    /** The number of temperatures. */
    const int ntemperatures;

    bool complete_swapping;

    /** The cumulative swap probabilities for each adjacent chain pair. */
    arma::vec cumulative_swap_probabilities;

    /** Counter of number of swap proposals for each chain pair. */
    arma::ivec nswap_proposals;

    PTChains(
        const arma::vec& inv_temperature_ladder,
        const arma::colvec& beta_start,
        const std::vector<arma::colvec>& gamma_start,
        const std::vector<arma::vec> proposal_variances,
        const arma::colvec& beta_mean_prior,
        const arma::mat& beta_prec_prior,
        const arma::ivec& ncategories,
        const int total_nobs,
        bool complete_swapping
    ) : ntemperatures(inv_temperature_ladder.n_elem), complete_swapping(complete_swapping) {

        // Create the chains for each temperature in the ladder.
        for (int i = 0; i < ntemperatures; i++) {
            chains.emplace_back(
                inv_temperature_ladder(i),
                beta_start,
                gamma_start,
                proposal_variances[i],
                beta_mean_prior,
                beta_prec_prior,
                ncategories,
                total_nobs
            );
        }

        cumulative_swap_probabilities = arma::vec(ntemperatures-1, arma::fill::zeros);
        nswap_proposals = arma::ivec(ntemperatures-1, arma::fill::zeros);
    }

    void simulate_step(const Data& data, gsl_rng* rng) {
        // Do within-temperature MH updates.
        for (auto& chain : chains) {
            chain.simulate_step(data, rng);
        }

        // Select pair of adjacent temperatures to swap.
        int swap_index = gsl_rng_uniform_int(rng, ntemperatures - 1);
        MspmChain& chain1 = chains[swap_index];
        MspmChain& chain2 = chains[swap_index + 1];

        double swap_probability = chain1.compute_swap_probability(chain2, data, rng);
        if (gsl_ran_flat(rng, 0.0, 1.0) <= swap_probability) {
            if (complete_swapping) {
                chain1.swap_beta(chain2);
            }
            chain1.swap_gammas(chain2);
            
        }
        cumulative_swap_probabilities(swap_index) += swap_probability;
        nswap_proposals(swap_index) += 1;
    }

    arma::vec get_inv_temperatures() {
        arma::vec temperatures(ntemperatures);
        for (int i = 0; i < ntemperatures; i++) {
            temperatures(i) = chains[i].inv_temperature;
        }
        return temperatures;
    }

    Rcpp::List get_proposal_variance() {
        Rcpp::List proposal_variances(ntemperatures);
        for (int i = 0; i < ntemperatures; i++) {
            proposal_variances[i] = chains[i].proposal_variance;
        }
        return proposal_variances;
    }
};

/**
 * A class for tuning the proposal variance of the MSPM parallel tempering sampler. This operates
 * on several chains.
 */
class PTProposalTuner {
public:

    std::vector<ProposalTuner> tuners;

    PTProposalTuner(
        double target_acceptance_rate,
        double target_epsilon,
        int window_size,
        int ntargets,
        int ntemperatures
    ) {
        for (int i = 0; i < ntemperatures; i++) {
            tuners.emplace_back(target_acceptance_rate, target_epsilon, window_size, ntargets);
        }
    }

    /**
     * Run one step of the tuning process, which includes updating the acceptance probabilities for 
     * the proposed gammas, and adjusting the proposal variance if the end of the window is reached.
     * 
     * Note that this function does not simulate a step of the chain itself. Simulating a step of
     * the chain should be done before calling this function.
     * 
     * @param chains The PTChains instance for which to perform the tuning step. 
     * @param data The data object containing the feature matrices and response vectors for each target.
     * @param rng The GSL random number generator to use for sampling.
     */
    void tune_step(PTChains& chains, const Data& data, gsl_rng* rng) {
        for (int i = 0; i < chains.ntemperatures; i++) {
            tuners[i].tune_step(chains.chains[i], data, rng);
        }
    }

    /**
     * Check whether the acceptance rates for the proposed gammas have reached the target acceptance 
     * rate within the specified epsilon threshold for all targets. This can be used for early stopping 
     * of the tuning process.
     * 
     * @return A boolean indicating whether the acceptance rates for all targets have reached the 
     * target acceptance rate within the epsilon threshold (true) or not (false).
     */
    bool has_reached_target() {
        for (auto& tuner : tuners) {
            if (!tuner.has_reached_target()) {
                return false;
            }
        }
        return true;
    }

    /**
     * Extract the acceptance rates for each chain.
     * 
     * @return A list of acceptance rates for each chain. Each element in the list is a vector 
     * of acceptance rates for the proposed gammas for each target in that chain.
     */
    Rcpp::List get_acceptance_rates() {
        Rcpp::List acceptance_rates(tuners.size());
        for (int i = 0; i < tuners.size(); i++) {
            acceptance_rates[i] = tuners[i].acceptance_rates;
        }
        return acceptance_rates;
    }
};

/**
 * A class for tuning the temperature ladder by equalizing the swap rates between the chains to 
 * a target rate.
 */
class TemperatureLadderTuner {
public:
    /**The target swap acceptance ratio that we want to achieve for each pair of adjacent 
     * temperatures. */
    const double target_swap_rate;

    /** The epsilon threshold for early stopping during tuning. */
    const double target_epsilon;

    /** The size of the window. This grow over time. */
    int window_size;

    const double window_growth_factor;

    /** The learning rate for adjusting the ladder gaps. */
    const double learning_rate;

    /** The minimum allowed gap between inverse temperatures in the ladder. */
    const double min_gap;

    /** The number of chains. */
    const int ntemperatures;

    /** The total number of iterations the tuner has run. */
    int step = 0;

    /** The number of iterations the tuner has taken within the current window. */
    int window_step = 0;

    /** The swap rate for each chain pair. */
    arma::colvec swap_rates;

    TemperatureLadderTuner(
        double target_swap_rate,
        double target_epsilon,
        int window_size,
        double window_growth_factor,
        double learning_rate,
        double min_gap,
        int ntemperatures
    ) : target_swap_rate(target_swap_rate), target_epsilon(target_epsilon), 
        window_size(window_size), window_growth_factor(window_growth_factor), 
        learning_rate(learning_rate), min_gap(min_gap), ntemperatures(ntemperatures) {
        
        swap_rates = arma::colvec(ntemperatures - 1, arma::fill::zeros);
    }

    void tune_step(PTChains& chains, const Data& data, gsl_rng* rng) {
        window_step++;
        if (window_step == window_size) {
            adjust_ladder(chains);

            // Increase window size.
            window_step = 0;
            window_size *= window_growth_factor;
            
            // Reset counters.
            chains.cumulative_swap_probabilities.zeros();
            chains.nswap_proposals.zeros();
        }
        step++;
    }

    /**
     * Check whether the swap rates between all pairs of adjacent temperatures have reached the 
     * target swap acceptance rate within the specified epsilon threshold.
     * 
     * @return A boolean indicating whether the swap acceptance rates for all pairs of adjacent 
     * temperatures have reached the target swap acceptance rate within the epsilon threshold 
     * (true) or not (false).
     */
    bool has_reached_target() {
        for (int i = 0; i < ntemperatures - 1; i++) {
            if (std::abs(swap_rates(i) - target_swap_rate) > target_epsilon) {
                return false;
            }
        }
        return true;
    }

private:
    /**
     * Adjust the temperature ladder based on the accumulated swap probabilities for each pair of
     * adjacent temperatures.
     * 
     * @param chains The PTChains instance representing the Markov chains at each temperature. The 
     * inverse temperatures (betas) in these chains will be updated based on the new ladder.
     */
    void adjust_ladder(PTChains& chains) {
        std::vector<double> ladder_gaps(ntemperatures - 1);

        // Compute old gaps.
        for(int i = 0; i < ntemperatures - 1; i++) {
            ladder_gaps[i] = chains.chains[i].inv_temperature - chains.chains[i+1].inv_temperature;
        }

        // Compute mean swap rate for each pair.
        for (int i = 0; i < ntemperatures - 1; i++) {
            swap_rates[i] = chains.cumulative_swap_probabilities[i] / 
                (chains.nswap_proposals[i] == 0 ? 1 : chains.nswap_proposals[i]);
        }

        // Update ladder gaps.
        for (int i = 0; i < ntemperatures - 1; i++) {
            ladder_gaps[i] *= std::exp(learning_rate * (swap_rates[i] - target_swap_rate));

            // Impose min.
            if (ladder_gaps[i] < min_gap) {
                ladder_gaps[i] = min_gap;
            }
        }

        // Normalize the gaps to ensure beta_K is fixed.
        double total_span = ladder_gaps[0];
        for (int i = 1; i < ntemperatures - 1; i++) {
            total_span += ladder_gaps[i];
        }
        double normalization_factor = (1.0 - chains.chains[ntemperatures-1].inv_temperature) / total_span;
        for (int i = 0; i < ntemperatures - 1; i++) {
            ladder_gaps[i] *= normalization_factor;
        }

        // Reconstruct ladder.
        chains.chains[0].inv_temperature = 1.0;
        for (int i = 1; i < ntemperatures-1; i++) {
            chains.chains[i].inv_temperature = chains.chains[i-1].inv_temperature - ladder_gaps[i-1];
        }
    }
};


/**
 * Perform the burn-in phase of the sampler. This is similar to the do_sampling function, but does 
 * not store any samples.
 * 
 * @param chains The PTChains instance representing the Markov chains at each temperature.
 * @param data The data object containing the feature matrices and response vectors for each target.
 * @param burnin The number of burn-in iterations to perform.
 * @param rng The GSL random number generator to use for sampling.
 * 
 * @return The total time taken for the burn-in phase in seconds.
 */
double do_burnin(
    PTChains& chains, 
    const Data& data,
    int burnin, 
    gsl_rng* rng,
    int verbose
) {
    // Measure sampling time.
    auto start_time = std::chrono::high_resolution_clock::now();

    // Burnin loop.
    for (unsigned int iter = 0; iter < burnin; iter++) {
        chains.simulate_step(data, rng);

        if (verbose > 0 && iter > 0 && (iter % verbose == 0 || iter == burnin - 1)) {
            Rcpp::Rcout << "Burn-in iteration " << (iter+1) << "/" << burnin << std::endl;
        }
    }

    // Measure sampling time in seconds.
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}

/**
 * Do the parallel tempering simulation and store samples to sample storage.
 * 
 * @param chains The PTChains instance representing the Markov chains at each temperature.
 * @param data The data object containing the feature matrices and response vectors for each target.
 * @param sample_storage The SampleStorage instance where the samples from the lowest temperature
 * chain will be stored.
 * @param iterations The total number of sampling iterations to perform.
 * @param thin The thinning interval for storing samples.
 * @param rng The GSL random number generator to use for sampling.
 * @param verbose The interval for printing progress updates to the console. If 0, no updates will
 * be printed.
 * 
 * @return The total time taken for the sampling phase in seconds.
 */
double do_sampling(
    PTChains& chains, 
    const Data& data,
    SampleStorage& sample_storage,
    int iterations,
    int thin,
    gsl_rng* rng,
    int verbose
) {
    // Measure sampling time.
    auto start_time = std::chrono::high_resolution_clock::now();

    // Sampling loop.
    for (unsigned int iter = 0; iter < iterations; iter++) {
        chains.simulate_step(data, rng);

        if (iter % thin == 0) {
            sample_storage.store_sample(chains.chains[0]);
        }

        if (verbose > 0 && iter > 0 && (iter % verbose == 0 || iter == iterations - 1)) {
            Rcpp::Rcout << "Sampling iteration " << (iter+1) << "/" << iterations << std::endl;
        }
    }

    // Measure sampling time in seconds.
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}


/**
 * Unpack the proposal variances from the Rcpp::List format to a std::vector of arma::vec.
 * 
 * @param proposal_variances The proposal variances in Rcpp::List format, where each element is a 
 * numeric vector of proposal variances for the gamma parameters for that chain.
 * @return A std::vector of arma::vec, where each element is a vector of proposal variances for 
 * the gamma parameters for that chain.
 */
std::vector<arma::vec> unpack_proposal_variances(const Rcpp::List& proposal_variances) {
    std::vector<arma::vec> proposal_variances_vec;
    for (int i = 0; i < proposal_variances.length(); i++) {
        proposal_variances_vec.push_back(Rcpp::as<arma::vec>(proposal_variances[i]));
    }
    return proposal_variances_vec;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List cpp_hprobit_pt(
    const Rcpp::List& xlist,
    const Rcpp::List& ylist,
    const arma::colvec& beta_mean_prior,
    const arma::mat& beta_prec_prior,
    const arma::ivec& ncategories,
    const Rcpp::List& gamma_start,
    const arma::colvec& beta_start,
    const Rcpp::List& proposal_variances,
    const arma::vec inv_temperature_ladder,
    const int iterations,
    const int burnin,
    const int thin,
    const int seed,
    const bool complete_swapping,
    const int verbose
) {
    if (verbose > 0){
        Rcpp::Rcout << "Starting parallel tempering sampler for multi-scale probit model..." << std::endl;
    }
    // Initialize GSL random number generator.
    gsl_rng_env_setup();                          // Read variable environnement
    const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
    gsl_rng* gen = gsl_rng_alloc(type);           // Rand generator allocation
    gsl_rng_set(gen, seed);

    // Define constants and unpack data.
    const Data data = unpack_data(xlist, ylist);
    const int nstore = iterations / thin;

    // Sample storages.
    SampleStorage sampling_storage(nstore, data.npredictors, ncategories);

    // Create chains.
    PTChains chains(
        inv_temperature_ladder,
        beta_start,
        unpack_gamma(gamma_start, ncategories),
        unpack_proposal_variances(proposal_variances),
        beta_mean_prior,
        beta_prec_prior,
        ncategories,
        data.nobs,
        complete_swapping
    );

    // Do burnin.
    double burnin_time = 0;
    if (burnin > 0) {
        if (verbose > 0) {
            Rcpp::Rcout << "Starting burn-in phase..." << std::endl;
        }

        burnin_time = do_burnin(chains, data, burnin, gen, verbose);

        if (verbose > 0) {
            Rcpp::Rcout << "Burn-in complete. Starting sampling phase..." << std::endl;
        }
    }

    // Do sampling.
    double sampling_time = do_sampling(chains, data, sampling_storage, iterations, thin, 
        gen, verbose);

    // Free pointers.
    gsl_rng_free(gen);

    // Return results.
    return Rcpp::List::create(
        _["storebeta"] = sampling_storage.store_beta,
        _["storegamma"] = sampling_storage.gamma_to_r_list(),
        _["sampling_time"] = sampling_time,
        _["burnin_time"] = burnin_time
    );
}


/**
 * Tune the proposal variance and temperature ladder for the parallel tempering sampler for the
 * multi-scale probit model. This function runs the sampler for a specified number of iterations, 
 * and adjusts the proposal variance and temperature ladder at specified intervals to try to reach 
 * target acceptance rates for the proposed gammas and the temperature swaps. This can be used to 
 * find good tuning parameters for the sampler before running the actual sampling function.
 * 
 * @param x_list The list of feature matrices for each target.
 * @param y_list The list of response vectors for each target.
 * @param mean_prior The mean vector for the normal prior on the beta coefficients.
 * @param prec_prior The precision matrix for the normal prior on the beta coefficients.
 * @param ncategories The number of categories for each target.
 * @param gamma_start The initial values for the gamma parameters for each target and chain.
 * @param beta_start The initial values for the beta coefficients for each chain.
 * @param proposal_variance The initial proposal variances for the gamma parameters for each target 
 * and chain. This should be a list of numeric vectors, where each element in the list corresponds 
 * to a chain, and contains a numeric vector of proposal variances for the gamma parameters for 
 * each target in that chain.
 * @param tune_proposal_variance A boolean indicating whether to tune the proposal variance during
 * the iterations.
 * @param target_acceptance_rate The target acceptance rate for the proposed gammas that the tuner
 * will try to achieve by adjusting the proposal variance.
 * @param target_acceptance_epsilon The epsilon threshold for checking whether the acceptance rates
 * have reached the target acceptance rate for the proposed gammas.
 * @param proposal_window_size The size of the window for tuning the proposal variance. The tuner 
 * will adjust the proposal variance after every window of iterations.
 * @param inv_temperature_ladder_start The initial inverse temperature ladder for the parallel 
 * tempering sampler.
 * @param tune_ladder A boolean indicating whether to tune the temperature ladder during the
 * iterations.
 * @param target_temp_swap_accept_rate The target acceptance rate for the temperature swaps that the
 * ladder tuner will try to achieve by adjusting the temperature ladder.
 * @param target_temp_swap_accept_epsilon The epsilon threshold for checking whether the swap
 * acceptance rates have reached the target acceptance rate for the temperature swaps.
 * @param temp_window_size The size of the window for tuning the temperature ladder. The tuner
 * will adjust the temperature ladder after every window of iterations.
 * @param temp_window_size_growth_factor The growth factor for the temperature ladder tuning window
 * size. After each adjustment of the temperature ladder, the window size will be multiplied by
 * this growth factor.
 * @param temp_ladder_learning_rate The learning rate for adjusting the temperature ladder gaps.
 * @param iterations The total number of iterations to run the tuning process for.
 * @param stop_early A boolean indicating whether to stop the tuning process early if the target
 * acceptance rates for both the proposed gammas and the temperature swaps have been reached within
 * the specified epsilon thresholds.
 * @param seed The random seed to use for the GSL random number generator.
 * @param complete_swapping A boolean indicating whether to perform complete swapping of the chain
 * states during the temperature swaps (i.e., swapping both beta and gamma parameters) or only
 * swap the gamma parameters.
 * @param verbose The interval for printing progress updates to the console. If 0, no updates will
 * be printed.
 * 
 * @return A list containing the final proposal variances for each chain, the acceptance rates for 
 * the proposed gammas for each chain, the final inverse temperature ladder, the swap acceptance 
 * rates for each pair of adjacent temperatures, and the final iteration number reached during tuning.
 */
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List cpp_hprobit_tune_pt(
    const Rcpp::List& x_list,
    const Rcpp::List& y_list,
    const arma::colvec& mean_prior,
    const arma::mat& prec_prior,
    const arma::ivec& ncategories,
    const Rcpp::List& gamma_start,
    const arma::colvec& beta_start,
    const Rcpp::List& proposal_variance,
    const bool tune_proposal_variance,
    const double target_acceptance_rate,
    const double target_acceptance_epsilon,
    const int proposal_window_size,
    const arma::vec inv_temperature_ladder_start,
    const bool tune_ladder,
    const double target_temp_swap_accept_rate,
    const double target_temp_swap_accept_epsilon,
    const int temp_window_size,
    const double temp_window_size_growth_factor,
    const double temp_ladder_learning_rate,
    const int iterations,
    const bool stop_early,
    const int seed,
    const bool complete_swapping,
    const int verbose
){
    // Initialize GSL random number generator.
    gsl_rng_env_setup();                          // Read variable environnement
    const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
    gsl_rng* gen = gsl_rng_alloc(type);           // Rand generator allocation
    gsl_rng_set(gen, seed);

    // Unpack data and define constants.
    const Data data = unpack_data(x_list, y_list);
    const int ntargets = data.ntargets;
    const int ntemperatures = inv_temperature_ladder_start.size();
    const int nobs = data.nobs;
    const double min_gap = (1.0 - inv_temperature_ladder_start[ntemperatures-1]) 
        / (10 * (ntemperatures - 1)); // The min gap between temperatures.
    std::vector<colvec> gamma_start_vec = unpack_gamma(gamma_start, ncategories);

    // Create chains.
    PTChains chains(
        inv_temperature_ladder_start,
        beta_start,
        gamma_start_vec,
        unpack_proposal_variances(proposal_variance),
        mean_prior,
        prec_prior,
        ncategories,
        nobs,
        complete_swapping
    );

    // Initialize tuners.
    PTProposalTuner proposal_tuner(
        target_acceptance_rate,
        target_acceptance_epsilon,
        proposal_window_size,
        ntargets,
        ntemperatures
    );
    TemperatureLadderTuner ladder_tuner(
        target_temp_swap_accept_rate,
        target_temp_swap_accept_epsilon,
        temp_window_size,
        temp_window_size_growth_factor,
        temp_ladder_learning_rate,
        min_gap,
        ntemperatures
    );
    
    // Tuning loop.
    int iter = 0;
    for (; iter < iterations; iter++) {
        chains.simulate_step(data, gen);
        
        // Tune the proposal variance.
        if (tune_proposal_variance) {
            proposal_tuner.tune_step(chains, data, gen);
        }
        
        // Tune the temperature ladder.
        if (tune_ladder) {
            ladder_tuner.tune_step(chains, data, gen);
        }

        // Print progress.
        if (verbose > 0 && iter % verbose == 0) {
            Rcpp::Rcout << "Tuning iteration " << (iter+1) << " of " << iterations << " ";
        }

        // Check for early stopping.
        if (stop_early && proposal_tuner.has_reached_target() && ladder_tuner.has_reached_target()) {
            if (verbose > 0) {
                Rcpp::Rcout << "Early stopping at iteration " << (iter+1) 
                            << " as target acceptance rates have been reached." << std::endl;
            }
            break;
        }
    }

    // Free pointer.
    gsl_rng_free(gen);

    return Rcpp::List::create(
        _["proposal_variance"] = chains.get_proposal_variance(),
        _["proposal_acceptance_rates"] = proposal_tuner.get_acceptance_rates(),
        _["inv_temperature_ladder"] = chains.get_inv_temperatures(),
        _["temp_swap_rates"] = ladder_tuner.swap_rates,
        _["final_iteration"] = iter
    );
}