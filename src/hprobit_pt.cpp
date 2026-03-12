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
 * Perform the burn-in phase of the sampler, without adapting the temperature ladder.
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
 * Run the burn-in phase of the sampler, which includes adapting the temperature ladder.
 * 
 * @param chains The vector of TemperatureChain objects representing the Markov chains at each 
 * temperature. The inverse temperatures (betas) in these chains will be updated based on the new 
 * ladder.
 * @param ntemperatures The number of temperatures (i.e., the size of the chains vector).
 * @param data The data object.
 * @param burnin The number of burn-in iterations to perform.
 * @param rng The GSL random number generator to use for sampling.
 * @param window_size The initial window size for computing the swap acceptance ratio during the
 * burn-in. This controls how frequently the temperature ladder is updated during the burn-in phase
 * to achieve the target swap acceptance ratio.
 * @param window_growth_factor The growth factor for the window size during the burn-in. This 
 * controls how the window size changes over time during the burn-in phase. A value greater than 
 * 1 will cause the window size to grow over time, which can help stabilize the temperature ladder
 * adaptation as more samples are collected.
 * @param target_swap_ratio The target swap acceptance ratio to achieve for each pair of adjacent
 * temperatures. This is used to compute the adjustment to the ladder gaps during the burn-in phase
 * to achieve the desired swap acceptance ratio.
 * @param ladder_adjust_learning_rate The learning rate for adjusting the ladder gaps. This controls
 * how aggressively the ladder is adjusted based on the difference between the observed swap rates
 * and the target swap ratio.
 * @param swap_rates A vector where the computed swap rates for each pair of adjacent temperatures 
 * will be stored. This should have a length of ntemperatures - 1, and is used to track the swap rates
 * during the burn-in phase.
 * @param save_burning_samples Whether to store the samples collected during the burn-in phase.
 * @param thin The thinning interval for storing burn-in samples.
 */
// void do_adaptive_burnin(
//     std::vector<TemperatureChain>& chains, 
//     int ntemperatures,
//     const Data& data,
//     int burnin, 
//     gsl_rng* rng,
//     int window_size,
//     double window_growth_factor,
//     double target_swap_ratio,
//     double ladder_adjust_learning_rate,
//     arma::vec& swap_rates,
//     bool save_burning_samples,
//     int thin,
//     int verbose
// ) {
//     std::vector<double> swap_probabilities(ntemperatures, 0);
//     std::vector<int> nswap_accepts(ntemperatures-1, 0);
//     std::vector<int> nswap_proposals(ntemperatures-1, 0);
//     int window_step = 0;
//     double min_gap = (1.0 - chains[ntemperatures-1].inv_temperature) / ((ntemperatures - 1) * 0.1);
    
//     // Do burnin iterations.
//     for (unsigned int iter = 0; iter < burnin; iter++) {
//         do_step(chains, ntemperatures, data, nswap_accepts, nswap_proposals, 
//                 &swap_probabilities, rng);
        
//         // Store samples.
//         if (save_burning_samples && (iter % thin == 0)) {
//             for (auto& chain : chains) {
//                 chain.store_burnin_sample();
//             }
//         }

//         // Do ladder adaptation at the end of each window.
//         window_step++;
//         if (window_step == window_size) {
            
//             adjust_ladder(
//                 chains, 
//                 ntemperatures, 
//                 swap_probabilities, 
//                 swap_rates, 
//                 nswap_proposals, 
//                 target_swap_ratio, 
//                 ladder_adjust_learning_rate, 
//                 min_gap
//             );

//             swap_probabilities.assign(ntemperatures, 0);
//             nswap_proposals.assign(ntemperatures-1, 0);
//             nswap_accepts.assign(ntemperatures-1, 0);
//             window_step = 0;
//             window_size *= window_growth_factor;
//         }

//         // Print progress.
//         if (verbose > 0 && iter > 0 && (iter % verbose == 0 || iter == burnin - 1)) {
//             double mean_swap_ratio = 0;
//             for (int i = 0; i < ntemperatures - 1; i++) {
//                 mean_swap_ratio += swap_probabilities[i] 
//                     / (nswap_proposals[i] == 0 ? 1 : nswap_proposals[i]);
//             }
//             mean_swap_ratio /= (ntemperatures - 1);

//             Rcpp::Rcout << "Burn-in iteration " << (iter+1) << "/" << burnin 
//                         << ", Mean swap acceptance ratio: " << mean_swap_ratio
//                         << std::endl;
//          }
//     }
// }

/**
 * Adjust the temperature ladder based on the accumulated swap probabilities for each pair of
 * adjacent temperatures.
 * 
 * @param chains The vector of TemperatureChain objects representing the Markov chains at each 
 * temperature. The inverse temperatures (betas) in these chains will be updated based on the 
 * new ladder.
 * @param ntemperatures The number of temperatures (i.e., the size of the chains vector).
 * @param swap_probabilities The vector of accumulated swap probabilities for each pair of adjacent
 * temperatures. This should have a length of ntemperatures - 1, where each element corresponds to 
 * the swap probability for the pair of chains at temperatures i and i+1.
 * @param swap_rates A vector where the computed swap rates for each pair of adjacent temperatures
 * will be stored. This should have a length of ntemperatures - 1.
 * @param nswap_proposals The total number of swap proposals that have been made for each pair of
 * adjacent temperatures. This should have a length of ntemperatures - 1, and is used to compute the
 * swap rates from the accumulated swap probabilities.
 * @param target_swap_ratio The target swap acceptance ratio that we want to achieve for each pair
 * of adjacent temperatures. This is used to compute the adjustment to the ladder gaps.
 * @param ladder_adjust_learning_rate The learning rate for adjusting the ladder gaps. This controls
 * how aggressively the ladder is adjusted based on the difference between the observed swap rates
 * and the target swap ratio.
 * @param min_gap The minimum allowed gap between inverse temperatures in the ladder. This is used
 * to prevent the ladder from collapsing and to ensure that there is sufficient separation between 
 * the temperatures for effective parallel tempering. The gap between inverse temperatures i and 
 * i+1 is defined as inv_temperature[i] - inv_temperature[i+1], and this function ensures that 
 * this gap does not become smaller than min_gap for any pair of adjacent temperatures after the 
 * adjustment.
 */
// void adjust_ladder(
//     PTChains& chains, 
//     int ntemperatures,
//     std::vector<double>& swap_probabilities,
//     arma::vec& swap_rates,
//     std::vector<int>& nswap_proposals,
//     double target_swap_ratio,
//     double ladder_adjust_learning_rate,
//     double min_gap
// ) {
//     std::vector<double> ladder_gaps(ntemperatures - 1);

//     // Compute old gaps.
//     for(int i = 0; i < ntemperatures - 1; i++) {
//         ladder_gaps[i] = chains.chains[i].inv_temperature - chains.chains[i+1].inv_temperature;
//     }

//     // Compute mean swap rate for each pair.
//     for (int i = 0; i < ntemperatures - 1; i++) {
//         swap_rates[i] = swap_probabilities[i] / (nswap_proposals[i] == 0 ? 1 : nswap_proposals[i]);
//     }

//     // Update ladder gaps.
//     for (int i = 0; i < ntemperatures - 1; i++) {
//         ladder_gaps[i] *= std::exp(ladder_adjust_learning_rate * 
//             (swap_rates[i] - target_swap_ratio));

//         // Impose min.
//         if (ladder_gaps[i] < min_gap) {
//             ladder_gaps[i] = min_gap;
//         }
//     }

//     // Normalize the gaps to ensure beta_K is fixed.
//     double total_span = ladder_gaps[0];
//     for (int i = 1; i < ntemperatures - 1; i++) {
//         total_span += ladder_gaps[i];
//     }
//     double normalization_factor = (1.0 - chains.chains[ntemperatures-1].inv_temperature) / total_span;
//     for (int i = 0; i < ntemperatures - 1; i++) {
//         ladder_gaps[i] *= normalization_factor;
//     }

//     // Reconstruct ladder.
//     chains.chains[0].inv_temperature = 1.0;
//     for (int i = 1; i < ntemperatures-1; i++) {
//         chains.chains[i].inv_temperature = chains.chains[i-1].inv_temperature - ladder_gaps[i-1];
//     }
// }


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
    const int ntargets = data.ntargets;
    const int ntemperatures = inv_temperature_ladder.n_elem;

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
 * Sample a probit model using parallel tempering.
 * 
 * This implementation uses Constant Swap-Acceptance Adaptation to determine the temperature 
 * ladder during the burn-in.
 * 
 * @param xlist The feature matrices for each target. Each element of the list should be a matrix 
 * of size (n_samples, n_features).
 * @param ylist The response vectors for each target. Each element of the list should be a vector
 * of size (n_samples,).
 * @param mean_prior The prior mean for the regression coefficients.
 * @param prec_prior The prior precision matrix for the regression coefficients.
 * @param fix_zero Whether to fix the first threshold to zero.
 * @param ncategories A vector containing the number of categories for each target.
 * @param gamma_start The starting values for the thresholds for each target. Each element of the 
 * list should be a vector of size (n_categories - 1,).
 * @param beta_start The starting values for the regression coefficients.
 * @param tune The tuning parameters for the Metropolis-Hastings updates for the thresholds. This 
 * list can either contain one vector of tuning parameters, or a one vector for each temperature. 
 * Each vector should be of size (n_categories - 1,).
 * @param ntemperatures The number of temperatures to use in the parallel tempering algorithm.
 * @param temperature_ladder The ladder of temperatures to use in the parallel tempering algorithm. 
 * This should be a vector of size (ntemperatures,).
 * @param target_temp_swap_accept_ratio The target swap acceptance ratio to achieve during the 
 * burn-in. Set this to -1 to disable temperature ladder adaptation and use the provided temperature 
 * ladder as is.
 * @param temp_window_size The window size for computing the swap acceptance ratio during the 
 * burn-in. This controls how frequently the temperature ladder is updated during the burn-in phase 
 * to achieve the target swap acceptance ratio.
 * @param temp_window_size_growth_factor The growth factor for the window size during the burn-in. 
 * This controls how the window size changes over time during the burn-in phase. A value greater 
 * than 1 will cause the window size to grow over time, which can help stabilize the temperature 
 * ladder adaptation as more samples are collected.
 * @param temp_ladder_learning_rate The learning rate for adjusting the temperature ladder during 
 * the burn-in. This controls how aggressively the temperature ladder is updated to achieve the 
 * target swap acceptance ratio.
 * @param iterations The total number of iterations to run the sampler for (excluding burn-in).
 * @param burnin The number of burn-in iterations to discard.
 * @param thin The thinning interval for storing samples.
 * @param seed The random seed for reproducibility.
 * @param complete_swapping Whether to perform complete swapping (i.e., swapping both beta and gamma 
 * parameters) or partial swapping (i.e., swapping only the gamma parameters).
 * @param verbose The frequency (in iterations) at which to print progress updates. Set to 0 to 
 * disable.
 * 
 * @return A list containing the stored samples for the regression coefficients and thresholds for 
 * each target. The regression coefficients will be stored in a matrix of size (n_samples, 
 * n_features) for each temperature, and the thresholds will be stored in a list of matrices, 
 * where each matrix is of size (n_samples, n_categories - 1) for each temperature.
 */
// Rcpp::List cpp_hprobit_pt(
//     const Rcpp::List& xlist,
//     const Rcpp::List& ylist,
//     const arma::colvec& mean_prior,
//     const arma::mat& prec_prior,
//     const arma::ivec& ncategories,
//     const Rcpp::List& gamma_start,
//     const arma::colvec& beta_start,
//     const Rcpp::List& tune,
//     const int ntemperatures,
//     const arma::vec temperature_ladder,
//     const double target_temp_swap_accept_ratio,
//     const int temp_window_size,
//     const double temp_window_size_growth_factor,
//     const double temp_ladder_learning_rate,
//     const int iterations,
//     const int burnin,
//     const int thin,
//     const int seed,
//     const bool complete_swapping,
//     const bool save_burning_samples,
//     const int verbose
// ) {
//     if (verbose != 0){
//         Rcpp::Rcout << "Starting parallel tempering sampler for multi-scale probit model..." << std::endl;
//     }
//     // Initialize GSL random number generator.
//     gsl_rng_env_setup();                          // Read variable environnement
//     const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
//     gsl_rng* gen = gsl_rng_alloc(type);           // Rand generator allocation
//     gsl_rng_set(gen, seed);

//     // Define constants and unpack data.
//     const int total_iterations = iterations + burnin;
//     const int nstore = iterations / thin;
//     const int nstore_bunin =  save_burning_samples ? burnin / thin : 0;
//     const int ntargets = xlist.size();
//     const Data data = unpack_data(xlist, ylist);

//     const int npredictors = data.Xall.n_cols;
//     const int nobs = data.Xall.n_rows;

//     const std::vector<colvec> gamma_start_vec = unpack_gamma(gamma_start, ncategories);

//     if (verbose != 0){
//         Rcpp::Rcout << "Data unpacked. Starting parallel tempering sampler with " << ntemperatures << " temperatures." << std::endl;
//     }

//     // Create chains.
//     std::vector<TemperatureChain> chains = create_temperature_chains(ntemperatures, 
//         temperature_ladder, beta_start, gamma_start_vec, ncategories, tune, mean_prior, prec_prior,
//         nstore, nobs, ntargets, save_burning_samples, nstore_bunin);

//     if (verbose != 0) {
//         Rcpp::Rcout << "Starting burn-in phase..." << std::endl;
//     }

//     arma::vec adaptation_swap_rates(ntemperatures-1, arma::fill::zeros);
//     // Measure burnin sampling time.
//     auto start_time_burnin = std::chrono::high_resolution_clock::now();
//     // Burn-in loop without temperature ladder adaptation.
//     if (target_temp_swap_accept_ratio == -1) {
//         do_burnin(chains, ntemperatures, data, burnin, save_burning_samples, thin, gen);
//     }
//     else { // Do with adaptation.
//         do_adaptive_burnin(
//             chains, 
//             ntemperatures, 
//             data, 
//             burnin, 
//             gen, 
//             temp_window_size, 
//             temp_window_size_growth_factor, 
//             target_temp_swap_accept_ratio,
//             temp_ladder_learning_rate,
//             adaptation_swap_rates,
//             save_burning_samples,
//             thin,
//             verbose
//         );
//     }
//     // Measure burnin time in seconds.
//     auto end_time_burnin = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed_burnin = end_time_burnin - start_time_burnin;
    
//     if (verbose != 0) {
//         Rcpp::Rcout << "Burn-in complete. Starting sampling phase..." << std::endl;
//     }

//     // Sampling loop.
//     std::vector<int> nswap_accepts(ntemperatures-1, 0);
//     std::vector<int> nswap_proposals(ntemperatures-1, 0);
//     std::vector<double> swap_probabilities(ntemperatures-1, 0);
//     // Measure sampling time.
//     auto start_time = std::chrono::high_resolution_clock::now();
//     for (unsigned int iter = 0; iter < iterations; iter++) {
//         do_step(chains, ntemperatures, data, nswap_accepts, nswap_proposals, 
//             &swap_probabilities, gen);

//         // Store samples.
//         if (iter % thin == 0) {
//             for (auto& chain : chains) {
//                 chain.store_sample();
//             }
//         }

//         // Print progress updates.
//         if (verbose > 0 && iter > 0 && (iter % verbose == 0 || iter == iterations - 1)) {

//             // Compute mean swap ratio.
//             double mean_swap_ratio = 0;
//             for (int i = 0; i < ntemperatures - 1; i++) {
//                 mean_swap_ratio += swap_probabilities[i] 
//                     / (nswap_proposals[i] == 0 ? 1 : nswap_proposals[i]);
//             }
//             mean_swap_ratio /= (ntemperatures - 1);

//             // Print progress.
//             Rcpp::Rcout << "Iteration " << (iter+1) << "/" << iterations 
//                         << ", Mean swap acceptance ratio: " << mean_swap_ratio
//                         << std::endl;
//         }
//     }
//     // Measure sampling time in seconds.
//     auto end_time = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed = end_time - start_time;

//     // Free pointers.
//     gsl_rng_free(gen);

//     // Wrap gammas into Rcpp::List for output. This prevents R from automatically converting the 
//     // matrices to vectors, which causes issues when we have multiple targets.
//     Rcpp::List storegamma_list;
//     for (size_t i = 0; i < chains[0].store_gamma.size(); ++i) {
//         storegamma_list.push_back(chains[0].store_gamma[i]);
//     }

//     // Wrap burnin gammas too.
//     Rcpp::List storegamma_burnin_list;
//     if (save_burning_samples) {
//         for (size_t i = 0; i < chains[0].store_gamma_burnin.size(); ++i) {
//             storegamma_burnin_list.push_back(chains[0].store_gamma_burnin[i]);
//         }
//     }

//     // Store inverse temperatures.
//     arma::vec inv_temps(ntemperatures);
//     arma::vec adapted_temps(ntemperatures);
//     for (int i = 0; i < ntemperatures; i++) {
//         inv_temps(i) = chains[i].inv_temperature;
//         adapted_temps(i) = 1.0 / chains[i].inv_temperature;
//     }

//     // Return results.
//     return Rcpp::List::create(
//         _["storebeta"] = chains[0].store_beta,
//         _["storegamma"] = storegamma_list,
//         _["nswap_accepts"] = nswap_accepts,
//         _["nswap_proposals"] = nswap_proposals,
//         _["adapted_inv_temps"] = inv_temps,
//         _["adapted_temps"] = adapted_temps,
//         _["adaptation_swap_rates"] = adaptation_swap_rates,
//         _["storebeta_burnin"] = chains[0].store_beta_burnin,
//         _["storegamma_burnin"] = storegamma_burnin_list,
//         _["sampling_time"] = elapsed.count(),
//         _["burnin_time"] = elapsed_burnin.count()
//     );
// }



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