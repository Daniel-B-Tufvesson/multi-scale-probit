/**
 * @file hprobit_pt.cpp
 * @author Daniel Tufvesson
 * @version 1.0
 * @brief A Metropolis-Hastings-withing-Gibbs sampler for the Multi-Scale Probit model, 
 * that uses parallel tempering to improve mixing.
 */

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "mspm_util.hpp"
#include "rtnorm.hpp"

/**
 * Struct to hold the data for the model. This includes the feature matrices and response vectors
 * for each target, as well as the combined feature matrix and precomputed X'X matrix for all 
 * targets.
 */
struct Data {
    std::vector<arma::mat> X;
    std::vector<arma::colvec> Y;
    arma::mat Xall;
    arma::mat XpX;
    std::vector<int> target_nobs;
};

/**
 * Class representing a single Markov chain at a given temperature in the parallel tempering 
 * algorithm. It supports both complete and partial swapping of the chain state.
 */
class TemperatureChain {
public:
    double inv_temperature; // 1/T, also known as beta in parallel tempering literature.
    arma::colvec beta;
    std::vector<arma::colvec> gamma;
    arma::colvec ystar;
    arma::ivec ncategories;
    arma::vec gamma_tune;
    arma::colvec mean_prior;
    arma::mat prec_prior;
    int npredictors;
    int ntargets;

    // Storage matrixes.
    arma::mat store_beta;
    std::vector<arma::mat> store_gamma;
    int nstored = 0; // Number of samples stored so far.
    int nsteps = 0; // Number of steps taken.

    /**
     * Constructor for the TemperatureChain class. Initializes the chain with the given starting 
     * values and sets up the storage matrixes for the samples.
     * 
     * @param inv_temperature The inverse temperature (1/T) for this chain. Higher values correspond 
     * to lower temperatures. This is also known as beta in the parallel tempering literature.
     * @param beta_start The starting values for the regression coefficients.
     * @param gamma_start The starting values for the threshold parameters for each target.
     * @param ncategories The number of categories for each target, which determines the number of
     * thresholds.
     * @param gamma_tune The tuning parameters for the proposal distribution for the thresholds. 
     * This controls the standard deviation of the truncated normal distribution used for proposing 
     * new threshold values.
     * @param mean_prior The prior mean for the regression coefficients.
     * @param prec_prior The prior precision matrix for the regression coefficients.
     * @param nstore The number of samples to store for this chain.
     * @param total_nobs The total number of observations across all targets. This is used to set 
     * the size of the ystar vector and the store_beta matrix.
     * @param ntargets The number of targets (i.e., the number of different response variables). 
     */
    TemperatureChain(
        double inv_temperature, 
        const arma::colvec beta_start, 
        const std::vector<arma::colvec> gamma_start,
        const arma::ivec& ncategories, 
        const arma::vec& gamma_tune,
        const arma::colvec& mean_prior,
        const arma::mat& prec_prior,
        const int nstore,
        const int total_nobs,
        const int ntargets
    ) : inv_temperature(inv_temperature), beta(beta_start), gamma(gamma_start), 
        ncategories(ncategories), gamma_tune(gamma_tune), mean_prior(mean_prior), 
        prec_prior(prec_prior), ntargets(ntargets) {
        
        
        // Initialize storage matrixes.
        store_beta = arma::mat(nstore, beta.n_rows, arma::fill::zeros);
        store_gamma = std::vector<arma::mat>(ntargets);
        for (unsigned int target = 0; target < ntargets; target++) {
            store_gamma[target] = arma::mat(nstore, ncategories(target)-1, arma::fill::zeros);
        }

        npredictors = beta.n_elem;
        ystar = arma::colvec(total_nobs, arma::fill::zeros);

        // Set extremes for gamma values. Note that this always overrides the start values for
        // the first and last gammas. Is this intentional in the original implementation?
        for (unsigned int target = 0; target < ntargets; target++) {
            gamma[target](0) = -std::numeric_limits<double>::max();
            gamma[target](ncategories[target]) = std::numeric_limits<double>::max();
        }
    }

    /**
     * Partially swap the state of this chain with another chain. This involves swapping only the 
     * gamma parameters between the two chains, while keeping the beta parameters unchanged.
     * Note that only the latest values of the parameters are swapped, and the stored samples in 
     * the storage matrixes are not swapped.
     * 
     * @param other_chain The other TemperatureChain instance with which to swap the gamma 
     * parameters.
     */
    void swap_gammas(TemperatureChain& other_chain) {
        for (unsigned int target = 0; target < gamma.size(); target++) {
            for (unsigned int j = 0; j < ncategories(target); j++) {
                double temp = gamma[target][j];
                gamma[target][j] = other_chain.gamma[target][j];
                other_chain.gamma[target][j] = temp;
            }
        }
    }

    /**
     * Partially swap the state of this chain with another chain. This involves swapping only the 
     * beta parameters between the two chains, while keeping the gamma parameters unchanged. Note
     * that only the latest values of the parameters are swapped, and the stored samples in the 
     * storage matrixes are not swapped.
     * 
     * @param other_chain The other TemperatureChain instance with which to swap the beta 
     * parameters.
     */
    void swap_beta (TemperatureChain& other_chain) {
        for (unsigned int j = 0; j < npredictors; j++) {
            double temp = beta[j];
            beta[j] = other_chain.beta[j];
            other_chain.beta[j] = temp;
        }
    }

    /**
     * Determine whether to swap the state of this chain with another chain based on the computed
     * acceptance ratio for the swap. The acceptance ratio is computed based on the likelihood of 
     * the data given the current parameter values in each chain, and the difference in inverse
     * temperatures between the two chains.
     * 
     * @param other_chain The other TemperatureChain instance with which to potentially swap states.
     * @param data The data object containing the feature matrices and response vectors for each 
     * target.
     * @param rng The GSL random number generator to use for sampling the uniform random variable 
     * for the acceptance step.
     */
    bool should_swap(
        const TemperatureChain& other_chain, 
        const Data& data,
        gsl_rng* rng
    ) {
        // Compute the log acceptance ratio for the swap.
        auto cdf = gsl_cdf_ugaussian_P;
        double log_swap_accept_ratio = 0;
        for (unsigned int target = 0; target < ntargets; target++) {
            const colvec ystar1 = data.X[target] * beta;
            const colvec ystar2 = data.X[target] * other_chain.beta;

            for (unsigned int i = 0; i < data.X[target].n_rows; i++) {
                log_swap_accept_ratio = log_swap_accept_ratio
                    + log(cdf(gamma[target](data.Y[target](i)) - ystar1[i]) - 
                          cdf(gamma[target](data.Y[target](i)-1) - ystar1[i]))
                    - log(cdf(gamma[target](data.Y[target](i)) - ystar2[i]) - 
                          cdf(gamma[target](data.Y[target](i)-1) - ystar2[i]));
            }
        }

        // Scale it by inv_temperature delta.
        log_swap_accept_ratio = log_swap_accept_ratio * (inv_temperature 
            - other_chain.inv_temperature);
        
        // Do the swap with the computed acceptance ratio.
        return gsl_ran_flat(rng, 0.0, 1.0) <= exp(log_swap_accept_ratio);
    }

    /**
     * Store the latest beta and gammas samples.
     */
    void storeSample() {
        // Store beta.
        for (unsigned int j = 0; j < npredictors; j++) {
            store_beta(nstored, j) = beta[j];
        }

        // Store gammas.
        for (unsigned int target = 0; target < ntargets; target++) {
            for (unsigned int j = 1; j < ncategories(target); j++) {
                store_gamma[target](nstored, j-1) = gamma[target](j);
            }
        }

        nstored++;
    }

    /**
     * Perform one step of the within-temperature Metropolis-Hastings updates for the thresholds 
     * and Gibbs update for the regression coefficients.
     * 
     * @param data The data object.
     * @param Xall The combined feature matrix for all targets.
     * @param rng The GSL random number generator to use for sampling.
     */
    void do_step(
        const Data& data,
        gsl_rng* rng
    ) {
        step_gamma(data, rng);
        step_beta(data, rng);
        nsteps++;
    }

private:

    /**
     * Perform the Metropolis-Hastings updates for the thresholds for each target. This includes
     * proposing new threshold values and doing the acceptance step.
     * 
     * @param data The data object.
     * @param rng The GSL random number generator to use for sampling.
     */
    void step_gamma(
        const Data& data,
        gsl_rng* rng
    ) {
        for (unsigned int t = 0; t < ntargets; t++) {
            int target = (nsteps+t) % ntargets;
            int ncats = ncategories(target); 

            // Make proposal.
            arma::colvec gamma_prop = propose_gamma(target, ncats, rng);
            double log_accept_ratio = compute_gamma_log_accept_ratio(data, gamma_prop, target);

            // Do MH acceptance step.
            if (gsl_ran_flat(rng, 0.0, 1.0) <= exp(log_accept_ratio)) {
                gamma[target] = gamma_prop;
            }
        }
    }

    /**
     * Propose new threshold values for a given target using a truncated normal distribution. 
     * 
     * @param target The index of the target for which to propose new threshold values.
     * @param ncats The number of categories for the target.
     * @param rng The GSL random number generator to use for sampling.
     * 
     * @return A vector containing the proposed new threshold values for the target.
     */
    arma::colvec propose_gamma(int target, int ncats, gsl_rng* rng) {
        arma::colvec gamma_prop(ncats + 1, arma::fill::zeros);
        gamma_prop.head(1) = -INFINITY;
        gamma_prop.tail(1) = INFINITY;

        if (ncats == 2) {
            // Draw new split point
            gamma_prop(1) = rtnorm(
                rng,
                -INFINITY,
                INFINITY,
                gamma[target](1),
                gamma_tune(target)
            ).first;
        } 
        else {
            for (unsigned int i = 1; i < ncats; i++) {
                if (i == 1) { // If first gamma
                    gamma_prop(i) = rtnorm(
                        rng,
                        gamma[target](i-1),
                        gamma[target](i+1),
                        gamma[target](i),
                        gamma_tune(target)
                    ).first;
                } 
                else { // If any other gamma
                    gamma_prop(i) = rtnorm(
                        rng,
                        gamma_prop(i-1),
                        gamma[target](i+1),
                        gamma[target](i),
                        gamma_tune(target)
                    ).first;
                }
            }
        }
        return gamma_prop;
    }

    /**
     * Compute the (log) likelihood acceptance ratio for the new proposal. The proposal is compared
     * to the last accepted value of the thresholds, and the likelihood is computed based on the 
     * current values of the regression coefficients and the data.
     * 
     * @param data The data object.
     * @param gamma_prop The proposed new threshold values for the target.
     * @param target The index of the target for which to compute the acceptance ratio.
     * 
     * @return The log acceptance ratio for the proposed new threshold values.
     */
    double compute_gamma_log_accept_ratio(
        const Data& data,
        const arma::colvec& gamma_prop,
        int target
    ) {
        double log_likelihood_ratio = 0;
        colvec y_star = data.X[target] * beta;
        auto cdf = gsl_cdf_ugaussian_P;

        for (unsigned int i = 0; i < data.X[target].n_rows; i++) {
            int y_val = data.Y[target](i);

            // Handle last category.
            if (y_val == ncategories(target)) {
                log_likelihood_ratio = log_likelihood_ratio
                    + log(1.0 - cdf(gamma_prop(y_val-1) - y_star[i]))
                    - log(1.0 - cdf(gamma[target](y_val-1) - y_star[i]));
            }
            // Handle first category.
            else if (y_val == 1) {
                log_likelihood_ratio = log_likelihood_ratio
                    + log(cdf(gamma_prop(y_val) - y_star[i]))
                    - log(cdf(gamma[target](y_val) - y_star[i]));
            }
            // Handle categories inbetween.
            else {
                log_likelihood_ratio = log_likelihood_ratio
                    + log(cdf(gamma_prop(y_val) - y_star[i]) -
                          cdf(gamma_prop(y_val-1) - y_star[i]))
                    - log(cdf(gamma[target](y_val) - y_star[i]) - 
                          cdf(gamma[target](y_val-1) - y_star[i]));
            }
        }
        return log_likelihood_ratio;
    }

    /**
     * Perform the Gibbs update for the regression coefficients. This involves first updating the
     * latent variables y* based on the current thresholds and regression coefficients, and then
     * drawing new regression coefficients from their full conditional distribution given the 
     * updated latent variables and thresholds.
     * 
     * @param data The data object.
     * @param rng The GSL random number generator to use for sampling.
     */
    void step_beta(
        const Data& data,
        gsl_rng* rng
    ) {
        // First update latent y*.
        int offset = 0;
        for (unsigned int target = 0; target < ntargets; target++) {
            if (target > 0) {
                int nobs = data.Y[target-1].n_elem;
                offset += nobs;
            }
            const colvec current_ystar = data.X[target] * beta;
            for (unsigned int i = 0; i < data.X[target].n_rows; i++) {
                int idx1 = data.Y[target](i) - 1;
                int idx2 = data.Y[target](i);
                ystar(offset + i) = rtnorm(
                    rng,
                    gamma[target](data.Y[target](i) - 1),
                    gamma[target](data.Y[target](i)),
                    current_ystar[i], 
                    1.0
                ).first;
            }
        }

        // Then update beta.
        arma::mat XpZ = arma::trans(data.Xall) * ystar;
        beta = mspm_util::NormNormregress_beta_draw(rng, data.XpX, XpZ, mean_prior, prec_prior, 1.0);
    }
};

// Function declarations. -----------------------------------------------------------------------

Data unpack_data(
    const Rcpp::List& xlist,
    const Rcpp::List& ylist
);

void do_step(
    std::vector<TemperatureChain>& chains, 
    int ntemperatures,
    const Data& data,
    int& nswap_accepts,
    int& nswap_proposals,
    gsl_rng* rng
);


// Function definitions. -------------------------------------------------------------------------

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
 * @param iterations The total number of iterations to run the sampler for (excluding burn-in).
 * @param burnin The number of burn-in iterations to discard.
 * @param thin The thinning interval for storing samples.
 * @param seed The random seed for reproducibility.
 * @param verbose The frequency (in iterations) at which to print progress updates. Set to 0 to 
 * disable.
 * 
 * @return A list containing the stored samples for the regression coefficients and thresholds for 
 * each target. The regression coefficients will be stored in a matrix of size (n_samples, 
 * n_features) for each temperature, and the thresholds will be stored in a list of matrices, 
 * where each matrix is of size (n_samples, n_categories - 1) for each temperature.
 */
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List cpp_hprobit_pt(
    const Rcpp::List& xlist,
    const Rcpp::List& ylist,
    const arma::colvec& mean_prior,
    const arma::mat& prec_prior,
    const int fix_zero,
    const arma::ivec& ncategories,
    const Rcpp::List& gamma_start,
    const arma::colvec& beta_start,
    const Rcpp::List& tune,
    const int ntemperatures,
    const int iterations,
    const int burnin,
    const int thin,
    const int seed,
    const int verbose
) {
    // Initialize GSL random number generator.
    gsl_rng_env_setup();                          // Read variable environnement
    const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
    gsl_rng* gen = gsl_rng_alloc(type);           // Rand generator allocation
    gsl_rng_set(gen, seed);

    // Define constants and unpack data.
    const int total_iterations = iterations + burnin;
    const int nstore = iterations / thin;
    const int ntargets = xlist.size();
    const Data data = unpack_data(xlist, ylist);

    const int npredictors = data.Xall.n_cols;
    const int nobs = data.Xall.n_rows;

    // Unpack gammas.
    std::vector<arma::colvec> gamma_start_vec(ntargets);
    for (unsigned int target = 0; target < ntargets; target++) {
        gamma_start_vec[target] = Rcpp::as<arma::colvec>(gamma_start[target]);
    }

    // Create chains.
    std::vector<TemperatureChain> chains;
    for (int i = 0; i < ntemperatures; i++) {
        double inv_temperature = 1.0 / std::pow(2.0, i); // Example temperature ladder: T = 1, 2, 4, 8,
        chains.emplace_back(
            inv_temperature,
            beta_start,
            gamma_start_vec,
            ncategories,
            arma::vec(ntargets, arma::fill::ones), // Todo: we need to unpack the tuning parameters here.
            mean_prior,
            prec_prior,
            nstore,
            nobs,
            ntargets
        );
    }

    // Burn-in loop.
    int nswap_accepts = 0;
    int nswap_proposals = 0;
    for (unsigned int iter = 0; iter < burnin; iter++) {
        do_step(chains, ntemperatures, data, nswap_accepts, nswap_proposals, gen);

        // Todo: adjust the temperature ladder based on acceptance ratio.
    }

    // Sampling loop.
    nswap_accepts = 0;
    nswap_proposals = 0;
    for (unsigned int iter = 0; iter < iterations; iter++) {
        do_step(chains, ntemperatures, data, nswap_accepts, nswap_proposals, gen);

        // Store samples.
        if (iter % thin == 0) {
            for (auto& chain : chains) {
                chain.storeSample();
            }
        }

        // Print progress updates.
        if (verbose > 0 && (iter % verbose == 0 || iter == iterations - 1)) {
            Rcpp::Rcout << "Iteration " << (iter+1) << "/" << iterations 
                        << ", Swap acceptance ratio: " 
                        << (nswap_proposals > 0 ? static_cast<double>(nswap_accepts) / nswap_proposals : 0.0)
                        << std::endl;
        }
    }

    // Free pointers.
    gsl_rng_free(gen);

    // Wrap gammas into Rcpp::List for output. This prevents R from automatically converting the 
    // matrices to vectors, which causes issues when we have multiple targets.
    Rcpp::List storegamma_list;
    for (size_t i = 0; i < chains[0].store_gamma.size(); ++i) {
        storegamma_list.push_back(chains[0].store_gamma[i]);
    }

    // Return results.
    return Rcpp::List::create(
        _["storebeta"] = chains[0].store_beta,
        _["storegamma"] = storegamma_list,
        _["nswap_accepts"] = nswap_accepts,
        _["nswap_proposals"] = nswap_proposals
    );
}

/**
 * Helper function to perform one step of the parallel tempering sampler, including 
 * within-temperature Metropolis-Hastings updates and between-temperature swaps.
 * 
 * @param chains The vector of TemperatureChain objects representing the Markov chains at each 
 * temperature.
 * @param ntemperatures The number of temperatures (i.e., the size of the chains vector).
 * @param data The data object.
 * @param nswap_accepts A reference to an integer that counts the number of accepted swaps between
 * temperatures.
 * @param nswap_proposals A reference to an integer that counts the number of proposed swaps 
 * between temperatures.
 * @param rng The GSL random number generator to use for sampling.
 */
void do_step(
    std::vector<TemperatureChain>& chains, 
    int ntemperatures,
    const Data& data,
    int& nswap_accepts,
    int& nswap_proposals,
    gsl_rng* rng
) {
    // Do within-temperature MH updates.
    for (auto& chain : chains) {
        chain.do_step(data, rng);
    }

    // Select pair of adjacent temperatures to swap.
    int swap_index = rand() % (ntemperatures - 1); // fix: we need to make sure the seed is set.
    TemperatureChain& chain1 = chains[swap_index];
    TemperatureChain& chain2 = chains[swap_index + 1];
    nswap_proposals++;
    if (chain1.should_swap(chain2, data, rng)) {
        // Todo: check if we only want to do partial chain swapping.
        chain1.swap_gammas(chain2);
        chain1.swap_beta(chain2);
        nswap_accepts++;
    }
}

/**
 * Unpack the data from the input lists into a Data object containing the feature matrices, 
 * response vectors, and precomputed cross-product of the combined feature matrix. This 
 * function assumes that the input lists are properly formatted and contain the expected types 
 * and dimensions. It also computes the combined feature matrix for all targets and the 
 * cross-product of this matrix, which will be used in the Gibbs updates for the regression 
 * coefficients.
 * 
 * @param xlist The list of feature matrices for each target.
 * @param ylist The list of response vectors for each target.
 * 
 * @return A Data object containing the unpacked feature matrices, response vectors, combined 
 * feature matrix, cross-product of the combined feature matrix, and the number of observations 
 * for each target.
 */
Data unpack_data(
    const Rcpp::List& xlist,
    const Rcpp::List& ylist
) {
    const int ntargets = xlist.size();
    Data data {
        .X = std::vector<arma::mat>(ntargets),
        .Y = std::vector<arma::colvec>(ntargets),
        .Xall = arma::mat(),
        .XpX = arma::mat(),
        .target_nobs = std::vector<int>(ntargets)
    };
    for (unsigned int target = 0; target < ntargets; ++target) {
        data.X[target] = Rcpp::as<arma::mat>(xlist[target]);
        data.Y[target] = Rcpp::as<arma::colvec>(ylist[target]);
        if (target == 0) {
            data.Xall = data.X[target];
        } 
        else {
            data.Xall = arma::join_vert(data.Xall, data.X[target]);
        }
        data.target_nobs[target] = data.Y[target].n_elem;
    }
    data.XpX = mspm_util::crossprod(data.Xall, data.Xall);
    return data;
}
