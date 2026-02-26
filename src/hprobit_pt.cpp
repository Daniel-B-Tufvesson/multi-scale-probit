/**
 * @file hprobit_pt.cpp
 * @author Daniel Tufvesson
 * @version 1.0
 * @brief A Metropolis-Hastings-withing-Gibbs sampler for the Multi-Scale Probit model, 
 * that uses parallel tempering to improve mixing.
 */

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "sampling.hpp"


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
    double compute_swap_probability(
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

            // Loop over all data points for target.
            for (unsigned int i = 0; i < data.X[target].n_rows; i++) {
                log_swap_accept_ratio = log_swap_accept_ratio
                    + log(cdf(gamma[target](data.Y[target](i)) - ystar2[i]) - 
                          cdf(gamma[target](data.Y[target](i)-1) - ystar2[i]))
                    - log(cdf(gamma[target](data.Y[target](i)) - ystar1[i]) - 
                          cdf(gamma[target](data.Y[target](i)-1) - ystar1[i]));
            }

        }

        // Scale it by inv_temperature delta.
        log_swap_accept_ratio *= (inv_temperature - other_chain.inv_temperature);
        
        // Do the swap with the computed acceptance ratio.
        return std::min(1.0, exp(log_swap_accept_ratio));
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
            double log_likelihood_ratio = compute_log_likelihood_ratio(data, gamma_prop, target);
            double log_proposal_ratio = compute_gamma_log_proposal_ratio(gamma_prop, target);
            double log_accept_ratio = inv_temperature * log_likelihood_ratio + log_proposal_ratio;

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
     * @return A column vector containing the proposed new threshold values for the given target. 
     * The first and last elements of the vector are set to -Inf and Inf, respectively, and the 
     * proposed values for the thresholds are in between these extremes.
     */
    arma::colvec propose_gamma(int target, int ncats, gsl_rng* rng) {
        arma::colvec gamma_prop(ncats + 1, arma::fill::zeros);
        gamma_prop.head(1) = -INFINITY;
        gamma_prop.tail(1) = INFINITY;
        // Note: probabilities are 0 for the edge thresholds.

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
                        gamma_prop(i-1),    // a
                        gamma[target](i+1), // b
                        gamma[target](i),   //mu
                        gamma_tune(target)  //sigma
                    ).first;
                }
            }
        }
        return gamma_prop;
    }

    /**
     * Compute the log probability ratio betweent the proposal distribution for the proposed gammas,
     * i.e. the log of q(ɣ|ɣ') / q(ɣ'|ɣ), where ɣ is the old gamma, ɣ' is the proposed gamma and 
     * q(.) is the truncated Gaussian proposal distribution.
     * 
     * Note that this is not the full proposal ratio, but only the part of the proposal ratio that
     * depends on the proposed gammas. The full proposal ratio also includes the probabilities of the
     * proposed gammas under the truncated normal distribution.
     * 
     * @param gamma_p The proposed new threshold values for the target.
     * @param target The index of the target for which to compute the log proposal ratio.
     * @return The log of the proposal ratio for the proposed gammas.
     */
    double compute_gamma_log_proposal_ratio(
        const colvec& gamma_p,
        int target
    ) {
        double log_proposal_ratio = 0;
        for (int i = 1; i < ncategories(target); i++) {
            // Normalization term Z(gamma) 
            double z_old = compute_log_trunc_gauss_norm_constant(
                gamma[target](i), 
                gamma[target](i-1), 
                gamma[target](i+1), 
                gamma_tune(target)
            );

            if (i == 1) { // If first gamma
                // Normalization term Z(gamma_p)
                double z_prop = compute_log_trunc_gauss_norm_constant(
                    gamma_p(i), 
                    gamma[target](i-1), 
                    gamma[target](i+1), 
                    gamma_tune(target)
                );
                
                log_proposal_ratio += z_old - z_prop;
            } 
            else { // If any other gamma
                // Normalization term Z(gamma_p)
                double z_prop = compute_log_trunc_gauss_norm_constant(
                    gamma_p(i), 
                    gamma_p(i-1), 
                    gamma[target](i+1), 
                    gamma_tune(target)
                );
                log_proposal_ratio += z_old - z_prop;
            }
        }
        return log_proposal_ratio;
    }

    /**
     * Compute the log of the normalization constant for the truncated Gaussian distribution used 
     * in the proposal distribution for the gammas.
     * 
     * @param x The value at which to compute the normalization constant.
     * @param a The lower truncation point for the truncated Gaussian distribution.
     * @param b The upper truncation point for the truncated Gaussian distribution.
     * @param sigma The standard deviation of the truncated Gaussian distribution.
     * @return The log of the normalization constant for the truncated Gaussian distribution at the 
     * given value x.
     */
    double compute_log_trunc_gauss_norm_constant(double x, double a, double b, double sigma) {
        auto cdf = gsl_cdf_ugaussian_P;
        double alpha = (a - x) / sigma;
        double beta = (b - x) / sigma;
        return log(cdf(beta) - cdf(alpha));
    }

    /**
     * Compute the unscaled log likelihood ratio between the current gammas and the proposed gammas.
     * Note that this is unscaled, i.e. it does not include the difference in inverse temperatures 
     * between the two chains.
     * 
     * @param data The data object.
     * @param gamma_p The proposed new threshold values for the target.
     * @param target The index of the target for which to compute the log likelihood ratio.
     * @return The unscaled log likelihood ratio between the current gammas and the proposed gammas 
     * for the given target.
     */
    double compute_log_likelihood_ratio(
        const Data& data,
        const colvec& gamma_p,
        int target
    ) {
        auto cdf = gsl_cdf_ugaussian_P;
        double log_likelihood_ratio = 0;
        colvec y_star = data.X[target] * beta;
        for (unsigned int i = 0; i < data.X[target].n_rows; i++) {
            int y_val = data.Y[target](i);

            // Handle last category.
            if (y_val == ncategories(target)) {
                log_likelihood_ratio = log_likelihood_ratio
                    + log(1.0 - cdf(gamma_p(y_val-1) - y_star[i]))
                    - log(1.0 - cdf(gamma[target](y_val-1) - y_star[i]));
            }
            // Handle first category.
            else if (y_val == 1) {
                log_likelihood_ratio = log_likelihood_ratio
                    + log(cdf(gamma_p(y_val) - y_star[i]))
                    - log(cdf(gamma[target](y_val) - y_star[i]));
            }
            // Handle categories inbetween.
            else {
                log_likelihood_ratio = log_likelihood_ratio
                    + log(cdf(gamma_p(y_val) - y_star[i]) -
                          cdf(gamma_p(y_val-1) - y_star[i]))
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



void do_step(
    std::vector<TemperatureChain>& chains, 
    int ntemperatures,
    const Data& data,
    std::vector<int>& nswap_accepts,
    std::vector<int>& nswap_proposals,
    std::vector<double>* swap_probabilities,
    gsl_rng* rng
);

void do_burnin(
    std::vector<TemperatureChain>& chains, 
    int ntemperatures,
    const Data& data,
    int burnin, 
    gsl_rng* rng
);

void do_adaptive_burnin(
    std::vector<TemperatureChain>& chains, 
    int ntemperatures,
    const Data& data,
    int burnin, 
    gsl_rng* rng,
    int window_size,
    double window_growth_factor,
    double target_swap_ratio,
    double ladder_adjust_learning_rate,
    arma::vec& swap_rates,
    int verbose
);

void adjust_ladder(
    std::vector<TemperatureChain>& chains, 
    int ntemperatures,
    std::vector<double>& swap_probabilities,
    arma::vec& swap_rates,
    std::vector<int>& nswap_proposals,
    double target_swap_ratio,
    double ladder_adjust_learning_rate,
    double min_gap
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
    const arma::ivec temperature_ladder,
    const double target_temp_swap_accept_ratio,
    const int temp_window_size,
    const double temp_window_size_growth_factor,
    const double temp_ladder_learning_rate,
    const int iterations,
    const int burnin,
    const int thin,
    const int seed,
    bool complete_swapping,
    const int verbose
) {
    if (verbose != 0){
        Rcpp::Rcout << "Starting parallel tempering sampler for multi-scale probit model..." << std::endl;
    }
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

    if (verbose != 0){
        Rcpp::Rcout << "Data unpacked. Starting parallel tempering sampler with " << ntemperatures << " temperatures." << std::endl;
    }

    // Create chains.
    std::vector<TemperatureChain> chains;
    for (int i = 0; i < ntemperatures; i++) {
        double inv_temperature = 1.0 / temperature_ladder(i);
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

    if (verbose != 0) {
        Rcpp::Rcout << "Starting burn-in phase..." << std::endl;
    }

    arma::vec adaptation_swap_rates(ntemperatures-1, arma::fill::zeros);
    // Burn-in loop without temperature ladder adaptation.
    if (target_temp_swap_accept_ratio == -1) {
        do_burnin(chains, ntemperatures, data, burnin, gen);
    }
    else { // Do with adaptation.
        do_adaptive_burnin(
            chains, 
            ntemperatures, 
            data, 
            burnin, 
            gen, 
            temp_window_size, 
            temp_window_size_growth_factor, 
            target_temp_swap_accept_ratio,
            temp_ladder_learning_rate,
            adaptation_swap_rates,
            verbose
        );
    }
    
    if (verbose != 0) {
        Rcpp::Rcout << "Burn-in complete. Starting sampling phase..." << std::endl;
    }

    // Sampling loop.
    std::vector<int> nswap_accepts(ntemperatures-1, 0);
    std::vector<int> nswap_proposals(ntemperatures-1, 0);
    std::vector<double> swap_probabilities(ntemperatures-1, 0);
    for (unsigned int iter = 0; iter < iterations; iter++) {
        do_step(chains, ntemperatures, data, nswap_accepts, nswap_proposals, 
            &swap_probabilities, gen);

        // Store samples.
        if (iter % thin == 0) {
            for (auto& chain : chains) {
                chain.storeSample();
            }
        }

        // Print progress updates.
        if (verbose > 0 && iter > 0 && (iter % verbose == 0 || iter == iterations - 1)) {

            // Compute mean swap ratio.
            double mean_swap_ratio = 0;
            for (int i = 0; i < ntemperatures - 1; i++) {
                mean_swap_ratio += swap_probabilities[i] 
                    / (nswap_proposals[i] == 0 ? 1 : nswap_proposals[i]);
            }
            mean_swap_ratio /= (ntemperatures - 1);

            // Print progress.
            Rcpp::Rcout << "Iteration " << (iter+1) << "/" << iterations 
                        << ", Mean swap acceptance ratio: " << mean_swap_ratio
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

    // Store inverse temperatures.
    arma::vec inv_temps(ntemperatures);
    arma::vec adapted_temps(ntemperatures);
    for (int i = 0; i < ntemperatures; i++) {
        inv_temps(i) = chains[i].inv_temperature;
        adapted_temps(i) = 1.0 / chains[i].inv_temperature;
    }

    // Return results.
    return Rcpp::List::create(
        _["storebeta"] = chains[0].store_beta,
        _["storegamma"] = storegamma_list,
        _["nswap_accepts"] = nswap_accepts,
        _["nswap_proposals"] = nswap_proposals,
        _["adapted_inv_temps"] = inv_temps,
        _["adapted_temps"] = adapted_temps,
        _["adaptation_swap_rates"] = adaptation_swap_rates
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
 * @param nswap_accepts A reference to a vector of integers that counts the number of accepted 
 * swaps for each pair of adjacent temperatures.
 * @param nswap_proposals A reference to a vector of integers that counts the number of proposed 
 * swaps for each pair of adjacent temperatures.
 * @param swap_probabilities A pointer to a vector of doubles where the swap probabilities for each
 * pair of adjacent temperatures will be accumulated. If this pointer is null, the swap probabilities 
 * will not be stored or accumulated.
 * @param rng The GSL random number generator to use for sampling.
 */
void do_step(
    std::vector<TemperatureChain>& chains, 
    int ntemperatures,
    const Data& data,
    std::vector<int>& nswap_accepts,
    std::vector<int>& nswap_proposals,
    std::vector<double>* swap_probabilities,
    gsl_rng* rng
) {
    // Do within-temperature MH updates.
    for (auto& chain : chains) {
        chain.do_step(data, rng);
    }

    // Select pair of adjacent temperatures to swap.
    int swap_index = gsl_rng_uniform_int(rng, ntemperatures - 1);
    TemperatureChain& chain1 = chains[swap_index];
    TemperatureChain& chain2 = chains[swap_index + 1];
    nswap_proposals[swap_index]++;
    double swap_probability = chain1.compute_swap_probability(chain2, data, rng);
    if (gsl_ran_flat(rng, 0.0, 1.0) <= swap_probability) {
        // Todo: check if we only want to do partial chain swapping.
        chain1.swap_gammas(chain2);
        chain1.swap_beta(chain2);
        nswap_accepts[swap_index]++;
    }
    // Store swap probability for temperature ladder adaptation.
    if (swap_probabilities != nullptr) {
        (*swap_probabilities)[swap_index] += swap_probability;
    }
}

/**
 * Perform the burn-in phase of the sampler, without adapting the temperature ladder.
 */
void do_burnin(
    std::vector<TemperatureChain>& chains, 
    int ntemperatures,
    const Data& data,
    int burnin, 
    gsl_rng* rng
) {
    std::vector<int> nswap_accepts(ntemperatures-1, 0);
    std::vector<int> nswap_proposals(ntemperatures-1, 0);
    for (unsigned int iter = 0; iter < burnin; iter++) {
        do_step(chains, ntemperatures, data, nswap_accepts, nswap_proposals, nullptr, rng);
    }
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
 */
void do_adaptive_burnin(
    std::vector<TemperatureChain>& chains, 
    int ntemperatures,
    const Data& data,
    int burnin, 
    gsl_rng* rng,
    int window_size,
    double window_growth_factor,
    double target_swap_ratio,
    double ladder_adjust_learning_rate,
    arma::vec& swap_rates,
    int verbose
) {
    std::vector<double> swap_probabilities(ntemperatures, 0);
    std::vector<int> nswap_accepts(ntemperatures-1, 0);
    std::vector<int> nswap_proposals(ntemperatures-1, 0);
    int window_step = 0;
    double min_gap = (1.0 - chains[ntemperatures-1].inv_temperature) / ((ntemperatures - 1) * 0.1);
    
    // Do burnin iterations.
    for (unsigned int iter = 0; iter < burnin; iter++) {
        do_step(chains, ntemperatures, data, nswap_accepts, nswap_proposals, 
                &swap_probabilities, rng);

        window_step++;
        if (window_step == window_size) {
            
            adjust_ladder(
                chains, 
                ntemperatures, 
                swap_probabilities, 
                swap_rates, 
                nswap_proposals, 
                target_swap_ratio, 
                ladder_adjust_learning_rate, 
                min_gap
            );

            swap_probabilities.assign(ntemperatures, 0);
            nswap_proposals.assign(ntemperatures-1, 0);
            nswap_accepts.assign(ntemperatures-1, 0);
            window_step = 0;
            window_size *= window_growth_factor;
        }

        // Print progress.
        if (verbose > 0 && iter > 0 && (iter % verbose == 0 || iter == burnin - 1)) {
            double mean_swap_ratio = 0;
            for (int i = 0; i < ntemperatures - 1; i++) {
                mean_swap_ratio += swap_probabilities[i] 
                    / (nswap_proposals[i] == 0 ? 1 : nswap_proposals[i]);
            }
            mean_swap_ratio /= (ntemperatures - 1);

            Rcpp::Rcout << "Burn-in iteration " << (iter+1) << "/" << burnin 
                        << ", Mean swap acceptance ratio: " << mean_swap_ratio
                        << std::endl;
         }
    }
}

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
void adjust_ladder(
    std::vector<TemperatureChain>& chains, 
    int ntemperatures,
    std::vector<double>& swap_probabilities,
    arma::vec& swap_rates,
    std::vector<int>& nswap_proposals,
    double target_swap_ratio,
    double ladder_adjust_learning_rate,
    double min_gap
) {
    std::vector<double> ladder_gaps(ntemperatures - 1);

    // Compute old gaps.
    for(int i = 0; i < ntemperatures - 1; i++) {
        ladder_gaps[i] = chains[i].inv_temperature - chains[i+1].inv_temperature;
    }

    // Compute mean swap rate for each pair.
    for (int i = 0; i < ntemperatures - 1; i++) {
        swap_rates[i] = swap_probabilities[i] / (nswap_proposals[i] == 0 ? 1 : nswap_proposals[i]);
    }

    // Update ladder gaps.
    for (int i = 0; i < ntemperatures - 1; i++) {
        ladder_gaps[i] *= std::exp(ladder_adjust_learning_rate * 
            (swap_rates[i] - target_swap_ratio));

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
    double normalization_factor = (1.0 - chains[ntemperatures-1].inv_temperature) / total_span;
    for (int i = 0; i < ntemperatures - 1; i++) {
        ladder_gaps[i] *= normalization_factor;
    }

    // Reconstruct ladder.
    chains[0].inv_temperature = 1.0;
    for (int i = 1; i < ntemperatures-1; i++) {
        chains[i].inv_temperature = chains[i-1].inv_temperature - ladder_gaps[i-1];
    }
}


