/**
 * @file sampling.hpp
 * @author Daniel Tufvesson
 * @version 1.0
 * @brief Header file for sampling functions for the Multi-Scale Probit model.
 */

#ifndef __SAMPLING_HPP
#define __SAMPLING_HPP

#include <RcppArmadillo.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "mspm_util.hpp"
#include "rtnorm.hpp"
#include <algorithm>
#include <cmath>

/**
 * Compute the acceptance rates for the proposed gammas for each target based on the observed 
 * acceptance probabilities for the proposed gammas.
 * 
 * @param acceptance_probabilities The vector containing the sum of the acceptance probabilities for 
 * the proposed gammas for each target.
 * @param nsamples The number of samples (iterations) over which the acceptance probabilities were
 * accumulated.
 * @param acceptance_rates The vector to store the computed acceptance rates for each target. This will
 * be updated in place with the computed acceptance rates for each target.
 */
void compute_acceptance_rate(
    const arma::vec& acceptance_probabilities,
    int nsamples,
    arma::vec& acceptance_rates
);

/**
 * Struct to hold the data for the model. This includes the feature matrices and response vectors
 * for each target, as well as the combined feature matrix and precomputed X'X matrix for all 
 * targets.
 */
struct Data {
public:
    /** The feature matrix for each target. */
    std::vector<arma::mat> X; 

    /** The integer category labels for each target */
    std::vector<arma::colvec> Y; 

    /** The combined (stacked) feature matrix for all targets. */
    arma::mat Xall;

    /** The precomputed X'X matrix for the combined feature matrix for all targets. */
    arma::mat XpX;

    /** The number of observations for each target. */
    arma::ivec target_nobs;

    /** The total number of observations. */
    int nobs;

    int ntargets;

    int npredictors;
};

class MspmChain {
public:
    /** The inverse temperature 1/T, also known as beta in parallel tempering literature. */
    double inv_temperature;

    const arma::colvec beta_start;

    const std::vector<arma::colvec> gammas_start;

    /** The current beta state. */
    arma::colvec beta;

    /** The current gammas state. */
    std::vector<arma::colvec> gammas;

    /** The current latent value state. */
    arma::colvec ystar;

    /** The proposal variance for the gamma proposal distribution. */
    arma::vec proposal_variance;

    /** The prior mean for the beta. */
    arma::colvec beta_mean_prior;

    /** The prior precision for the beta. */
    arma::mat beta_prec_prior;

    /** The number of categories for each target. */
    arma::ivec ncategories;

    /** The number of steps the chain has progressed so far. */
    int nsteps = 0;

    /** The MH proposal probabilities will be added to this for each target. */
    arma::vec cumulative_acceptance_probabilities;

    /** Storage of gamma proposals. Mostly to avoid reallocating the same colvecs each iteration. */
    std::vector<arma::colvec> gamma_proposals;

    /** The number of predictor covariates. */
    const int npredictors;

    /** The number of target gamma groups. */
    const int ntargets;

    /** The cached log likelihood for each target for the current state. These are unscaled
     * by temperature. */
    arma::vec log_likelihood;

    /** How many times the compute_log_likelihood function has been called. */
    int nlikelihood_calls = 0;

    /** A label for tracking the swapping paths. */
    int replica_id;

    MspmChain(
        double inv_temperature,
        const arma::colvec beta_start,
        const std::vector<arma::colvec> gammas_start,
        const arma::vec& proposal_variance,
        const arma::colvec& beta_mean_prior,
        const arma::mat& beta_prec_prior,
        const arma::ivec& ncategories,
        const int total_nobs,
        int replica_id
    ) : inv_temperature(inv_temperature), beta_start(beta_start), gammas_start(gammas_start),
        beta(beta_start), gammas(gammas_start), 
        proposal_variance(proposal_variance), beta_mean_prior(beta_mean_prior), 
        beta_prec_prior(beta_prec_prior), ncategories(ncategories), npredictors(beta_start.n_elem), 
        ntargets(gammas_start.size()), replica_id(replica_id) {
        
        ystar = arma::colvec(total_nobs, arma::fill::zeros);
        cumulative_acceptance_probabilities = arma::vec(gammas.size(), arma::fill::zeros);
        log_likelihood = arma::vec(gammas.size(), arma::fill::zeros);
        
        // Initialize gamma proposal storage.
        gamma_proposals = std::vector<arma::colvec>(ntargets);
        for (int target = 0; target < ntargets; target++) {
            // We include -inf and +inf edge gammas too.
            gamma_proposals[target] = arma::colvec(ncategories(target)+1, arma::fill::zeros);
        }
    }

    /**
     * Reset the chain to starting state. 
     */
    void reset() {
        beta = beta_start;
        for (int target = 0; target < ntargets; target++) {
            gammas[target] = gammas_start[target];
        }
        ystar.zeros();
        log_likelihood.zeros();
        reset_statistics();
    }

    /**
     * Reset only the statistics of the chain, such as the cumulative acceptance probabilities and 
     * the number of likelihood calls, but keep the current state of the chain.
     */
    void reset_statistics() {
        cumulative_acceptance_probabilities.zeros();
        nlikelihood_calls = 0;
    }

    /**
     * Perform one step of the within-temperature Metropolis-Hastings update for the thresholds 
     * and Gibbs update for the regression coefficients.
     * 
     * @param data The data object.
     * @param rng The GSL random number generator to use for sampling.
     */
    void simulate_step(const Data& data, gsl_rng* rng) {
        if (nsteps == 0) {
            // Compute initial log likelihoods.
            refresh_log_likelihood(data);
        }
        step_gamma(data, rng);
        step_beta(data, rng);
        nsteps++;
    }

    void refresh_log_likelihood(const Data& data) {
        for (int target = 0; target < ntargets; target++) {
            log_likelihood(target) = compute_log_likelihood(gammas[target], data, target);
        }
    }

private:

    /**
     * Perform the Metropolis-Hastings updates for the thresholds for each target. This includes
     * proposing new threshold values and doing the acceptance step.
     * 
     * @param data The data object.
     * @param rng The GSL random number generator to use for sampling.
     */
    void step_gamma(const Data& data, gsl_rng* rng) {
        for (int t = 0; t < ntargets; t++) {
            int target = (nsteps + t) % ntargets;
            int ncats = ncategories(target);

            // Propose new gamma values for the target.
            propose_gamma(gammas[target], gamma_proposals[target], ncats, 
                proposal_variance[target], rng);

            double log_likelihood_proposed = compute_log_likelihood(gamma_proposals[target], 
                data, target);

            double log_likelihood_ratio = log_likelihood_proposed - log_likelihood(target);
            
            double log_proposal_ratio = compute_gamma_log_proposal_ratio(gammas[target], 
                gamma_proposals[target], ncats, proposal_variance[target]);

            double log_accept_ratio = inv_temperature * log_likelihood_ratio + log_proposal_ratio;

            // Do MH acceptance step.
            double acceptance_probability = std::min(1.0, std::exp(log_accept_ratio));
            if (gsl_ran_flat(rng, 0.0, 1.0) <= acceptance_probability) {
                gammas[target] = gamma_proposals[target];
                log_likelihood(target) = log_likelihood_proposed;
            }
            cumulative_acceptance_probabilities[target] += acceptance_probability;
        }
    }

    /**
     * Propose new threshold values for a given target using a truncated normal distribution. 
     * 
     * @param gamma The current gamma thresholds for the target.
     * @param gamma_prop The colvec to store the proposed gamma thresholds for the target. This 
     * will be updated in place with the proposed gamma values.
     * @param ncategories The number of categories for the target.
     * @param sigma The standard deviation of the truncated Gaussian distribution used for proposing new
     * gamma values. This controls the tuning of the proposal distribution.
     * @param rng A pointer to a GSL random number generator object, which is used to draw random
     * samples from the truncated normal distribution.
     * 
     * @return A colvec containing the proposed new gamma thresholds for the target. The length of the
     * returned colvec will be equal to ncategories - 1, since there are ncategories - 1 thresholds 
     * for a target with ncategories.
     */
    void propose_gamma(
        const colvec& gamma,
        arma::vec& gamma_prop,
        int ncategories, 
        double sigma,
        gsl_rng* rng
    ) {
        gamma_prop.head(1) = -INFINITY;
        gamma_prop.tail(1) = INFINITY;
        // Note: probabilities are 0 for the edge thresholds.

        if (ncategories == 2) {
            // Draw new split point
            gamma_prop(1) = rtnorm(rng, -INFINITY, INFINITY, gamma(1), sigma).first;
        } 
        else {
            for (unsigned int i = 1; i < ncategories; i++) {
                if (i == 1) { // If first gamma
                    gamma_prop(i) = rtnorm(rng, gamma(i-1), gamma(i+1), gamma(i), sigma).first;
                } 
                else { // If any other gamma
                    gamma_prop(i) = rtnorm(rng, gamma_prop(i-1), gamma(i+1), gamma(i), sigma).first;
                }
            }
        }
    }

    /**
     * Compute the log likelihood for the proposed gamma target.
     */
    double compute_log_likelihood(
        const arma::colvec& gamma_p,
        const Data& data,
        int target
    ) {
        int ncats = ncategories[target];
        double log_likelihood = 0.0;
        auto cdf = gsl_cdf_ugaussian_P;
        arma::colvec y_star = data.X[target] * beta;
        for (int i = 0; i < data.target_nobs[target]; i++) {
            int y_val = data.Y[target][i];

            // Handle last category.
            if (y_val == ncats){
                log_likelihood += log(1.0 - cdf(gamma_p(y_val-1) - y_star[i]));
            }
            // Handle first category.
            else if (y_val == 1) {
                log_likelihood += log(cdf(gamma_p(y_val) - y_star[i]));
            }
            // Handle categories inbetween.
            else {
                log_likelihood += log(cdf(gamma_p(y_val) - y_star[i]) 
                    - cdf(gamma_p(y_val-1) - y_star[i]));
            }
        }
        nlikelihood_calls++;
        return log_likelihood;
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
     * Compute the log probability ratio betweent the proposal distribution for the proposed gammas,
     * i.e. the log of q(ɣ|ɣ') / q(ɣ'|ɣ), where ɣ is the old gamma, ɣ' is the proposed gamma and 
     * q(.) is the truncated Gaussian proposal distribution.
     * 
     * Note that this is not the full proposal ratio, but only the part of the proposal ratio that
     * depends on the proposed gammas. The full proposal ratio also includes the probabilities of the
     * proposed gammas under the truncated normal distribution.
     * 
     * @param gamma The current gamma thresholds for the target.
     * @param gamma_p The proposed new gamma thresholds for the target.
     * @param ncategories The number of categories for the target.
     * @param sigma The standard deviation of the truncated Gaussian distribution used for proposing 
     * new gamma values.
     * 
     * @return The log of the proposal ratio for the proposed gammas.
     */
    double compute_gamma_log_proposal_ratio(
        const colvec& gamma,
        const colvec& gamma_p,
        int ncategories,
        double sigma
    ) {
        double log_proposal_ratio = 0;
        for (int i = 1; i < ncategories; i++) {
            // Normalization term Z(gamma) 
            double z_old = compute_log_trunc_gauss_norm_constant(
                gamma(i), 
                gamma(i-1), 
                gamma(i+1), 
                sigma
            );

            if (i == 1) { // If first gamma
                // Normalization term Z(gamma_p)
                double z_prop = compute_log_trunc_gauss_norm_constant(
                    gamma_p(i), 
                    gamma(i-1), 
                    gamma(i+1), 
                    sigma
                );
                
                log_proposal_ratio += z_old - z_prop;
            } 
            else { // If any other gamma
                // Normalization term Z(gamma_p)
                double z_prop = compute_log_trunc_gauss_norm_constant(
                    gamma_p(i), 
                    gamma_p(i-1), 
                    gamma(i+1), 
                    sigma
                );
                log_proposal_ratio += z_old - z_prop;
            }
        }
        return log_proposal_ratio;
    }

    void step_beta(const Data& data, gsl_rng* rng) {
        // Draw latent y*.
        arma::colvec ystar = arma::colvec(data.nobs, arma::fill::zeros);
        int offset = 0;
        for (unsigned int target = 0; target < data.ntargets; target++) {
            if (target > 0) {
                int nobs = data.Y[target-1].n_elem;
                offset += nobs;
            }
            const colvec current_ystar = data.X[target] * beta;
            for (unsigned int i = 0; i < data.X[target].n_rows; i++) {
                ystar(offset + i) = rtnorm(
                    rng,
                    gammas[target](data.Y[target](i) - 1),
                    gammas[target](data.Y[target](i)),
                    current_ystar[i], 
                    1.0
                ).first;
            }
        }

        // Draw new beta.
        arma::mat XpZ = arma::trans(data.Xall) * ystar;
        beta = mspm_util::NormNormregress_beta_draw(rng, data.XpX, XpZ, beta_mean_prior, 
            beta_prec_prior, 1.0);
    }

public:

    /**
     * Adjust the proposal variance for the gamma parameters based on the observed acceptance 
     * probabilities for the proposed gamma values. This is done by comparing the observed acceptance 
     * rates to a target acceptance rate, and adjusting the standard deviation of the truncated normal 
     * proposal distribution for the gammas accordingly.
     * 
     * @param acceptance_rates A vector of observed acceptance rates for the proposed gamma values 
     * for each target.
     * @param target_acceptance_rate The target acceptance rate for the proposed gammas.
     * @param learning_rate The learning rate for adjusting the proposal variance.
     */
    void adjust_proposal_variance(
        const arma::vec& acceptance_rates,
        const double target_acceptance_rate,
        const double learning_rate
    ) {
        for (int i = 0; i < ntargets; i++) {
            double acceptance_rate = acceptance_rates(i);
            double log_sigma = std::log(proposal_variance(i));
            log_sigma += learning_rate * (acceptance_rate - target_acceptance_rate);
            proposal_variance(i) = std::exp(log_sigma);
        }
    }



    // Parallel tempering logic. -----------------------------------------------------------------

    void swap(MspmChain& other_chain, const Data& data, bool complete_swap) {
        // Swap gammas.
        for (int target = 0; target < gammas.size(); target++) {
            for (int j = 0; j < gammas[target].n_elem; j++) {
                double temp = gammas[target][j];
                gammas[target][j] = other_chain.gammas[target][j];
                other_chain.gammas[target][j] = temp;
            }
        }

        // Swap betas.
        if (complete_swap) {
            for (int j = 0; j < npredictors; j++) {
                double temp = beta[j];
                beta[j] = other_chain.beta[j];
                other_chain.beta[j] = temp;
            }
        }

        // Swap replica labels.
        int temp_label = replica_id;
        replica_id = other_chain.replica_id;
        other_chain.replica_id = temp_label;

        // New state means new likelihoods, so recompute them.
        refresh_log_likelihood(data);
    }

    /**
     * Determine whether to swap the state of this chain with another chain based on the computed
     * acceptance ratio for the swap. The acceptance ratio is computed based on the likelihood of 
     * the data given the current parameter values in each chain, and the difference in inverse
     * temperatures between the two chains.
     * 
     * @param other_chain The other MspmChain instance with which to potentially swap states.
     * @param data The data object containing the feature matrices and response vectors for each 
     * target.
     * @param rng The GSL random number generator to use for sampling the uniform random variable 
     * for the acceptance step.
     */
    double compute_swap_probability(
        const MspmChain& other_chain, 
        const Data& data,
        gsl_rng* rng
    ) {
        double log_swap_accept_ratio = arma::sum(other_chain.log_likelihood) 
            - arma::sum(log_likelihood);

        // Scale it by inv_temperature delta.
        log_swap_accept_ratio *= (inv_temperature - other_chain.inv_temperature);
        
        // Do the swap with the computed acceptance ratio.
        return std::min(1.0, exp(log_swap_accept_ratio));
    }
};

class SampleStorage {
public:

    /** Storage matrix for sampled betas. */
    arma::mat store_beta;

    /** Vector of storage matrixes for each target gamma group. */
    std::vector<arma::mat> store_gammas;

    /** The number of samples stored so far. */
    int nstored = 0;

    /** The number of target gamma groups. */
    const int ntargets;

    /** The number of samples to store. */
    const int nstore_size;

    SampleStorage(
        int nstore_size,
        int npredictors,
        arma::ivec ncategories
    ) : ntargets(ncategories.n_elem), nstore_size(nstore_size) {

        // Initialize storage for beta and gammas.
        store_beta = arma::mat(nstore_size, npredictors, arma::fill::zeros);
        store_gammas = std::vector<arma::mat>(ntargets);
        for (int target = 0; target < ntargets; ++target) {
            store_gammas[target] = arma::mat(nstore_size, ncategories(target)-1, arma::fill::zeros);
        }
    }

    /**
     * Store the current state of the chain.
     */
    void store_sample(MspmChain& chain) {
        // Store beta.
        for (unsigned int j = 0; j < chain.beta.n_rows; j++) {
            store_beta(nstored, j) = chain.beta[j];
        }
        // Store gamma.
        for (unsigned int target = 0; target < chain.gammas.size(); target++) {
            // Note: the first and last gammas should not be saved.
            for (unsigned int j = 1; j < chain.gammas[target].n_rows - 1; j++) {
                store_gammas[target](nstored, j-1) = chain.gammas[target](j);
            }
        }
        nstored++;
    }

    Rcpp::List gamma_to_r_list() {
        Rcpp::List gammas = Rcpp::List::create();
        for (unsigned int target = 0; target < ntargets; ++target) {
            gammas.push_back(store_gammas[target]);
        }
        return gammas;
    }
};

/**
 * A class for tuning the proposal variance of the MSPM sampler.
 */
class ProposalTuner {
public:

    /** The target acceptance-rejection rate for the proposed gammas. */
    const double target_acceptance_rate;

    /** The epsilon threshold for early stopping. */
    const double target_epsilon;

    /** The window size for computing acceptance rates and adjusting the proposal variance 
     * during tuning */
    const int window_size;

    /** The total number of iterations the tuner has run. */
    int step = 0;

    /** The number of iterations the tuner has taken within the current window. */
    int window_step = 0;

    /** The acceptance rates for the current window. Is reset at the beginning of each window. */
    arma::vec acceptance_rates;

    ProposalTuner(
        double target_acceptance_rate,
        double target_epsilon,
        int window_size,
        int ntargets
    ) : target_acceptance_rate(target_acceptance_rate), target_epsilon(target_epsilon), 
        window_size(window_size) {

        acceptance_rates = arma::vec(ntargets, arma::fill::zeros);
    }

    /**
     * Run one step of the tuning process, which includes updating the acceptance probabilities for 
     * the proposed gammas, and adjusting the proposal variance if the end of the window is reached.
     * 
     * Note that this function does not simulate a step of the chain itself. Simulating a step of
     * the chain should be done before calling this function.
     * 
     * @param chain The MspmChain instance for which to perform the tuning step. 
     * @param data The data object containing the feature matrices and response vectors for each target.
     * @param rng The GSL random number generator to use for sampling.
     */
    void tune_step(MspmChain& chain, const Data& data, gsl_rng* rng) {

        window_step++;
        if (window_step == window_size){
            window_step = 0;

            // Do adjustment.
            double learning_rate = 1.0 / std::sqrt(step);
            compute_acceptance_rate(chain.cumulative_acceptance_probabilities, window_size, 
                acceptance_rates);
            chain.adjust_proposal_variance(acceptance_rates, target_acceptance_rate, learning_rate);
            chain.cumulative_acceptance_probabilities.zeros();
        }

        step++;
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
        for (int target = 0; target < acceptance_rates.n_elem; target++) {
            if (std::abs(acceptance_rates(target) - target_acceptance_rate) > target_epsilon) {
                return false;
            }
        }
        return true;
    }
};

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
);

/**
 * Unpack gamma thresholds from Rcpp::List and set boundary values.
 * 
 * @param gammaStart Rcpp::List of initial gamma vectors for each target.
 * @param ncat arma::ivec with number of categories for each target.
 * @return std::vector<colvec> with boundary values set.
 */
std::vector<colvec> unpack_gamma(const Rcpp::List& gammaStart, const arma::ivec& ncat);



#endif // __SAMPLING_HPP