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
 * Struct to hold the data for the model. This includes the feature matrices and response vectors
 * for each target, as well as the combined feature matrix and precomputed X'X matrix for all 
 * targets.
 */
struct Data {

    /** The feature matrix for each target. */
    std::vector<arma::mat> X; 

    /** The integer category labels for each target */
    std::vector<arma::colvec> Y; 

    /** The combined (stacked) feature matrix for all targets. */
    arma::mat Xall;

    /** The precomputed X'X matrix for the combined feature matrix for all targets. */
    arma::mat XpX;

    /** The number of observations for each target. */
    std::vector<int> target_nobs;

    /** The total number of observations. */
    int nobs;

    int ntargets;

    int npredictors;
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


/**
 * Propose a new gamma via the Metropolis-Hastings method.
 * 
 * @param gamma The current gamma thresholds for the target, which will be updated in place if the
 * proposed gammas are accepted.
 * @param beta The current regression coefficients for the target.
 * @param X The feature matrix for the target.
 * @param Y The response vector for the target.
 * @param ncategories The number of categories for the target.
 * @param sigma The standard deviation of the truncated Gaussian distribution used for proposing new
 * gamma values. This controls the tuning of the proposal distribution.
 * @param inv_temperature The inverse temperature (1/T) for the current chain. This is used to scale
 * the log likelihood ratio in the acceptance probability for the proposed gammas.
 * @param acceptance_probability The resulting acceptance probability for the proposed gammas.
 * @param rng The GSL random number generator to use for sampling.
 * 
 * @return A boolean indicating whether the proposed gammas were accepted (true) or rejected (false).
 */
bool mh_update_gamma(
    colvec& gamma,
    const colvec& beta,
    const mat& X,
    const colvec& Y,
    int ncategories,
    double sigma,
    double inv_temperature,
    double& acceptance_probability,
    gsl_rng* rng
);

/**
 * Perform the Gibbs update for the regression coefficients. This involves first updating the
 * latent variables y* based on the current thresholds and regression coefficients, and then
 * drawing new regression coefficients from their full conditional distribution given the 
 * updated latent variables and thresholds.
 * 
 * @param beta The current regression coefficients, which will be updated in place with the new 
 * sampled values.
 * @param gamma The current threshold parameters for each target.
 * @param data The data object.
 * @param beta_mean_prior The prior mean for the regression coefficients.
 * @param beta_prec_prior The prior precision matrix for the regression coefficients.
 * @param rng The GSL random number generator to use for sampling.
 */
void gibbs_update_beta(
    colvec& beta,
    const std::vector<colvec>& gamma,
    const Data& data,
    const arma::colvec& beta_mean_prior,
    const arma::mat& beta_prec_prior,
    gsl_rng* rng
);


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
);

#endif // __SAMPLING_HPP