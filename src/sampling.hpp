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
 * Compute the log likelihood ratio between the current gammas and the proposed gammas.
 * 
 * @param gamma The current gamma thresholds for the target.
 * @param gamma_p The proposed gamma thresholds for the target.
 * @param X The feature matrix for the target.
 * @param Y The response vector for the target.
 * @param beta The current regression coefficients.
 * @param ncategories The number of categories for the target.
 * 
 * @return The log likelihood ratio for the proposed gamma thresholds compared to the current 
 * gamma thresholds.
 */
double compute_log_likelihood_ratio(
    const colvec& gamma,
    const colvec& gamma_p,
    const mat& X,
    const colvec& Y,
    const colvec& beta,
    int ncategories
);

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
);

/**
 * Propose new threshold values for a given target using a truncated normal distribution. 
 * 
 * @param gamma The current gamma thresholds for the target.
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
arma::colvec propose_gamma(
    const colvec& gamma,
    int ncategories, 
    double sigma,
    gsl_rng* rng
);

/**
 * Propose a new gamma via the Metropolis-Hastings method.
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
    gsl_rng* rng
);


#endif // __SAMPLING_HPP