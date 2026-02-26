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
};


#endif // __SAMPLING_HPP