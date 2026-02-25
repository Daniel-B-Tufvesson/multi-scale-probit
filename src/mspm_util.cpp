/**
 * @file mspm_util.cpp
 * @author Johan Falkenjack
 * @author Daniel Tufvesson
 * @version 1.0
 * @brief Utility functions for mspm.
 */

#include "mspm_util.hpp"
namespace mspm_util {


colvec NormNormregress_beta_draw(
    const gsl_rng * gen,
    const mat &XpX,
    const mat &XpY,
    const colvec &b0,
    const mat &B0,
    const double &sigma2
) {
    // this function gets the cross-product matrix X'X and the matrix X'Y
    // to minimize the amount of computation within the function
    const int k = XpX.n_cols;
    const double sig2_inv = 1.0 / sigma2;
    const mat sig_beta = inv_sympd(B0 + XpX * sig2_inv);
    const mat C = chol(sig_beta, "lower");
    const mat betahat = sig_beta * (B0*b0 + XpY*sig2_inv);
    std::vector<double> sample(k);
    for (unsigned int i; i<k; ++i) {
        sample[i] = gsl_ran_ugaussian(gen);
    }
    colvec vec_sample(sample);
    return(C * vec_sample + betahat);
}

mat gamma2alpha(const mat& gamma){
    const int m = gamma.n_rows - 2;
    mat alpha(m, 1);
    alpha[0] = std::log(gamma[1]);
    for (int j=1; j< m ; ++j){
        alpha[j] = std::log(gamma[j+1] - gamma[j]);
    }
    return alpha;
}

mat alpha2gamma(const mat& alpha){
    const int m = alpha.n_rows;
    mat gamma(m+2, 1);
    gamma[0] = -300;
    gamma[m+1] = 300;
    double gamma_sum;
    for (int j=1; j<m+1 ; ++j){
        gamma_sum = 0.0;
        for(int i=0;i<j; ++i){
            gamma_sum += exp(alpha[i]);
        }
        gamma[j] = gamma_sum;
    }
    return gamma;
}

}