/**
 * @file mspm_util.hpp
 * @author Johan Falkenjack
 * @author Daniel Tufvesson
 * @version 1.0
 * @brief Utility functions for mspm.
 */

#ifndef __MSPM_UTIL_HPP
#define __MSPM_UTIL_HPP

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace arma;
using namespace Rcpp;

namespace mspm_util {

inline mat rcpp2arma(Rcpp::NumericMatrix &m) {
    return mat(m.begin(), m.nrow(), m.ncol(), false);
}

inline vec rcpp2arma(Rcpp::NumericVector &v) {
    return vec(v.begin(), v.size(), false);
}

inline mat rmvnorm(int n, const vec &mu, const mat &Sigma) {
    int ncols = Sigma.n_cols;
    mat Y = randn(n, ncols);
    return trans(repmat(mu, 1, n)) + Y * chol(Sigma);
}

inline colvec rmvnorm_fast(const gsl_rng * gen, const colvec mu, const mat Sigma) {
    std::vector<double> sample(Sigma.n_cols);
    for (unsigned int i = 0; i < Sigma.n_cols; ++i) {
        sample[i] = gsl_ran_ugaussian(gen);
    }
    colvec vec_sample(sample);
    return mu + trans(chol(Sigma)) * vec_sample;
}

inline double rscinvchisq_safe(const gsl_rng * gen, const double df, const double scale) {
    if (df <= 0) {
        Rcerr << "df must be greater than zero\n";
        return 0;
    }
    if (scale <= 0) {
        Rcerr << "scale must be greater than zero\n";
        return 0;
    } 
    return (df * scale)/gsl_ran_chisq(gen, df);
}

inline double rscinvchisq(const gsl_rng * gen, const double df, const double scale) {
    return (df * scale)/gsl_ran_chisq(gen, df);
}

inline double dscinvchisq(double x, double df, double scale, bool logprob) {
    if (df <= 0) {
        Rcerr << "df must be greater than zero\n";
        return 0;
    }
    if (scale <= 0) {
        Rcerr << "scale must be greater than zero\n";
        return 0;
    }
    if (x < 0) {
        return 0;
    }
    double nu = df/2;
    if (logprob) {
        return nu * log(nu) - lgamma(nu) + 
            nu * log(scale) - (nu + 1) * log(x) - (nu * scale/x);
    } else {
        return ((pow(nu, nu))/tgamma(nu)) * pow(scale, nu) * 
            pow(x, (-(nu + 1))) * exp(-nu * scale/x);
    }
}

inline double dgamma(double x, double df, double scale) {
    return gsl_ran_gamma_pdf(x, df, 1/scale);
}

inline double rgamma(const gsl_rng * gen, const double df, const double scale) {
    return gsl_ran_gamma(gen, df, 1/scale);
}

inline double log_dgamma(double x, double k, double theta) {
    return -lgamma(k) - k * (-log(theta)) + (k-1)*log(x) - x*theta;
}

inline double log_dnorm2(double x, double mean, double sd) {
    return -log(2*datum::pi*pow(sd, 2))/2 - (pow(x-mean, 2)/(2*pow(sd, 2)));
}

inline double log_dnorm(double x) {
    return log_dnorm2(x, 0, 1);
}

inline double log_pnorm(double x) {
    return log1p(-exp(1.4*x)) - log(-x) - pow(x, 2)/2 - 1.04557;
}

// WARNING! UNSAFE
inline double log_subtract(double x, double y) {
    double a = -(y-x);
    if (a < 0.693) {
        return x + log(-expm1(-a));
    } else {
        return x + log1p(-exp(-a));
    }
}

// inline double dscinvchisq(double x, double df, double scale) {
//   if (df <= 0) {
//     Rcerr << "df must be greater than zero\n";
//     return 0;
//   }
//   if (scale <= 0) {
//     Rcerr << "scale must be greater than zero\n";
//     return 0;
//   }
//   if (x < 0) {
//     return 0;
//   }
//   double nu = df/2;
//   return ((pow(nu, nu))/gamma(nu)) * pow(scale, nu) * 
//       pow(x, (-(nu + 1))) * exp(-nu * scale/x);
// }

inline mat crossprod(const mat &m1, const mat &m2) {
    return(trans(m1) * m2);
}

/**
 * Linear regression with Gaussian errors beta draw (multivariate Normal prior). The regression 
 * model is y = X * beta + epsilon,  epsilon ~ N(0,sigma2). 
 * 
 * @param gen a GSL random number generator.
 * @param XpX the cross-product matrix X'X.
 * @param XpY the matrix X'y.
 * @param b0 the prior mean of beta.
 * @param B0 the prior precision (the inverse variance) of beta.
 * @param sigma2 the variance of the errors.
 * @return a draw of beta from the posterior distribution.
 */
colvec NormNormregress_beta_draw (
    const gsl_rng * gen,
    const mat &XpX,
    const mat &XpY,
    const colvec &b0,
    const mat &B0,
    const double &sigma2
);

// For use by ordered and hybrid probit, see Albert & Chib 

// Transform gamma to alpha
mat gamma2alpha(const mat& gamma);

// Transform alpha to gamma 
mat alpha2gamma(const mat& alpha);

}
#endif