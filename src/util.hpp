#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace arma;
using namespace Rcpp;

namespace util {

mat rcpp2arma(Rcpp::NumericMatrix &m) {
  return mat(m.begin(), m.nrow(), m.ncol(), false);
}

vec rcpp2arma(Rcpp::NumericVector &v) {
  return vec(v.begin(), v.size(), false);
}

inline mat rmvnorm(int n, const vec &mu, const mat &Sigma) {
  int ncols = Sigma.n_cols;
  mat Y = randn(n, ncols);
  return trans(repmat(mu, 1, n)) + Y * chol(Sigma);
}

// inline mat rnmvnorm_fast(int n, const vec &mu, const mat &Sigma) {
//   return trans(repmat(mu, 1, n)) + randn(n, Sigma.n_cols) * chol(Sigma);
// }

inline colvec rmvnorm_fast(const gsl_rng * gen, const colvec mu, const mat Sigma) {
  std::vector<double> sample(Sigma.n_cols);
  for (unsigned int i; i<Sigma.n_cols; ++i) {
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

// linear regression with Gaussian errors beta draw 
// (multivariate Normal prior)
// regression model is y = X * beta + epsilon,  epsilon ~ N(0,sigma2)
// XpX is X'X
// XpY is X'y
// b0 is the prior mean of beta
// B0 is the prior precision (the inverse variance) of beta
colvec
NormNormregress_beta_draw (const gsl_rng * gen,
                           const mat &XpX,
                           const mat &XpY,
                           const colvec &b0,
                           const mat &B0,
                           const double &sigma2){

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

// For use by ordered and hybrid probit, see Albert & Chib 

// Transform gamma to alpha
mat gamma2alpha(const mat& gamma){
  const int m = gamma.n_rows - 2;
  mat alpha(m, 1);
  alpha[0] = std::log(gamma[1]);
  for (int j=1; j< m ; ++j){
    alpha[j] = std::log(gamma[j+1] - gamma[j]);
  }
  return alpha;
}

// Transform alpha to gamma 
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