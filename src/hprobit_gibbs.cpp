/**
 * @file hprobit_gibbs.cpp
 * @author Johan Falkenjack
 * @author Daniel Tufvesson
 * @version 1.0
 * @brief A Metropolis-Hastings-withing-Gibbs sampler for the Multi-Scale Probit model.
 */

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "rtnorm.hpp"
#include <utility> 
#include "util.hpp"

#define TESTBROKE 1
// #define TESTCOMP 1

using namespace Rcpp;
using namespace R;
using namespace arma;
using namespace util;


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List cpp_hprobit(const Rcpp::List Xlist, // Todo: pass by reference instead?
                       const Rcpp::List Ylist,
                       const arma::colvec meanPrior,
                       const arma::mat precPrior,
                       const int fixZero,
                       const arma::ivec ncat,
                       const Rcpp::List gammaStart,
                       const arma::colvec betaStart,
                       const arma::vec tune,
                       const int iterations,
                       const int burnin,
                       const int thin,
                       const int seed,
                       const int verbose) {
  
  //--- GSL random init ---
  // Used to take efficient samples from truncated normal.
  gsl_rng_env_setup();                          // Read variable environnement
  const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
  gsl_rng *gen = gsl_rng_alloc (type);          // Rand generator allocation
  if (seed == 0) {
    unsigned int seed = floor(randu()*100000);
  }
  // Rcout << "Setting gsl seed to " << seed << "\n";
  gsl_rng_set(gen, seed);
  
  // Define constants
  const int tot_iter = iterations+burnin;
  const int nstore = iterations/thin;
  const int ntargets = Xlist.size();
  double gammamax = std::numeric_limits<double>::max();
  double gammamin = -std::numeric_limits<double>::max();
  

  // Unpack Xlist, Ylist and gammaStart
  std::vector<mat> X(ntargets);
  std::vector<colvec> Y(ntargets);
  std::vector<colvec> gamma(ntargets);
  std::vector<colvec> new_gamma(ntargets);
  mat Xall;
  colvec Yall;
  std::vector<int> sizes(ntargets);
  for (unsigned int target = 0; target < ntargets; ++target) {
    X[target] = Rcpp::as<mat>(Xlist[target]);
    Y[target] = Rcpp::as<colvec>(Ylist[target]);
    gamma[target] = Rcpp::as<colvec>(gammaStart[target]);
    if (target == 0) {
      Xall = X[target];
      Yall = Y[target];
    } else {
      Xall = join_vert(Xall, X[target]);
      Yall = join_cols(Yall, Y[target]);
    }
    sizes[target] = Y[target].n_elem;
  }
  // Define constants
  const int k = Xall.n_cols;
  const int N = Xall.n_rows;
  const mat XpX = crossprod(Xall, Xall);
  
  // Update gammas with system double extremes
  for (unsigned int target = 0; target < ntargets; ++target) {
    gamma[target](0) = gammamin;
    gamma[target](ncat[target]) = gammamax;
  }
  
  
  // Storage matrices
  mat storebeta(nstore, Xall.n_cols, arma::fill::zeros);
  std::vector<mat> storegamma(ntargets);
  
  for (unsigned int target = 0; target < ntargets; ++target) {
    storegamma[target] = mat(nstore, ncat(target)-1, arma::fill::zeros);
  }
  
  
  // #ifdef TESTBROKE
  //   Rcout << "Setting starting points. \n";
  //   sleep(1);
  // #endif
  
  // Set starting points
  colvec beta = betaStart;
  
  // #ifdef TESTBROKE
  //   Rcout << "gamma[0]: "<< gamma[0] <<"\n";
  //   Rcout << "gamma[1]: "<< gamma[1] <<"\n";
  //   Rcout << "gamma[2]: "<< gamma[1] <<"\n";
  //   Rcout << "gamma[3]: "<< gamma[1] <<"\n";
  //   sleep(1);
  // #endif
  
  // Set Z vector starting point as OLS estimates

  colvec Z(N, arma::fill::zeros);// = X * beta;
  colvec Xbeta;
  
  // Bookkeeping
  unsigned int count = 0;
  ivec accepts(ntargets+1, arma::fill::zeros);
  int offset;

  // Gibbs loop
  for (unsigned int iter = 0; iter<tot_iter; ++iter) {
    // Rcout << "Starting iteration. \n";
    // sleep(1);

    // Step 1: Update gammas: (gamma | u, beta)
    // Cowles update of gamma
    offset = 0;
    // double loglikerat = 0.0;
    for (unsigned int t = 0; t < ntargets; ++t) {
      int target = (iter+t) % ntargets;
      colvec gamma_p(ncat(target)+1, arma::fill::zeros);
      gamma_p.head(1) = -INFINITY;
      gamma_p.tail(1) = INFINITY;
      if (ncat(target) == 2) {
          // Draw new split point
          gamma_p(1) = rtnorm(gen,
                              -INFINITY,
                              INFINITY,
                              gamma[target](1),
                              tune(target)
                              ).first;
        // }
      } else {
        for (int i=1; i<(ncat(target)); ++i){
            if (i == 1) { // If first gamma
              gamma_p(i) = rtnorm(gen,
                                  gamma[target](i-1),
                                  gamma[target](i+1),
                                  gamma[target](i),
                                  tune(target)
                                  ).first;
            } else { // If any other gamma
              gamma_p(i) = rtnorm(gen,
                                  gamma_p(i-1),
                                  gamma[target](i+1),
                                  gamma[target](i),
                                  tune(target)
                                  ).first;
          }//if first gamma
        }//for each category of target
      }//if-else ncat>2
      
      // loop over observations and construct the acceptance ratio
      double loglikerat = 0.0;
      Xbeta = X[target] * beta;
      for (unsigned int i=0; i<X[target].n_rows; ++i){
        int y_val = Y[target](i);
        if (y_val == ncat(target)){
          loglikerat = loglikerat
          + log(1.0  - gsl_cdf_ugaussian_P(gamma_p(y_val-1) - Xbeta[i]))
          - log(1.0 - gsl_cdf_ugaussian_P(gamma[target](y_val-1) - Xbeta[i]));
        }
        else if (y_val == 1){
          loglikerat = loglikerat + log(gsl_cdf_ugaussian_P(gamma_p(y_val) - Xbeta[i]))
          - log(gsl_cdf_ugaussian_P(gamma[target](y_val) - Xbeta[i]));
        }
        else{
          loglikerat = loglikerat
          + log(gsl_cdf_ugaussian_P(gamma_p(y_val) - Xbeta[i]) -
            gsl_cdf_ugaussian_P(gamma_p(y_val-1) - Xbeta[i]))
          - log(gsl_cdf_ugaussian_P(gamma[target](y_val) - Xbeta[i]) -
            gsl_cdf_ugaussian_P(gamma[target](y_val-1) - Xbeta[i]));
        }
      }
      new_gamma[target] = gamma_p; // new_gamma is not used.
      if (gsl_ran_flat(gen, 0.0, 1.0) <= exp(loglikerat)){
         gamma[target] = gamma_p;
        if (iter >= burnin) {
          ++accepts(target);
        }
      }
    }//for each target
    // if (gsl_ran_flat(gen, 0.0, 1.0) <= exp(loglikerat)){
    //   ++accepts(ntargets);
    //   gamma = new_gamma;
    // }
    
    // Step 2: Update Z: (Z | gamma, beta, y)
    offset = 0;
    for (unsigned int target = 0; target < ntargets; ++target) {
      if (target > 0) {
        offset = offset + sizes[target-1];
      }
      Xbeta = X[target] * beta;
      for (unsigned int i=0; i<X[target].n_rows; ++i){
        Z(offset+i) = rtnorm(gen,
          gamma[target](Y[target](i)-1),
          gamma[target](Y[target](i)),
          Xbeta[i], 
               1.0
        ).first;
      }
    }
    
    
    // Step 3: Update beta (beta | Z, gamma)
    arma::mat XpZ = arma::trans(Xall)*Z;
    beta = NormNormregress_beta_draw(gen, XpX, XpZ, meanPrior, precPrior, 1.0);
    
    // print output to stdout
    if(verbose > 0 && iter % verbose == 0){
      Rprintf("\n\noprobit_gibbs iteration %i of %i \n", (iter+1), tot_iter);
      Rprintf("beta = \n");
      // for (unsigned int target = 0; target < ntargets; ++target) {
      //   Rprintf("%10.5f\n", beta[0]-gamma[target](1));
      // }
      for (unsigned int j=0; j<k; ++j)
        Rprintf("%10.5f\n", beta[j]);
      for (unsigned int target = 0; target < ntargets; ++target) {
        Rprintf("Metropolis acceptance rate for gamma = %3.5f\n",
                static_cast<double>(accepts[target])/static_cast<double>(iter+1));
      }
      
    }
    
    // store values in matrices
    if (iter >= burnin && ((iter % thin)==0)){
      for (unsigned int j=0; j<k; ++j) {
        // Rcout << "["<< storemat.n_rows << ", "<< storemat.n_cols << "]\n";
        // Rcout << count << " " << j<< "\n";
        // sleep(1);
        storebeta(count, j) = beta[j];
      }
      for (unsigned int target = 0; target < ntargets; ++target) {
        for (unsigned int j=1; j<(ncat[target]); ++j){
          // Rcout << "["<< storemat.n_rows << ", "<< storemat.n_cols << "]\n";
          // Rcout << count << " " << j+k-2 << "\n";
          // sleep(1);
          storegamma[target](count, j-1) = gamma[target](j);
        }
      }
      ++count;
    }
    
  }
  
  // for (unsigned int target = 0; target < ntargets; ++target) {
  //   Rprintf("Metropolis acceptance rate for gamma = %3.5f\n",
  //           static_cast<double>(accepts[target])/static_cast<double>(tot_iter));
  // }
  
  Rcpp::List gammas = Rcpp::List::create();
  for (unsigned int target = 0; target < ntargets; ++target) {
    gammas.push_back(storegamma[target]);
  }
  
  return List::create(
    _["storebeta"] = storebeta,
    // _["store_beta"] = trans(store_beta),
    _["storegamma"] = gammas,
    _["accepts"] = accepts,
    _["total_iter"] = tot_iter
  );
}


/*** R
#//print(This is R code)
*/
