/**
 * @file hprobit_gibbs.cpp
 * @author Johan Falkenjack
 * @author Daniel Tufvesson
 * @version 1.0
 * @brief A Metropolis-Hastings-withing-Gibbs sampler for the Multi-Scale Probit model.
 */

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "sampling.hpp"

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List cpp_hprobit(
    const Rcpp::List& Xlist,
    const Rcpp::List& Ylist,
    const arma::colvec& meanPrior,
    const arma::mat& precPrior,
    const int fixZero,
    const arma::ivec& ncat,
    const Rcpp::List& gammaStart,
    const arma::colvec& betaStart,
    const arma::vec& tune,
    const int iterations,
    const int burnin,
    const int thin,
    const int seed,
    const int verbose
) {
  
    //--- GSL random init ---
    // Used to take efficient samples from truncated normal.
    gsl_rng_env_setup();                          // Read variable environnement
    const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
    gsl_rng *gen = gsl_rng_alloc (type);          // Rand generator allocation
    gsl_rng_set(gen, seed);
  
    // Unpack data and define constants.
    const Data data = unpack_data(Xlist, Ylist);
    const int ntargets = data.ntargets;
    const int tot_iter = iterations+burnin;
    const int nstore = iterations/thin;

    // Unpack gamma.
    std::vector<colvec> gamma(ntargets);
    for (unsigned int target = 0; target < ntargets; ++target) {
        gamma[target] = Rcpp::as<colvec>(gammaStart[target]);
        gamma[target](0) = -std::numeric_limits<double>::max();
        gamma[target](ncat[target]) = std::numeric_limits<double>::max();
    }
  
  
    // Storage matrices
    mat storebeta(nstore, betaStart.n_rows, arma::fill::zeros);
    std::vector<mat> storegamma(ntargets);
    for (unsigned int target = 0; target < ntargets; ++target) {
        storegamma[target] = mat(nstore, ncat(target)-1, arma::fill::zeros);
    }
  
    // Set starting points
    colvec beta = betaStart;
  
    // Set Z vector starting point as OLS estimates
    colvec Z(data.nobs, arma::fill::zeros);// = X * beta;
    colvec Xbeta;
  
    // Bookkeeping
    unsigned int count = 0;
    ivec accepts(ntargets+1, arma::fill::zeros);
    int offset;

    // Gibbs loop
    for (unsigned int iter = 0; iter<tot_iter; ++iter) {

        // Step 1: Update gammas: (gamma | u, beta)
        // Cowles update of gamma
        for (unsigned int t = 0; t < ntargets; ++t) {
            int target = (iter+t) % ntargets;

            // Make proposal for gamma.
            bool accepted = mh_update_gamma(
                gamma[target],
                beta,
                data.X[target],
                data.Y[target],
                ncat(target),
                tune(target),
                1.0, // inv_temperature is 1.0 for the Gibbs sampler.
                gen
            );
            if (accepted && iter >= burnin) {
                ++accepts(target);
            }
        }
    
        // Step 2: Update Z: (Z | gamma, beta, y)
        offset = 0;
        for (unsigned int target = 0; target < ntargets; ++target) {
            if (target > 0) {
                offset = offset + data.target_nobs[target-1];
            }
            Xbeta = data.X[target] * beta;
            for (unsigned int i=0; i<data.X[target].n_rows; ++i){
                Z(offset+i) = rtnorm(gen,
                    gamma[target](data.Y[target](i)-1),
                    gamma[target](data.Y[target](i)),
                    Xbeta[i], 
                    1.0
                ).first;
            }
        }
    
    
        // Step 3: Update beta (beta | Z, gamma)
        arma::mat XpZ = arma::trans(data.Xall)*Z;
        beta = mspm_util::NormNormregress_beta_draw(gen, data.XpX, XpZ, meanPrior, precPrior, 1.0);
    
        // print output to stdout
        if(verbose > 0 && iter % verbose == 0){
            Rprintf("\n\noprobit_gibbs iteration %i of %i \n", (iter+1), tot_iter);
            Rprintf("beta = \n");
            for (unsigned int j=0; j < data.npredictors; ++j){
                Rprintf("%10.5f\n", beta[j]);
            }
            for (unsigned int target = 0; target < ntargets; ++target) {
                Rprintf("Metropolis acceptance rate for gamma = %3.5f\n",
                        static_cast<double>(accepts[target])/static_cast<double>(iter+1));
            }
        }
    
        // store values in matrices
        if (iter >= burnin && ((iter % thin)==0)){
            for (unsigned int j=0; j < data.npredictors; ++j) {
                storebeta(count, j) = beta[j];
            }
            for (unsigned int target = 0; target < ntargets; ++target) {
                for (unsigned int j=1; j<(ncat[target]); ++j){
                    storegamma[target](count, j-1) = gamma[target](j);
                }
            }
            ++count;
        }
    }
  
    Rcpp::List gammas = Rcpp::List::create();
    for (unsigned int target = 0; target < ntargets; ++target) {
        gammas.push_back(storegamma[target]);
    }
  
    return List::create(
        _["storebeta"] = storebeta,
        _["storegamma"] = gammas,
        _["accepts"] = accepts,
        _["total_iter"] = tot_iter
    );
}

// colvec propose_gamma(
//     const colvec& gamma,
//     int target,
//     double sigma,
//     int ncategories,
//     gsl_rng *gen
// ) {
//     colvec gamma_p(ncategories+1, arma::fill::zeros);
//     gamma_p.head(1) = -INFINITY;
//     gamma_p.tail(1) = INFINITY;
//     if (ncategories == 2) {
//         // Draw new split point
//         gamma_p(1) = rtnorm(
//             gen,
//             -INFINITY,
//             INFINITY,
//             gamma(1),
//             sigma
//         ).first;
//     } 
//     else {
//         for (int i=1; i < ncategories; ++i){
//             if (i == 1) { // If first gamma
//                 gamma_p(i) = rtnorm(gen,
//                     gamma(i-1),
//                     gamma(i+1),
//                     gamma(i),
//                     sigma
//                 ).first;
//             } 
//             else { // If any other gamma
//                 gamma_p(i) = rtnorm(gen,
//                     gamma_p(i-1),
//                     gamma(i+1),
//                     gamma(i),
//                     sigma
//                 ).first;
//             }//if first gamma
//         }//for each category of target
//     }//if-else ncat>2
//     return gamma_p;
// }
