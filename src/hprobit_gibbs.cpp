/**
 * @file hprobit_gibbs.cpp
 * @author Johan Falkenjack
 * @author Daniel Tufvesson
 * @version 2.0
 * @brief A Metropolis-Hastings-withing-Gibbs sampler for the Multi-Scale Probit model.
 */

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "sampling.hpp"
#include <chrono>

void do_step(
    int iter,
    int ntargets,
    const arma::ivec& ncat,
    colvec& beta,
    std::vector<colvec>& gamma,
    const arma::vec& tune,
    const Data& data,
    ivec& accepts,
    const arma::colvec& meanPrior,
    const arma::mat& precPrior,
    gsl_rng* gen
) {
    // Step 1: Update gammas with Metropolis-Hastings.
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
        if (accepted) {
            ++accepts(target);
        }
    }

    // Step 2: Update beta with Gibbs.
    gibbs_update_beta(
        beta,
        gamma,
        data,
        meanPrior,
        precPrior,
        gen
    );
}

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
    const bool save_burnin_samples,
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
    // const int tot_iter = iterations+burnin;
    const int nstore = iterations/thin;
    const int nstore_burnin = save_burnin_samples ? burnin : 0;

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

    // Storage matrices for burnin.
    mat storebeta_burnin(nstore_burnin, betaStart.n_rows, arma::fill::zeros);
    std::vector<mat> storegamma_burnin(ntargets);
    if (save_burnin_samples) {
        for (unsigned int target = 0; target < ntargets; ++target) {
            storegamma_burnin[target] = mat(burnin, ncat(target)-1, arma::fill::zeros);
        }
    }
  
    // Set starting points
    colvec beta = betaStart;
  
    // Set Z vector starting point as OLS estimates
    colvec Z(data.nobs, arma::fill::zeros);// = X * beta;
    colvec Xbeta;
  
    // Bookkeeping
    unsigned int count = 0;
    ivec accepts(ntargets+1, arma::fill::zeros);

    // Measure burnin sampling time.
    auto start_time_burnin = std::chrono::high_resolution_clock::now();

    // Burnin loop.
    for (int iter = 0; iter < burnin; iter++) {
        do_step(
            iter,
            ntargets,
            ncat,
            beta,
            gamma,
            tune,
            data,
            accepts,
            meanPrior,
            precPrior,
            gen
        );

        // Print progress.
        if(verbose > 0 && iter % verbose == 0){
            double avg_acceptance_rate = 0.0;
            Rcpp::Rcout << "\n\noprobit_gibbs burnin iteration " << (iter+1) 
                << " of " << burnin << std::endl;

            // Print acceptance rates for each target.
            for (unsigned int target = 0; target < ntargets; ++target) {
                double acceptance_rate = static_cast<double>(accepts[target]) / static_cast<double>(iter+1);
                avg_acceptance_rate += acceptance_rate;
                Rcpp::Rcout << "Metropolis acceptance rate for gamma (target " << target << ") = " 
                    << acceptance_rate << std::endl;
            }
        }

        // Store burnin sample.
        if (save_burnin_samples && ((iter % thin)==0)) {
            for (unsigned int j=0; j < data.npredictors; ++j) {
                storebeta_burnin(count, j) = beta[j];
            }
            for (unsigned int target = 0; target < ntargets; ++target) {
                for (unsigned int j=1; j<(ncat[target]); ++j){
                    storegamma_burnin[target](count, j-1) = gamma[target](j);
                }
            }
            ++count;
        }
    }
    // Measure burnin time in seconds.
    auto end_time_burnin = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_burnin = end_time_burnin - start_time_burnin;

    // Reset bookkeeping for main sampling loop.
    count = 0;
    accepts.fill(0);

    // Measure sampling time.
    auto start_time = std::chrono::high_resolution_clock::now();

    // Sampling loop.
    for (int iter = 0; iter<iterations; ++iter) {

        do_step(
            iter,
            ntargets,
            ncat,
            beta,
            gamma,
            tune,
            data,
            accepts,
            meanPrior,
            precPrior,
            gen
        );
        
        // Print progress.
        if(verbose > 0 && iter % verbose == 0){
            double avg_acceptance_rate = 0.0;
            Rcpp::Rcout << "\n\noprobit_gibbs iteration " << (iter+1) 
                << " of " << iterations << std::endl;

            // Print acceptance rates for each target.
            for (unsigned int target = 0; target < ntargets; ++target) {
                double acceptance_rate = static_cast<double>(accepts[target]) / static_cast<double>(iter+1);
                avg_acceptance_rate += acceptance_rate;
                Rcpp::Rcout << "Metropolis acceptance rate for gamma (target " << target << ") = " 
                    << acceptance_rate << std::endl;
            }
        }
    
        // Store sample in matrices
        if ((iter % thin)==0) {
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

    // Measure sampling time in seconds.
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
  
    // Pack stored gammas.
    Rcpp::List gammas = Rcpp::List::create();
    for (unsigned int target = 0; target < ntargets; ++target) {
        gammas.push_back(storegamma[target]);
    }

    // Pack stored burnin gammas.
    Rcpp::List storegamma_burnin_list = Rcpp::List::create();
    if (save_burnin_samples) {
        for (unsigned int target = 0; target < ntargets; ++target) {
            storegamma_burnin_list.push_back(storegamma_burnin[target]);
        }
    }
  
    return List::create(
        _["storebeta"] = storebeta,
        _["storegamma"] = gammas,
        _["accepts"] = accepts,
        _["total_iter"] = iterations + burnin,
        _["storebeta_burnin"] = storebeta_burnin,
        _["storegamma_burnin"] = storegamma_burnin_list,
        _["sampling_time"] = elapsed.count(),
        _["burnin_time"] = elapsed_burnin.count()
    );
}
