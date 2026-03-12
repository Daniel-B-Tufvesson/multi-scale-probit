
#include "sampling.hpp"
// [[Rcpp::depends(RcppArmadillo)]]

Data unpack_data(
    const Rcpp::List& xlist,
    const Rcpp::List& ylist
) {
    const int ntargets = xlist.size();
    Data data {
        .X = std::vector<arma::mat>(ntargets),
        .Y = std::vector<arma::colvec>(ntargets),
        .Xall = arma::mat(),
        .XpX = arma::mat(),
        .target_nobs = arma::ivec(ntargets, arma::fill::zeros),
        .nobs = 10,
        .ntargets = ntargets,
        .npredictors = 10
    };
    for (unsigned int target = 0; target < ntargets; ++target) {
        data.X[target] = Rcpp::as<arma::mat>(xlist[target]);
        data.Y[target] = Rcpp::as<arma::colvec>(ylist[target]);
        if (target == 0) {
            data.Xall = data.X[target];
        } 
        else {
            data.Xall = arma::join_vert(data.Xall, data.X[target]);
        }
        data.target_nobs[target] = data.Y[target].n_elem;
    }
    data.XpX = mspm_util::crossprod(data.Xall, data.Xall);
    data.nobs = data.Xall.n_rows;
    data.npredictors = data.Xall.n_cols;
    return data;
}

std::vector<colvec> unpack_gamma(const Rcpp::List& gammaStart, const arma::ivec& ncat) {
    int ntargets = gammaStart.size();
    std::vector<colvec> gamma(ntargets);
    for (int target = 0; target < ntargets; ++target) {
        gamma[target] = Rcpp::as<colvec>(gammaStart[target]);
        gamma[target](0) = -std::numeric_limits<double>::max();
        gamma[target](ncat[target]) = std::numeric_limits<double>::max();
    }
    return gamma;
}


void compute_acceptance_rate(
    const arma::vec& acceptance_probabilities,
    int nsamples,
    arma::vec& acceptance_rates
) {
    for (unsigned int target = 0; target < acceptance_probabilities.n_rows; ++target) {
        acceptance_rates(target) = acceptance_probabilities(target) / static_cast<double>(nsamples);
    }
}



