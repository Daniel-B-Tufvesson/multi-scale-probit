
#include "sampling.hpp"

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
) {
    const int ntargets = xlist.size();
    Data data {
        .X = std::vector<arma::mat>(ntargets),
        .Y = std::vector<arma::colvec>(ntargets),
        .Xall = arma::mat(),
        .XpX = arma::mat(),
        .target_nobs = std::vector<int>(ntargets),
        .nobs = 10
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
    return data;
}