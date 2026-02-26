
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
        .target_nobs = std::vector<int>(ntargets),
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
) {
    double log_likelihood_ratio = 0.0;
    colvec y_star = X * beta;
    auto cdf = gsl_cdf_ugaussian_P;
    for (unsigned int i = 0; i < X.n_rows; ++i){
        int y_val = Y(i);
        // Handle last category.
        if (y_val == ncategories){
            log_likelihood_ratio = log_likelihood_ratio
            + log(1.0 - cdf(gamma_p(y_val-1) - y_star[i]))
            - log(1.0 - cdf(gamma(y_val-1) - y_star[i]));
        }
        // Handle first category.
        else if (y_val == 1){
            log_likelihood_ratio = log_likelihood_ratio 
            + log(cdf(gamma_p(y_val) - y_star[i]))
            - log(cdf(gamma(y_val) - y_star[i]));
        }
        // Handle categories inbetween.
        else {
            log_likelihood_ratio = log_likelihood_ratio
            + log(cdf(gamma_p(y_val) - y_star[i]) - 
                  cdf(gamma_p(y_val-1) - y_star[i]))
            - log(cdf(gamma(y_val) - y_star[i]) - 
                  cdf(gamma(y_val-1) - y_star[i]));
        }
    }
    return log_likelihood_ratio;
}

/**
 * Compute the log of the normalization constant for the truncated Gaussian distribution used 
 * in the proposal distribution for the gammas.
 * 
 * @param x The value at which to compute the normalization constant.
 * @param a The lower truncation point for the truncated Gaussian distribution.
 * @param b The upper truncation point for the truncated Gaussian distribution.
 * @param sigma The standard deviation of the truncated Gaussian distribution.
 * @return The log of the normalization constant for the truncated Gaussian distribution at the 
 * given value x.
 */
double compute_log_trunc_gauss_norm_constant(double x, double a, double b, double sigma) {
    auto cdf = gsl_cdf_ugaussian_P;
    double alpha = (a - x) / sigma;
    double beta = (b - x) / sigma;
    return log(cdf(beta) - cdf(alpha));
}

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
) {
    double log_proposal_ratio = 0;
    for (int i = 1; i < ncategories; i++) {
        // Normalization term Z(gamma) 
        double z_old = compute_log_trunc_gauss_norm_constant(
            gamma(i), 
            gamma(i-1), 
            gamma(i+1), 
            sigma
        );

        if (i == 1) { // If first gamma
            // Normalization term Z(gamma_p)
            double z_prop = compute_log_trunc_gauss_norm_constant(
                gamma_p(i), 
                gamma(i-1), 
                gamma(i+1), 
                sigma
            );
            
            log_proposal_ratio += z_old - z_prop;
        } 
        else { // If any other gamma
            // Normalization term Z(gamma_p)
            double z_prop = compute_log_trunc_gauss_norm_constant(
                gamma_p(i), 
                gamma_p(i-1), 
                gamma(i+1), 
                sigma
            );
            log_proposal_ratio += z_old - z_prop;
        }
    }
    return log_proposal_ratio;
}

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
) {
    arma::colvec gamma_prop(ncategories + 1, arma::fill::zeros);
    gamma_prop.head(1) = -INFINITY;
    gamma_prop.tail(1) = INFINITY;
    // Note: probabilities are 0 for the edge thresholds.

    if (ncategories == 2) {
        // Draw new split point
        gamma_prop(1) = rtnorm(
            rng,
            -INFINITY,
            INFINITY,
            gamma(1),
            sigma
        ).first;
    } 
    else {
        for (unsigned int i = 1; i < ncategories; i++) {
            if (i == 1) { // If first gamma
                gamma_prop(i) = rtnorm(
                    rng,
                    gamma(i-1),
                    gamma(i+1),
                    gamma(i),
                    sigma
                ).first;
            } 
            else { // If any other gamma
                gamma_prop(i) = rtnorm(
                    rng,
                    gamma_prop(i-1),    // a
                    gamma(i+1), // b
                    gamma(i),   //mu
                    sigma
                ).first;
            }
        }
    }
    return gamma_prop;
}


bool mh_update_gamma(
    colvec& gamma,
    const colvec& beta,
    const mat& X,
    const colvec& Y,
    int ncategories,
    double sigma,
    double inv_temperature,
    gsl_rng* rng
) {
    // Make proposal.
    arma::colvec gamma_prop = propose_gamma(
        gamma,
        ncategories,
        sigma,
        rng
    );
    double log_likelihood_ratio = compute_log_likelihood_ratio(
        gamma,
        gamma_prop,
        X,
        Y,
        beta,
        ncategories
    );
    double log_proposal_ratio = compute_gamma_log_proposal_ratio(
        gamma,
        gamma_prop,
        ncategories,
        sigma
    );

    double log_accept_ratio = inv_temperature * log_likelihood_ratio + log_proposal_ratio;

    // Do MH acceptance step.
    if (gsl_ran_flat(rng, 0.0, 1.0) <= exp(log_accept_ratio)) {
        gamma = gamma_prop;
        return true;
    }
    return false;
}

void gibbs_update_beta(
    colvec& beta,
    const std::vector<colvec>& gamma,
    const Data& data,
    const arma::colvec& beta_mean_prior,
    const arma::mat& beta_prec_prior,
    gsl_rng* rng
) {
    // Draw latent y*.
    arma::colvec ystar = arma::colvec(data.nobs, arma::fill::zeros);
    int offset = 0;
    for (unsigned int target = 0; target < data.ntargets; target++) {
        if (target > 0) {
            int nobs = data.Y[target-1].n_elem;
            offset += nobs;
        }
        const colvec current_ystar = data.X[target] * beta;
        for (unsigned int i = 0; i < data.X[target].n_rows; i++) {
            ystar(offset + i) = rtnorm(
                rng,
                gamma[target](data.Y[target](i) - 1),
                gamma[target](data.Y[target](i)),
                current_ystar[i], 
                1.0
            ).first;
        }
    }

    // Draw new beta.
    arma::mat XpZ = arma::trans(data.Xall) * ystar;
    beta = mspm_util::NormNormregress_beta_draw(
        rng, 
        data.XpX, 
        XpZ, 
        beta_mean_prior, 
        beta_prec_prior, 
        1.0
    );
}