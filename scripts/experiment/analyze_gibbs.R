# Load experiment data from "experiment_gibbs.R" and analyze it.

source("scripts/experiment/diagnostics.R")

cv_res = readRDS("results/cv_gibbs_5000.rds")
print("Loaded results")


# Unpack betas from all splits into a single matrix.
n_splits = length(cv_res$allBurninBetas)
n_samples = nrow(cv_res$allBurninBetas[[1]])
n_beta_params = ncol(cv_res$allBurninBetas[[1]])
all_betas <- lapply(1:n_beta_params, function(j) {
    do.call(cbind, lapply(cv_res$allBurninBetas, function(mat) mat[, j, drop = FALSE]))
})

# Unpack gammas from all the splits into a list matrixes (grouped by target).
all_gammas <- list()
n_targets <- length(cv_res$allBurninGammas[[1]])
n_gammas <- lapply(1:n_targets, function(target_idx) ncol(cv_res$allBurninGammas[[1]][[target_idx]]))
all_gammas <- lapply(1:n_targets, function(target_idx) {
    lapply(1:n_gammas[[target_idx]], function(gamma_idx) {
        do.call(cbind, lapply(cv_res$allBurninGammas, function(split_list) split_list[[target_idx]][, gamma_idx, drop = FALSE]))
    })
})

# for (j in 1:n_targets) {
#     gammas_j <- list()
#     for (i in 1:n_splits) {
#         gammas_j[[i]] <- cv_res$allBurninGammas[[i]][[j]]
#     }
#     all_gammas[[j]] <- gammas_j
#     n_gammas[j] = ncol(gammas_j)
# }

###################################################################################################
# Compute the cumulative R-hat.
rhat_for_every = 50
cum_rhat_burnin_beta <- list()
for (j in 1:n_beta_params) {
    cum_rhat_burnin_beta[[j]] <- cumulative_gelman_rubin_rhat(all_betas[[j]], for_every = rhat_for_every)
}

cum_rhat_burnin_gammas <- list()
for (j in 1:n_targets) {
    cum_rhat <- list()
    for (i in 1:n_gammas[[j]]) {
        cum_rhat[[i]] <- cumulative_gelman_rubin_rhat(all_gammas[[j]][[i]], for_every = rhat_for_every)
    }
    cum_rhat_burnin_gammas[[j]] <- cum_rhat
}

##################################################################################################
# Plot rhat values.
plot_cum_rhat <- function(){
    x_vals <- seq_along(cum_rhat_burnin_beta[[1]]) * rhat_for_every
    plot(x_vals, cum_rhat_burnin_beta[[1]], type = "l", ylim = c(1, 1.2), 
        xlab = "Cumulative Samples", ylab = "R-hat", main = "Cumulative R-hat for Burn-in Beta")
    abline(h = 1.05, col = "red", lty = 2, lwd = 2)  # Not converged.
    abline(h = 1.01, col = "orange", lty = 2, lwd = 2)  # Maybe converged.
}


###################################################################################################
# Print the final diagnostic report.
print_diagnostic_report <- function() {
    cat("Diagnostic Report for Gibbs Sampling:\n")
    cat("1. Cumulative R-hat Analysis for Burn-in Beta:\n")
    for (j in 1:n_beta_params) {
        final_rhat <- tail(cum_rhat_burnin_beta[[j]], n = 1)
        converged_at <- x_vals[which(cum_rhat_burnin_beta[[j]] <= 1.05)[1]]
        converged_at <- ifelse(is.na(converged_at), "N/A", converged_at)
        if (final_rhat <= 1.01) {
            rhat_ok <- "✅"
        } 
        else if (final_rhat <= 1.05) {
            rhat_ok <- "🟡"
        } 
        else {
            rhat_ok <- "❌"
        }
        cat(paste0("   - ", rhat_ok, " Beta ", j, " converged at ", converged_at, "\n"))
    }

    cat("2. Cumulative R-hat Analysis for Burn-in Gammas:\n")
    for (j in 1:n_targets) {
        cat(paste0("   - Gammas for target ", j, ":\n"))
        for (i in 1:n_gammas[[j]]) {
            final_rhat <- tail(cum_rhat_burnin_gammas[[j]][[i]], n = 1)
            converged_at <- x_vals[which(cum_rhat_burnin_gammas[[j]][[i]] <= 1.05)[1]]
            converged_at <- ifelse(is.na(converged_at), "N/A", converged_at)
            if (final_rhat <= 1.01) {
                rhat_ok <- "✅"
            } 
            else if (final_rhat <= 1.05) {
                rhat_ok <- "🟡"
            } 
            else {
                rhat_ok <- "❌"
            }
            cat(paste0("      - ", rhat_ok, " Gamma ", i, " converged at ", converged_at, "\n"))
        }
    }
}

print_diagnostic_report()