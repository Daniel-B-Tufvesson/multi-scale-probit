# Load experiment data from "experiment_gibbs.R" and analyze it.

source("scripts/experiment/diagnostics.R")

cv_res = readRDS("results/cv_gibbs_5000.rds")
print("Loaded results")


# Unpack burnin betas from all splits into a single matrix.
n_splits = length(cv_res$allBurninBetas)
n_samples = nrow(cv_res$allBurninBetas[[1]])
n_beta_params = ncol(cv_res$allBurninBetas[[1]])
all_burnin_betas <- lapply(1:n_beta_params, function(j) {
    do.call(cbind, lapply(cv_res$allBurninBetas, function(mat) mat[, j, drop = FALSE]))
})

# Unpack burnin gammas from all the splits into a list matrixes (grouped by target).
n_targets <- length(cv_res$allBurninGammas[[1]])
n_gammas <- lapply(1:n_targets, function(target_idx) ncol(cv_res$allBurninGammas[[1]][[target_idx]]))
all_burnin_gammas <- lapply(1:n_targets, function(target_idx) {
    lapply(1:n_gammas[[target_idx]], function(gamma_idx) {
        do.call(cbind, lapply(cv_res$allBurninGammas, function(split_list) split_list[[target_idx]][, gamma_idx, drop = FALSE]))
    })
})

# Unpack sampled betas.
all_sampled_betas <- lapply(1:n_beta_params, function(j) {
    do.call(cbind, lapply(cv_res$allBetas, function(mat) mat[, j, drop = FALSE]))
})

# Unpack sampled gammas.
all_sampled_gammas <- lapply(1:n_targets, function(target_idx) {
    lapply(1:n_gammas[[target_idx]], function(gamma_idx) {
        do.call(cbind, lapply(cv_res$allGammas, function(split_list) split_list[[target_idx]][, gamma_idx, drop = FALSE]))
    })
})


###################################################################################################
# Compute the cumulative R-hat.
rhat_for_every = 50
cum_rhat_burnin_beta <- list()
for (j in 1:n_beta_params) {
    cum_rhat_burnin_beta[[j]] <- cumulative_gelman_rubin_rhat(all_burnin_betas[[j]], for_every = rhat_for_every)
}

cum_rhat_burnin_gammas <- list()
for (j in 1:n_targets) {
    cum_rhat <- list()
    for (i in 1:n_gammas[[j]]) {
        cum_rhat[[i]] <- cumulative_gelman_rubin_rhat(all_burnin_gammas[[j]][[i]], for_every = rhat_for_every)
    }
    cum_rhat_burnin_gammas[[j]] <- cum_rhat
}

# Compute regular R-hat for sampled values.
rhat_sampled_betas <- list()
for (j in 1:n_beta_params) {
    rhat_sampled_betas[[j]] <- gelman_rubin_rhat_single_param(all_sampled_betas[[j]])
}
rhat_sampled_gammas <- list()
for (j in 1:n_targets) {
    rhat_gammas <- list()
    for (i in 1:n_gammas[[j]]) {
        rhat_gammas[[i]] <- gelman_rubin_rhat_single_param(all_sampled_gammas[[j]][[i]])
    }
    rhat_sampled_gammas[[j]] <- rhat_gammas
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

##################################################################################################
# Compute min ESS/second.
beta_min_ess_per_second <- list() # One for each split.
for (j in 1:n_splits) {
    ess_values <- effective_sample_size(cv_res$allBetas[[j]], method = "acf")
    time_seconds <- cv_res$samplingTime[j]
    beta_min_ess_per_second[[j]] <- min(ess_values) / time_seconds
}

###################################################################################################
# Some helper functions.
is_rhat_ok <- function(rhat) {
    if (rhat <= 1.01) {
        return("✅")
    } else if (rhat <= 1.05) {
        return("🟡")
    } else {
        return("❌")
    }
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
        rhat_ok <- is_rhat_ok(final_rhat)
        cat(paste0("   - ", rhat_ok, " Beta ", j, " converged at ", converged_at, "\n"))
    }

    cat("2. Cumulative R-hat Analysis for Burn-in Gammas:\n")
    for (j in 1:n_targets) {
        cat(paste0("   - Gammas for target ", j, ":\n"))
        for (i in 1:n_gammas[[j]]) {
            final_rhat <- tail(cum_rhat_burnin_gammas[[j]][[i]], n = 1)
            converged_at <- x_vals[which(cum_rhat_burnin_gammas[[j]][[i]] <= 1.05)[1]]
            converged_at <- ifelse(is.na(converged_at), "N/A", converged_at)
            rhat_ok <- is_rhat_ok(final_rhat)
            cat(paste0("      - ", rhat_ok, " Gamma ", i, " converged at ", converged_at, "\n"))
        }
    }
    cat("3. R-hat Analysis for Sampled Betas\n")
    for (j in 1:n_beta_params) {
        final_rhat <- rhat_sampled_betas[[j]]
        rhat_ok <- is_rhat_ok(final_rhat)
        cat(paste0("   - ", rhat_ok, " Beta ", j, ": R-hat = ", round(final_rhat, 4), "\n"))
    }

    cat("4. R-hat Analysis for Sampled Gammas\n")
    for (j in 1:n_targets) {
        cat(paste0("   - Gammas for target ", j, ":\n"))
        for (i in 1:n_gammas[[j]]) {
            final_rhat <- rhat_sampled_gammas[[j]][[i]]
            rhat_ok <- is_rhat_ok(final_rhat)
            cat(paste0("      - ", rhat_ok, " Gamma ", i, ": R-hat = ", round(final_rhat, 4), "\n"))
        }
    }

    cat("5. Minimum ESS per second for Sampled Betas\n")
    for (j in 1:n_splits) {
        cat(paste0("   - Split ", j, ": ", round(beta_min_ess_per_second[[j]], 4), " ESS/s\n"))
    }
}

print_diagnostic_report()