# Load experiment data from "experiment_gibbs.R" and analyze it.

source("scripts/experiment/diagnostics.R")
source("scripts/experiment/plotting.R")

name <- "gibbs_5000"
cv_res = readRDS(paste0("results/cv_", name, ".rds"))
print("Loaded results")


# Unpack burnin betas from all splits into a single matrix.
n_splits = length(cv_res$allBurninBetas)
n_burnin_sample = nrow(cv_res$allBurninBetas[[1]])
n_beta_params = ncol(cv_res$allBurninBetas[[1]])
all_burnin_betas <- lapply(1:n_beta_params, function(j) {
    do.call(cbind, lapply(cv_res$allBurninBetas, function(mat) mat[, j, drop = FALSE]))
})

# Unpack burnin gammas from all the splits into a list matrixes (grouped by target).
n_targets <- length(cv_res$allBurninGammas[[1]])
n_gammas <- lapply(1:n_targets, function(target_idx) ncol(cv_res$allBurninGammas[[1]][[target_idx]]))
n_gamma_params <- sum(unlist(n_gammas))
all_burnin_gammas <- lapply(1:n_targets, function(target_idx) {
    lapply(1:n_gammas[[target_idx]], function(gamma_idx) {
        do.call(cbind, lapply(cv_res$allBurninGammas, function(split_list) split_list[[target_idx]][, gamma_idx, drop = FALSE]))
    })
})

# Unpack sampled betas.
n_samples <- nrow(cv_res$allBetas[[1]])
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
# Compute the cumulative R-hat for burn-in samples.
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

# Compute convergence time for burn-in beta samples.
burnin_convergence_time_beta <- list()
n_burnin_beta_convergences = 0 # Number of betas that converged.
for (j in 1:n_beta_params) {
    convergence_time <- first_convergence_index(cum_rhat_burnin_beta[[j]], step_size = rhat_for_every)
    convergence_time_beta[[j]] <- convergence_time

    if (!is.na(convergence_time)) {
        n_burnin_beta_convergences <- n_burnin_beta_convergences + 1
    }
}
# Compute convergence time for burn-in beta samples.
burnin_convergence_time_gammas <- list()
n_burnin_gamma_convergences = 0 # Number of gammas that converged.
for (j in 1:n_targets) {
    convergence_times <- list()
    for (i in 1:n_gammas[[j]]) {
        convergence_time <- first_convergence_index(cum_rhat_burnin_gammas[[j]][[i]], step_size = rhat_for_every)
        convergence_times[[i]] <- convergence_time

        if (!is.na(convergence_time)) {
            n_burnin_gamma_convergences <- n_burnin_gamma_convergences + 1
        }
    }
    burnin_convergence_time_gammas[[j]] <- convergence_times
}

# Compute median convergence time for burn-in betas and gammas.
median_convergence_time_beta <- median(unlist(convergence_time_beta), na.rm = TRUE)
median_convergence_time_gammas <- median(unlist(burnin_convergence_time_gammas), na.rm = TRUE)

# Compute regular R-hat for sampled values.
rhat_sampled_betas <- list()
n_sampled_betas_convergences = 0 # Number of sampled betas that converged.
for (j in 1:n_beta_params) {
    rhat_sampled_betas[[j]] <- gelman_rubin_rhat_single_param(all_sampled_betas[[j]])

    if (rhat_sampled_betas[[j]] <= 1.01) {
        n_sampled_betas_convergences <- n_sampled_betas_convergences + 1
    }
}
rhat_sampled_gammas <- list()
n_sampled_gammas_convergences = 0 # Number of sampled gammas that converged.
for (j in 1:n_targets) {
    rhat_gammas <- list()
    for (i in 1:n_gammas[[j]]) {
        rhat_gammas[[i]] <- gelman_rubin_rhat_single_param(all_sampled_gammas[[j]][[i]])

        if (rhat_gammas[[i]] <= 1.01) {
            n_sampled_gammas_convergences <- n_sampled_gammas_convergences + 1
        }
    }
    rhat_sampled_gammas[[j]] <- rhat_gammas
}

##################################################################################################
# Plot convergence time distribution for burn-in betas and gammas.
density_plot(unlist(convergence_time_beta),
             title = "Density of Convergence Times for Burn-in Betas", 
             xlabel = "Convergence Time (samples)")
density_plot(unlist(burnin_convergence_time_gammas),
             title = "Density of Convergence Times for Burn-in Gammas", 
             xlabel = "Convergence Time (samples)")

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
# Compute min ESS
beta_min_ess <- list() # One for each split.
beta_min_ess_per_second <- list() 
beta_min_ess_per_iteration <- list()
for (j in 1:n_splits) {
    ess_values <- effective_sample_size(cv_res$allBetas[[j]], method = "acf")
    time_seconds <- cv_res$samplingTime[j]
    min_ess <- min(ess_values)
    beta_min_ess[[j]] <- min_ess
    beta_min_ess_per_second[[j]] <- min_ess / time_seconds
    beta_min_ess_per_iteration[[j]] <- min_ess / n_samples
}

gamma_min_ess <- list() # One for each split.
gamma_min_ess_per_second <- list()
gamma_min_ess_per_iteration <- list()
for (j in 1:n_splits) {
    gammas_split <- cv_res$allGammas[[j]]
    ess_values <- unlist(lapply(gammas_split, function(gamma_mat) effective_sample_size(gamma_mat, method = "acf")))
    time_seconds <- cv_res$samplingTime[j]
    min_ess <- min(ess_values)
    gamma_min_ess[[j]] <- min_ess
    gamma_min_ess_per_second[[j]] <- min_ess / time_seconds
    gamma_min_ess_per_iteration[[j]] <- min_ess / n_samples
}

# Cpmpute median min ESS.
median_min_ess_beta <- median(unlist(beta_min_ess))
median_min_ess_gamma <- median(unlist(gamma_min_ess))
median_min_ess_per_second_beta <- median(unlist(beta_min_ess_per_second))
median_min_ess_per_second_gamma <- median(unlist(gamma_min_ess_per_second))
median_min_ess_per_iteration_beta <- median(unlist(beta_min_ess_per_iteration))
median_min_ess_per_iteration_gamma <- median(unlist(gamma_min_ess_per_iteration))

###################################################################################################
# Plot distributions of ESS/per second.
density_plot(unlist(beta_min_ess_per_second),
             title = "Density of Minimum ESS/s for Sampled Betas", 
             xlabel = "Minimum ESS/s")

density_plot(unlist(gamma_min_ess_per_second),
             title = "Density of Minimum ESS/s for Sampled Gammas", 
             xlabel = "Minimum ESS/s")

# Plot distributions of ESS/iteration.
density_plot(unlist(beta_min_ess_per_iteration),
             title = "Density of Minimum ESS/iteration for Sampled Betas", 
             xlabel = "Minimum ESS/iteration")

density_plot(unlist(gamma_min_ess_per_iteration),
             title = "Density of Minimum ESS/iteration for Sampled Gammas", 
             xlabel = "Minimum ESS/iteration")

###################################################################################################
# Save diagnostic results.
diagnostic_results <- list(
    cum_rhat_burnin_beta = cum_rhat_burnin_beta,
    cum_rhat_burnin_gammas = cum_rhat_burnin_gammas,
    burnin_convergence_time_beta = convergence_time_beta,
    burnin_convergence_time_gammas = burnin_convergence_time_gammas,
    rhat_sampled_betas = rhat_sampled_betas,
    rhat_sampled_gammas = rhat_sampled_gammas,
    beta_min_ess = beta_min_ess,
    gamma_min_ess = gamma_min_ess,
    beta_min_ess_per_second = beta_min_ess_per_second,
    gamma_min_ess_per_second = gamma_min_ess_per_second,
    beta_min_ess_per_iteration = beta_min_ess_per_iteration,
    gamma_min_ess_per_iteration = gamma_min_ess_per_iteration
)
saveRDS(diagnostic_results, file = paste0("results/diagnostics_", name, ".rds"))
cat("Saved diagnostic results to ", paste0("results/diagnostics_", name, ".rds"), "\n")

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

to_percent <- function(value) {
    return(paste0(round(value * 100, 2), "%"))
}


###################################################################################################
# Print the final diagnostic report.
print_diagnostic_report <- function() {
    cat("Diagnostic Report for Gibbs Sampling:\n")
    cat("1. Cumulative R-hat Analysis for Burn-in Beta:\n")
    cat("   - ", n_burnin_beta_convergences, "/", n_beta_params, "(",
        to_percent(n_burnin_beta_convergences / n_beta_params), ") burn-in betas converged.\n")
    for (j in 1:n_beta_params) {
        final_rhat <- tail(cum_rhat_burnin_beta[[j]], n = 1)
        converged_at <- convergence_time_beta[[j]]
        rhat_ok <- is_rhat_ok(final_rhat)
        cat(paste0("   - ", rhat_ok, " Beta ", j, " converged at ", converged_at, "\n"))
    }

    cat("2. Cumulative R-hat Analysis for Burn-in Gammas:\n")
    cat(paste0("      - ", n_burnin_gamma_convergences, "/", n_gamma_params, " (",
            to_percent(n_burnin_gamma_convergences / n_gamma_params), ") burn-in gammas converged.\n"))
    for (j in 1:n_targets) {
        cat(paste0("   - Gammas for target ", j, ":\n"))
        for (i in 1:n_gammas[[j]]) {
            final_rhat <- tail(cum_rhat_burnin_gammas[[j]][[i]], n = 1)
            converged_at <- burnin_convergence_time_gammas[[j]][[i]]
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

    cat("5. ESS Analysis for joint Betas and Gammas\n")
    cat("   - Median minimum ESS for sampled betas: ", round(median_min_ess_beta, 4), " ESS\n")
    cat("   - Median minimum ESS for sampled gammas: ", round(median_min_ess_gamma, 4), " ESS\n")
    cat("   - Median minimum ESS/s for sampled betas: ", round(median_min_ess_per_second_beta, 4), " ESS/s\n")
    cat("   - Median minimum ESS/s for sampled gammas: ", round(median_min_ess_per_second_gamma, 4), " ESS/s\n")
    cat("   - Median minimum ESS/iteration for sampled betas: ", round(median_min_ess_per_iteration_beta, 4), " ESS/iteration\n")
    cat("   - Median minimum ESS/iteration for sampled gammas: ", round(median_min_ess_per_iteration_gamma, 4), " ESS/iteration\n")

    cat("6. Summary\n")
    cat("   - ", n_burnin_beta_convergences, "/", n_beta_params, "(",
        to_percent(n_burnin_beta_convergences / n_beta_params), ") burn-in betas converged.\n")
    cat("   - ", n_burnin_gamma_convergences, "/", n_gamma_params, "(",
        to_percent(n_burnin_gamma_convergences / n_gamma_params), ") burn-in gammas converged.\n")
    cat("   - ", n_sampled_betas_convergences, "/", n_beta_params, "(",
        to_percent(n_sampled_betas_convergences / n_beta_params), ") sampled betas converged.\n")
    cat("   - ", n_sampled_gammas_convergences, "/", n_gamma_params, "(",
        to_percent(n_sampled_gammas_convergences / n_gamma_params), ") sampled gammas converged.\n")
    cat("   - Median convergence time for burn-in betas: ", median_convergence_time_beta, " samples\n")
    cat("   - Median convergence time for burn-in gammas: ", median_convergence_time_gammas, " samples\n")

}

print_diagnostic_report()