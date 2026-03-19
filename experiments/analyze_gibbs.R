# Load samples from the gibbs experiment.

source("experiments/diagnostics.R")

# Load the samples from the gibbs runs.
name <- "gibbs-samples-big"
all_runs <- readRDS(paste0("experiments/results/", name, ".rds"))
print("Loaded results")

# Run diagnostics on the loaded results.
cat("Running diagnostics...\n")
dia_res <- run_full_diagnostic(all_runs)

# # Save diagnostic results.
# saveRDS(dia_res, file = paste0("experiments/results/diagnostics-", name, ".rds"))
# cat("Saved diagnostic results to ", paste0("results/diagnostics-", name, ".rds"), "\n")

# # Plot and print diagnostics.
plot_diagnostic_plots(dia_res)
print_diagnostic_report(dia_res)

plot_cum_rhat(dia_res$cum_rhat_burnin_beta[[1]], 100)
plot_cum_rhat(dia_res$cum_rhat_burnin_beta[[2]], 100)
plot_cum_rhat(dia_res$cum_rhat_burnin_beta[[3]], 100)


