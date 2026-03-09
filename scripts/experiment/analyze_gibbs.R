# Load experiment data from "experiment_gibbs.R" and analyze it.

source("scripts/experiment/diagnostics.R")

# Load the cross-validation results for the Gibbs sampling experiment.
name <- "gibbs_50000_10"
cv_res = readRDS(paste0("results/cv_", name, ".rds"))
print("Loaded results")

# Run diagnostics on the loaded results.
cat("Running diagnostics...\n")
dia_res <- run_full_diagnostic(cv_res)

# Save diagnostic results.
saveRDS(dia_res, file = paste0("results/diagnostics_", name, ".rds"))
cat("Saved diagnostic results to ", paste0("results/diagnostics_", name, ".rds"), "\n")

# Plot and print diagnostics.
#plot_diagnostic_plots(dia_res)
print_diagnostic_report(dia_res)




