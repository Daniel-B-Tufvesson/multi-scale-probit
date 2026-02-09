
library(ggplot2)

source("R/util.R")
source("R/internal.R")

# Plot all posterior beta distributions. This creates a separate plot for each beta 
# parameter.
# 
# Arguments:
# fit            : Fitted roprobit model object.
# bw             : Bandwidth for density estimation (default: 0.15).
# color_gradient: Color gradient for density fill based on sign certainty (default: c("red
#                  , "green3")).
# filename       : If provided, saves the plots to a PDF file with this name (default: NULL).
plot_posteriors_beta <- function(
    fit,
    ...,
    bw = 0.15,
    color_gradient = c("red3", "green3"),
    filename = NULL
) {
    # Create a color palette for the density fill
    pal <- colorRampPalette(color_gradient)(100)

    dat <- data.frame(fit$beta)
    for (i in 1:ncol(dat)) {
        # Calculate the proportion of samples greater than zero for the parameter
        posneg <- mean(dat[, i] > 0)

        # Determine the fill color and label based on the sign certainty
        if (posneg > 0.5) {
            val <- floor(200*(posneg-0.5))
            fillcolor <- pal[max(1, val)]
        }
        else if (posneg < 0.5) {
            val <- floor(200*(0.5-posneg))
            fillcolor <- pal[max(1, val)]
        }
        else {
            fillcolor <- pal[1]
        }

        name <- paste(colnames(dat)[i])

        # Create the density plot for the parameter
        gg <- ggplot(dat, aes(x=!!dat[, i])) +
            geom_density(bw=bw, fill=fillcolor, alpha=.6) + # Density plot with fill
            labs(x = name) + # X-axis label with certainty
            geom_vline(xintercept = 0) + # Vertical line at zero
            theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill="white", color="white"),
                axis.title.x = element_text(size = 30),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                legend.position = "none",
                axis.text.x = element_text(size=30))
        
        plot(gg) # Display the plot

        # If a filename is provided, save the plot to a PDF file
        if (!is.null(filename)) {
            dev.copy2pdf(file=paste0(filename, "_", val, "_", gsub("\\.", "dot", 
                         colnames(values)[i]), ".pdf"), out.type = "pdf")
        }
    }
}

# Plot all posterior gamma distributions in the same graph.
#
# Arguments:
# fit    : Fitted roprobit model object.
# colors : Vector of colors for different gamma groups (default: NULL).
# xlim   : Limits for x-axis (default: NULL).
# filename: If provided, saves the plot to a PDF file with this name (default: NULL).
plot_posterior_gammas <- function(
    fit,
    ...,
    colors = NULL,
    xlim = NULL,
    filename = NULL
) {
    #dat <- data.frame(fit$gammas)
    n_groups <- length(fit$gammas)

    # Create a color palette if colors are not provided
    if (is.null(colors)) {
        colors <- rainbow(n_groups)
    }

    # Reshape all gamma samples into a long-format data frame
    gamma_long <- data.frame()
    for (i in 1:n_groups) {
        gamma_samples <- data.frame(fit$gammas[[i]])
        n_gamma_i <- ncol(gamma_samples)
        for (j in 1:n_gamma_i) {
            temp <- data.frame(
                value = gamma_samples[, j],
                group = factor(paste0("Group ", i)),
                gamma = factor(paste0("Gamma ", j))
            )
            gamma_long <- rbind(gamma_long, temp)
        }
    }

    # Plot all densities in one ggplot call, grouped and colored by group and gamma
    gg <- ggplot(gamma_long, aes(x = value, color = group, fill = group, 
                 group = interaction(group, gamma))) +
        geom_density(alpha = 0.2) +
        xlab(expression(paste('x'^'T',beta))) + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill="white", color="white"),
              axis.title.x = element_text(size = 15),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              legend.title = element_blank(),
              axis.text.x = element_text(size=15))

    # Set x-axis limits if provided
    if (!is.null(xlim)) {
        gg <- gg + xlim(xlim)
    }

    plot(gg) # Display the plot

    if (!is.null(filename)) {
        dev.copy2pdf(file=paste0(filename, ".pdf"), out.type = "pdf")
    }
}

# Plot evaluation results for each draw. Creates a density plot of the evaluation metric
# values across all draws.
#
# Arguments:
# evals  : Either a single mspm_labeled_evaluation object or a list of such objects.
# metrics: The evaluation metric to plot (default: NULL, which plots all the metrics).
# plotMean: Whether to plot the mean of the metric distribution (default: FALSE).
# plotMedian: Whether to plot the median of the metric distribution (default: FALSE).
plot_eval_draws <- function(
    evals,
    ...,
    metrics = NULL,
    plotMean = FALSE,
    plotMedian = FALSE,
    xlim = c(0,1)
) {
    # Todo: add option for plotting each target instead of just the mean over them.

    # Package single eval object as a list.
    if (class(evals) == "mspm_labeled_evaluation") {
        evals <- list(evals)
    }

    # Validate metrics.
    if (is.null(metrics)) {
        .validateMetricsForEvalObjects2(evals)
        metrics <- evals[[1]]$metrics
    }
    else {
        .validateMetricsForEvalObjects1(evals, metrics)
    }

    # For each metric, collect drawMeans from all evals and plot together
    for (metric in metrics) {
        # Collect drawMeans for this metric from all evals
        dat_list <- lapply(seq_along(evals), function(i) {
            data.frame(drawMeans = evals[[i]]$drawMeans[[metric]],
                       eval = factor(paste0("Eval ", i)))
        })
        dat <- do.call(rbind, dat_list)

        # Plot the density of drawMeans for this metric
        gg <- ggplot(dat, aes(x=drawMeans, fill=eval, color=eval)) +
            geom_density(alpha=.6) +
            labs(x = paste0(metric, " score")) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill="white", color="white"),
                  axis.title.x = element_text(size = 15),
                  axis.title.y = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  legend.title = element_blank(),
                  legend.text = element_text(size=12),
                  axis.text.x = element_text(size=15))

        if (!is.null(xlim)) {
            gg <- gg + xlim(xlim)
        }

        # Plot median lines.
        if (plotMedian) {
            medians <- tapply(dat$drawMeans, dat$eval, median)
            for (i in seq_along(medians)) {
                gg <- gg + geom_vline(xintercept = medians[i], 
                                      color=scales::hue_pal()(length(medians))[i], 
                                      linetype="dashed", size=0.5)
            }
        }
        # Plot mean lines.
        if (plotMean) {
            means <- tapply(dat$drawMeans, dat$eval, mean)
            for (i in seq_along(means)) {
                gg <- gg + geom_vline(xintercept = means[i], 
                                      color=scales::hue_pal()(length(means))[i], 
                                      linetype="solid", size=0.5)
            }
        }

        plot(gg) # Display the plot.
    }
}

# Verify that the provided metrics vectors contain at least all the provided metrics.
#
# Arguments:
# metricsList  : A list of character vectors for the metrics.
# requiredMetrics: A character vector of metric names to check for.
.validateMetrics1 <- function(metricsList, requiredMetrics) {
    for (metrics in metricsList) {
        for (requiredMetric in metrics) {
            if (!(requiredMetric %in% metrics)) {
                stop(paste("Evaluation object is missing metric:", requiredMetric))
            }
        }
    }
}

# Verify that the all character vectors contain the same metrics.
#
# Arguments:
# evals  : A list of character vectors containing the metrics.
.validateMetrics2 <- function(metricsList) {
    base_metrics <- sort(metricsList[[1]])
    for (metrics in metricsList) {
        if (!identical(base_metrics, sort(metrics))) {
            stop("Evaluation objects contain different metrics.")
        }
    }
}


# Plot the distribution of differences between two evaluations over multiple draws. 
#
# Arguments:
# eval1  : The first mspm_labeled_evaluation object.
# eval2  : The second mspm_labeled_evaluation object.
# metrics: The evaluation metric to plot (default: NULL, which plots all the metrics).
# label1 : Label for the first evaluation in the legend (default: "1").
# label2 : Label for the second evaluation in the legend (default: "2").
# plotData: Which data to plot the differences from. Options are "drawMeans" (default), 
#             "allDraws", or "metricMeans".
plot_eval_draws_diff <- function(
    eval1,
    eval2,
    ...,
    metrics = NULL,
    label1 = "1",
    label2 = "2",
    plotData = "drawMeans",
    addXLabel = TRUE
) {
    # Todo: add option for plotting each target instead of just the mean over them.

    # Validate data.
    supportedData = c("drawMeans", "allDraws", "metricMeans")
    if (!(plotData %in% supportedData)) {
        stop(paste("plotData must be one of:", paste(supportedData, collapse=", ")))
    }

    # Validate eval objects.
    if (class(eval1) != "mspm_labeled_evaluation" || 
        class(eval2) != "mspm_labeled_evaluation") {
        stop("Both eval1 and eval2 must be mspm_labeled_evaluation objects.")
    }

    # Validate metrics.
    if (is.null(metrics)) {
        .validateMetrics2(list(evalMetrics(eval1), evalMetrics(eval2)))
        metrics <- evalMetrics(eval1)
    }
    else {
        .validateMetrics1(list(evalMetrics(eval1), evalMetrics(eval2)), metrics)
    }

    # Plot drawMeans.
    if (plotData == "drawMeans") {
        # For each metric, compute differences and plot
        for (metric in metrics) {
            xlabel <- ifelse(addXLabel, paste0("Difference in ", metric), "")
            .plot_dist_diff(
                eval1$drawMeans[[metric]], 
                eval2$drawMeans[[metric]],
                label1,
                label2,
                xlabel
            )
        }
    }
    # Plot metricMeans.
    else if (plotData == "metricMeans") {
        # For each metric, compute differences and plot
        ntargets = ntargets(eval1)
        for (target in 1:ntargets) {
            xlabel <- ifelse(addXLabel, paste0("Difference in mean over metrics for ", target), "")
            .plot_dist_diff(
                eval1$metricMeans[, target], 
                eval2$metricMeans[, target],
                label1,
                label2,
                xlabel
            )
        }

        # Plot mean as well if more than 1 target.
        if (ntargets > 1) {
            xlabel <- ifelse(addXLabel, "Difference in mean over metrics (harmonic mean)", "")
            .plot_dist_diff(
                eval1$metricMeans[, ntargets + 1], 
                eval2$metricMeans[, ntargets + 1],
                label1,
                label2,
                xlabel
            )
        }
    }
    # Plot allDraws.
    else if (plotData == "allDraws") {
        # For each metric, compute differences and plot
        for (metric in metrics) {
            # Get the results matrices for this metric
            res1 <- eval1$drawResults[[metric]]
            res2 <- eval2$drawResults[[metric]]

            # Plot for each target.
            ntargets <- ntargets(eval1)
            for (target in 1:ntargets) {
                xlabel <- ifelse(addXLabel, 
                                 paste0("Difference in ", metric, " for target ", target), "")
                .plot_dist_diff(
                    res1[, target], 
                    res2[, target],
                    label1,
                    label2,
                    xlabel
                )
            }

            # Plot mean as well if more than 1 target.
            if (ntargets > 1) {
                xlabel <- ifelse(addXLabel, 
                                 paste0("Difference in ", metric, " (harmonic mean)"), "")
                .plot_dist_diff(
                    res1[, ntargets + 1], 
                    res2[, ntargets + 1],
                    label1,
                    label2,
                    xlabel
                )
            }
        }
    }
    else {
        stop("Unsupported dataToPlot option.")
    }

}


.plot_dist_diff <- function(
    values1,
    values2,
    label1,
    label2,
    xlabel = NULL
) {
    diffs <- values1 - values2
    dat <- data.frame(differences = diffs)

    # Plot the density of differences for this metric
    # Compute density manually to split left/right of zero
    dens <- density(dat$differences, na.rm = TRUE)
    dens_df <- data.frame(x = dens$x, y = dens$y)
    dens_df$side <- factor(ifelse(dens_df$x < 0, "left", "right"), levels = c("left", "right"))

    # Calculate area under the density curve for each side
    total_area <- sum(dens_df$y)
    left_area <- sum(dens_df$y[dens_df$side == "left"]) / total_area
    right_area <- sum(dens_df$y[dens_df$side == "right"]) / total_area

    # Format as percentages
    left_pct <- sprintf("%.1f%%", 100 * left_area)
    right_pct <- sprintf("%.1f%%", 100 * right_area)

    # Add labels for legend with percentages
    side_labels <- c(
        left = paste0(label1, " (", left_pct, ")"),
        right = paste0(label2, " (", right_pct, ")")
    )

    # Create the plot
    gg <- ggplot(dens_df, aes(x = x, y = y, fill = side)) +
        geom_area(alpha = .6) +
        scale_fill_manual(
            name = "",
            values = c(left = "#D55E00", right = "#009E73"),
            labels = side_labels
        ) +
        geom_line(color = "black", inherit.aes = FALSE, data = dens_df, aes(x = x, y = y)) +
        geom_vline(xintercept = 0, linetype = "solid", color = "black") +
        theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill="white", color="white"),
                axis.title.x = element_text(size = 15),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.text.x = element_text(size=15),
                legend.title = element_blank(),
                legend.text = element_text(size=12),
                legend.position = c(0.02, 0.98),
                legend.justification = c(0, 1),
                legend.background = element_blank(),
                legend.key.size = unit(1.5, 'lines'))
    
    gg <- gg + labs(x = ifelse(is.null(xlabel), "", xlabel))

    plot(gg) # Display the plot
}


# Plot CV distributions --------------------------------------------------------------------------


#' Plot the distribution of differences between two cross-validation results. 
#'
#' @param cv_res1 The first mspm_cv_result object.
#' @param cv_res2 The second mspm_cv_result object.
#' @param metrics The evaluation metric to plot (default: NULL, which plots all the metrics).
#' @param plotHarmonicMeans Whether to plot the harmonic means when multiple metrics are present (default: TRUE).
#' @param label1 Label for the first CV result in the legend (default: "CV 1").
#' @param label2 Label for the second CV result in the legend (default: "CV 2").
#' @param plotData Which data to plot the differences from. Options are "drawMeans" (default), or "allDraws".
#' @param addXLabel Whether to add an x-axis label indicating the metric (default: TRUE).
plot_cv_diff <- function(
    cv_res1,
    cv_res2,
    ...,
    metrics = NULL,
    plotHarmonicMeans = TRUE,
    label1 = "CV 1",
    label2 = "CV 2",
    plotData = "means",
    addXLabel = TRUE
) {
    # Validate plot option.
    supportedData = c("means", "allDraws")
    if (!(plotData %in% supportedData)) {
        stop(paste("plotData must be one of:", paste(supportedData, collapse=", ")))
    }

    # Validate eval objects.
    if (class(cv_res1) != "mspm_cv_result" || 
        class(cv_res2) != "mspm_cv_result") {
        stop("Both cv_res1 and cv_res2 must be mspm_cv_result objects.")
    }

    # Validate metrics.
    if (is.null(metrics)) {
        .validateMetrics2(list(evalMetrics(cv_res1), evalMetrics(cv_res2)))
        metrics <- evalMetrics(cv_res1)
    }
    else {
        .validateMetrics1(list(evalMetrics(cv_res1), evalMetrics(cv_res2)), metrics)
    }

    # Plot means.
    if (plotData == "means") {

        # For each metric, compute differences and plot
        for (metric in metrics) {
            for (target in 1:ntargets(cv_res1)) {
                xlabel <- ifelse(addXLabel, paste0("Difference in mean ", metric, " for target ", target))
                .plot_dist_diff(
                    cvMeans(cv_res1)[[target]][, metric],
                    cvMeans(cv_res2)[[target]][, metric],
                    label1,
                    label2,
                    xlabel
                )
            }
        }

        # Plot harmonic mean if needed.
        if (plotHarmonicMeans && length(metrics) > 1) {
            for (target in 1:ntargets(cv_res1)) {
                xlabel <- ifelse(addXLabel, paste0("Difference in harmonic mean for target ", target))
                .plot_dist_diff(
                    cvMeans(cv_res1)[[target]][, "HarmonicMean"],
                    cvMeans(cv_res2)[[target]][, "HarmonicMean"],
                    label1,
                    label2,
                    xlabel
                )
            }
        }
    }
    # Todo: implment this.
    # Plot allDraws.
    # else if (plotData == "allDraws") {
    #     # For each metric, compute differences and plot
    #     for (metric in metrics) {
    #         # Get the results matrices for this metric
    #         res1 <- cvDrawResults(cv_res1)[[metric]]
    #         res2 <- cvDrawResults(cv_res2)[[metric]]

    #         # Plot for each target.
    #         ntargets <- ncol(res1)
    #         for (target in 1:ntargets) {
    #             xlabel <- ifelse(addXLabel, 
    #                              paste0("Difference in ", metric, " for target ", target), "")
    #             .plot_dist_diff(
    #                 res1[, target], 
    #                 res2[, target],
    #                 label1,
    #                 label2,
    #                 xlabel
    #             )
    #         }
    #     }
    # }
    # else {
    #     stop("Unsupported dataToPlot option.")
    # }
}


# Plot MCMC chains -------------------------------------------------------------------------------


#' Plot the MCMC chains for the beta of the fit.
#'
#' @param fit An object of class 'mspm' containing the fitted model.
#' @param title An optional title for the plot (default: "Beta").
plot_beta_chains <- function(fit, title = "Beta") {
    chains <- as.matrix(beta(fit))
    n_iter <- nrow(chains)
    n_param <- ncol(chains)
    param_labels <- paste0("beta_", seq_len(n_param))
    math_labels <- lapply(seq_len(n_param), function(i) bquote(beta[.(i)]))
    chain_df <- data.frame(
        Iteration = rep(1:n_iter, times = n_param),
        Value = as.vector(chains),
        Parameter = factor(rep(param_labels, each = n_iter), levels = param_labels)
    )

    # Plot the chains.
    gg <- ggplot(chain_df, aes(x = Iteration, y = Value, color = Parameter)) +
        geom_line() +
        labs(x = "Iteration", y = "", color = "Parameter") +
        scale_color_manual(
            values = scales::hue_pal()(n_param),
            labels = math_labels
        ) +
        theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill="white", color="white"),
                axis.title.x = element_text(size = 15),
                axis.title.y = element_text(size = 15),
                legend.title = element_blank(),
                legend.text = element_text(size=12),
                axis.text.x = element_text(size=15),
                axis.text.y = element_text(size=15))

    # Plot x-axis.
    y_min <- min(chain_df$Value, na.rm = TRUE)
    y_max <- max(chain_df$Value, na.rm = TRUE)
    if (y_min <= 0 && y_max >= 0) {
        gg <- gg + geom_hline(yintercept = 0, color = "black", linetype = "dashed")
    }

    # Add title.
    if(!is.null(title)) {
        gg <- gg + ggtitle(title) + theme(plot.title = element_text(size=20))
    }

    plot(gg)
}

#' Plot the MCMC chains for all the gammas of the fit. Each gamma is plotted in a separate graph.
#'
#' @param fit An object of class 'mspm' containing the fitted model.
plot_gamma_chains <- function(
    fit,
    ...,
    seperateGraphs = FALSE
) {
    gammas <- gammas(fit)

    if (seperateGraphs) {
        for (i in seq_along(gammas)) {
            .plot_gamma_chains(gammas[[i]], paste0("Gammas Scale ", i))
        }
    }
    else {
        .plot_gamma_chains_grouped_by_scale(gammas, title = "All gammas")
    }
}

#' Plot the MCMC chains for the gamma of the fit.
#'
#' @param gamma An mcmc object containing the gamma samples for one scale.
#' @param title An optional title for the plot (default: NULL).
.plot_gamma_chains <- function(gamma, title = NULL) {
    chains <- as.matrix(gamma)
    n_iter <- nrow(chains)
    n_param <- ncol(chains)
    param_labels <- paste0("gamma_", seq_len(n_param))
    math_labels <- lapply(seq_len(n_param), function(i) bquote(gamma[.(i)]))
    chain_df <- data.frame(
        Iteration = rep(1:n_iter, times = n_param),
        Value = as.vector(chains),
        Parameter = factor(rep(param_labels, each = n_iter), levels = param_labels)
    )

    gg <- ggplot(chain_df, aes(x = Iteration, y = Value, color = Parameter)) +
        geom_line() +
        labs(x = "Iteration", y = "", color = "Parameter") +
        scale_color_manual(
            values = scales::hue_pal()(n_param),
            labels = math_labels
        ) +
        theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill="white", color="white"),
                axis.title.x = element_text(size = 15),
                axis.title.y = element_text(size = 15),
                legend.title = element_blank(),
                legend.text = element_text(size=12),
                axis.text.x = element_text(size=15),
                axis.text.y = element_text(size=15))
    
    # Plot x-axis.
    y_min <- min(chain_df$Value, na.rm = TRUE)
    y_max <- max(chain_df$Value, na.rm = TRUE)
    if (y_min <= 0 && y_max >= 0) {
        gg <- gg + geom_hline(yintercept = 0, color = "black", linetype = "dashed")
    }

    # Add title.
    if(!is.null(title)) {
        gg <- gg + ggtitle(title) + theme(plot.title = element_text(size=20))
    }

    plot(gg)
}

#' Plot the MCMC chains for all gammas grouped by scale in the same graph.
#'
#' @param gammas A list of mcmc objects containing the gamma samples for each scale.
#' @param title An optional title for the plot (default: NULL).
.plot_gamma_chains_grouped_by_scale <- function(gammas, title = NULL) {
    # Combine all gamma chains into one data frame
    chain_list <- list()
    scale_labels <- list()
    for (scale_idx in seq_along(gammas)) {
        chains <- as.matrix(gammas[[scale_idx]])
        n_iter <- nrow(chains)
        n_param <- ncol(chains)
        for (param_idx in seq_len(n_param)) {
            chain_list[[length(chain_list) + 1]] <- data.frame(
                Iteration = 1:n_iter,
                Value = chains[, param_idx],
                Scale = scale_idx,
                Gamma = param_idx
            )
            scale_labels[[length(scale_labels) + 1]] <- bquote(gamma[.(param_idx)]^{(.(scale_idx))})
        }
    }
    chain_df <- do.call(rbind, chain_list)
    chain_df$Label <- factor(
        paste0("gamma_", chain_df$Gamma, "_group_", chain_df$Scale),
        levels = unique(paste0("gamma_", chain_df$Gamma, "_group_", chain_df$Scale))
    )

    # Plot the chains.
    gg <- ggplot(chain_df, aes(x = Iteration, y = Value, color = Label)) +
        geom_line() +
        labs(x = "Iteration", y = "", color = "Parameter") +
        scale_color_manual(
            values = scales::hue_pal()(length(scale_labels)),
            labels = scale_labels
        ) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill="white", color="white"),
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15),
              legend.title = element_blank(),
              legend.text = element_text(size=12),
              axis.text.x = element_text(size=15),
              axis.text.y = element_text(size=15))

    # Plot x-axis.
    y_min <- min(chain_df$Value, na.rm = TRUE)
    y_max <- max(chain_df$Value, na.rm = TRUE)
    if (y_min <= 0 && y_max >= 0) {
        gg <- gg + geom_hline(yintercept = 0, color = "black", linetype = "dashed")
    }

    # Plot title.
    if(!is.null(title)) {
        gg <- gg + ggtitle(title) + theme(plot.title = element_text(size=20))
    }

    plot(gg)
}