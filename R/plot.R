
library(ggplot2)

source("R/util.R")

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
            dev.copy2pdf(file=paste0(filename, "_", val, "_", gsub("\\.", "dot", colnames(values)[i]), ".pdf"), out.type = "pdf")
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
    gg <- ggplot(gamma_long, aes(x = value, color = group, fill = group, group = interaction(group, gamma))) +
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
    # Package single eval object as a list.
    if (class(evals) == "mspm_labeled_evaluation") {
        evals <- list(evals)
    }

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

        if (plotMedian) {
            medians <- tapply(dat$drawMeans, dat$eval, median)
            for (i in seq_along(medians)) {
                gg <- gg + geom_vline(xintercept = medians[i], color=scales::hue_pal()(length(medians))[i], linetype="dashed", size=0.5)
            }
        }
        if (plotMean) {
            means <- tapply(dat$drawMeans, dat$eval, mean)
            for (i in seq_along(means)) {
                gg <- gg + geom_vline(xintercept = means[i], color=scales::hue_pal()(length(means))[i], linetype="solid", size=0.5)
            }
        }

        plot(gg) # Display the plot.
    }
}

# Verify that the eval objects contain at least all the provided metrics.
#
# Arguments:
# evals  : A list of mspm_labeled_evaluation objects.
# metrics: A character vector of metric names to check for.
.validateMetricsForEvalObjects1 <- function(evals, metrics) {
    for (eval in evals) {
        for (metric in metrics) {
            if (!(metric %in% eval$metrics)) {
                stop(paste("Evaluation object is missing metric:", metric))
            }
        }
    }
}

# Verify that the all eval objects contain the same metrics.
#
# Arguments:
# evals  : A list of mspm_labeled_evaluation objects.
.validateMetricsForEvalObjects2 <- function(evals) {
    base_metrics <- evals[[1]]$metrics
    for (eval in evals) {
        if (!identical(sort(base_metrics), sort(eval$metrics))) {
            stop("Evaluation objects contain different metrics.")
        }
    }
}