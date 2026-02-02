
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
# eval  : The mspm_labeled_evaluation object containing evaluation results.
# metrics: The evaluation metric to plot (default: NULL, which plots all the metrics).
plot_eval_draws <- function(
    eval,
    ...,
    metrics = NULL,
    plotMean = TRUE,
    plotMedian = TRUE
) {
    if (is.null(metrics)) {
        metrics <- eval$metrics
    }

    # Draw the mean over the targets.
    for (metric in metrics) {
        drawMeans = eval$drawMeans[[metric]]
        dat <- data.frame(drawMeans)

        gg <- ggplot(dat, aes(x=drawMeans)) +
            geom_density(fill="blue", alpha=.6) + # Density plot with fill
            labs(x = paste0(metric, " score")) + # X-axis label
            theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill="white", color="white"),
                axis.title.x = element_text(size = 30),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                legend.position = "none",
                axis.text.x = element_text(size=30))

        if (plotMedian) {
            median <- median(drawMeans)
            gg <- gg + geom_vline(xintercept = median, color="black", linetype="dashed", size=0.5)
        }
        if (plotMean) {
            mean <- mean(drawMeans)
            gg <- gg + geom_vline(xintercept = mean, color="black", linetype="solid", size=0.5)
        }
        
        plot(gg) # Display the plot.
    }
}