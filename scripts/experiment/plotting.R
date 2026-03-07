library(ggplot2)

#' Create a Density Plot
#'
#' Generates a density plot for a numeric vector using ggplot2.
#'
#' @param x A numeric vector to plot
#' @param title Optional title for the plot. Default is NULL.
#' @param xlabel Label for the x-axis. Default is "Value".
#' @param ylabel Label for the y-axis. Default is "Density".
#' @param color Color of the density line. Default is "steelblue".
#' @param fill Color of the density fill. Default is "steelblue" with alpha = 0.3.
#' @param alpha Transparency of the fill. Default is 0.3.
#' @param size Line size for the density curve. Default is 1.
#'
#' @return A ggplot object that can be further customized or printed.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' 
#' # Simple example
#' x <- rnorm(1000, mean = 0, sd = 1)
#' p <- density_plot(x, title = "Distribution of Values")
#' print(p)
#' 
#' # Customize appearance
#' p <- density_plot(x, 
#'                   title = "Custom Density Plot",
#'                   xlabel = "Parameter Value",
#'                   color = "darkblue",
#'                   fill = "lightblue",
#'                   alpha = 0.5)
#' print(p)
#' }
#'
#' @export
density_plot <- function(x, 
                         title = NULL,
                         xlabel = "Value",
                         ylabel = "Density",
                         color = "steelblue",
                         fill = "steelblue",
                         alpha = 0.3,
                         size = 1) {
  
  # Ensure input is numeric
  x <- as.numeric(x)
  
  # Remove NA values
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    stop("Input vector contains no valid (non-NA) values")
  }
  
  # Create a data frame for ggplot
  df <- data.frame(value = x)
  
  # Calculate median
  med <- median(x)
  
  # Create the density plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = value)) +
    ggplot2::geom_density(color = color, 
                          fill = fill, 
                          alpha = alpha,
                          size = size) +
    ggplot2::geom_vline(xintercept = med, color = "black", linetype = "dashed", size = 0.5) +
    ggplot2::labs(x = xlabel,
                  y = ylabel,
                  title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"))

  # Add some horizontal spacing.
  xlim <- range(x)
  x_range <- xlim[2] - xlim[1]
  xlim <- c(xlim[1] - 0.6 * x_range, xlim[2] + 0.6 * x_range)
  p <- p + ggplot2::xlim(xlim)
  
  return(p)
}