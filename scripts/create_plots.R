# MOVE FILE TO scripts/create_plots.R

library(ggplot2)
#install.packages("tikzDevice")
library(grDevices)
library(reshape2)
library(plyr)
library(coda)

source("R/util.R")

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(1,2,4,5,6,7,8,3)]

# The palette with black:
cbbPalette1 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbbPalette2 <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
                "#CC79A7", "#F0E442")

plot_stats_gg <- function(separate, 
                        combined,
                        xlim=c(0, 1), 
                        ylim=c(0,4),
                        bw = "nrd0", 
                        set_names = NULL,
                        filename = NULL, 
                        model.type = "Binary",
                        custom.palette = NULL,
                        monochrome = FALSE) {
  # Set colors
  if (is.null(custom.palette)) {
    active.palette <- if (model.type == "Binary") cbbPalette2 else cbbPalette1
  } else {
    active.palette <- custom.palette
  }
  # Prepare data
  dat <- data.frame(combined, separate)
  colnames(dat) <- c("Multi-Scale Probit", paste(model.type, "Probit"))
  dat <- melt(dat)
  colnames(dat) <- c("model", "draw")
  means <- ddply(dat[dat$draw >= 0 & dat$draw <= 1, ], "model", summarise, grp.mean=mean(draw, na.rm = TRUE))
  # Plot
  if (monochrome) {
    gg <- ggplot(dat, aes(x=draw, linetype=model)) +
      scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
      geom_vline(data=means, aes(xintercept=grp.mean, linetype=model))
  } else {
    gg <- ggplot(dat, aes(x=draw, fill=model)) +
      geom_vline(data=means, aes(xintercept=grp.mean, color=model)) +
      scale_fill_manual(values=active.palette) +
      scale_color_manual(values=active.palette)
  }
  # Finalize plot
  gg <- gg + theme(panel.background = element_rect(fill="white", color="white")) +
    coord_cartesian(xlim=xlim) +
    geom_density(bw=bw, alpha=.5) +
    xlim(xlim) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size=30),
          legend.title = element_blank(),
          legend.text=element_text(size=30),
          legend.justification=c(0, 1),
          legend.position.inside=c(0, 1),
          legend.background = element_blank(),
          legend.key.size = unit(2, 'lines'))
  plot(gg)
  # Uncomment to save plot as PDF
  # if (!is.null(filename)) {
  #   dev.copy2pdf(file=filename, out.type = "pdf")
  # }
}

# combined_sim = readRDS("eval_sim_48_all_2018-09-18.RDS")
# separate_sim1 = readRDS("eval_sim_48_1_2018-09-18.RDS")
# separate_sim2 = readRDS("eval_sim_48_2_2018-09-18.RDS")
# separate_sim3 = readRDS("eval_sim_48_3_2018-09-18.RDS")
# separate_sim4 = readRDS("eval_sim_48_4_2018-09-18.RDS")

# combined_sim = readRDS("eval_sim_48_3sets_all_2018-10-02.RDS")
# separate_sim1 = readRDS("eval_sim_48_3sets_1_2018-10-02.RDS")
# separate_sim2 = readRDS("eval_sim_48_3sets_2_2018-10-02.RDS")
# separate_sim3 = readRDS("eval_sim_48_3sets_3_2018-10-02.RDS")

# combined_sim = readRDS("eval_sim_8_all_2018-09-14.RDS")
# separate_sim1 = readRDS("eval_sim_8_1_2018-09-14.RDS")
# separate_sim2 = readRDS("eval_sim_8_2_2018-09-14.RDS")
# separate_sim3 = readRDS("eval_sim_8_3_2018-09-14.RDS")
# separate_sim4 = readRDS("eval_sim_8_4_2018-09-14.RDS")




print_totals <- function(model_totals, name, sleeptime = 1.0, monochrome = FALSE) {
  if (monochrome) {
    folder <- "images_mc"
  } else {
    folder <- "images"
  }
  for (i in c(1,2,4,6)) {
    print(metric_names[i])
    plot_stats_gg(model_totals[[1]][, i],
                  model_totals[[4]][, i],
                  bw = 0.05,
                  filename = paste0(folder,"/",metric_names[i],"_",name[1],".pdf"),
                  monochrome = monochrome)
    Sys.sleep(sleeptime)
    plot_stats_gg(model_totals[[2]][, i],
                  model_totals[[5]][, i],
                  bw = 0.05,
                  filename = paste0(folder,"/",metric_names[i],"_",name[2],".pdf"),
                  model.type = "Ordered",
                  monochrome = monochrome)
    Sys.sleep(sleeptime)
    plot_stats_gg(model_totals[[3]][, i],
                  model_totals[[6]][, i],
                  bw = 0.05,
                  filename = paste0(folder,"/",metric_names[i],"_",name[3],".pdf"),
                  model.type = "Ordered",
                  monochrome = monochrome)
    Sys.sleep(sleeptime)
  }
}


metric_names <- strsplit("F1_train,
                         F1_test,
                         kendall_train,
                         kendall_test,
                         harm_train,
                         harm_test", ",
                         ")[[1]]
simnames <- c("binary", "ordered1", "ordered2")
# corpus_names <- c("LL_binary", "Legimus_ordered1", "Hegas_ordered2")
trial60 <- readRDS("sim_60_dataset_all500_finished2_reorganized.RDS")
# trial600noshrink <- readRDS("sim_600_dataset_all500_noshrinkage_finished2_reorganized.RDS")
# trial_text <- readRDS("text_2thirds_dataset_all500_reorganized.RDS")
# trial_text_small <- readRDS("text_1thirds_dataset_all500_reorganized.RDS")


# Plot the synthetic trials.
print_totals(trial60$metric_totals, name = paste0("60_sim_", simnames), monochrome = TRUE)
print_totals(trial60$metric_totals, name = paste0("60_sim_", simnames), monochrome = FALSE)

# Commented out: Plot the other datasets. 
# print_totals(trial600noshrink$metric_totals, name = paste0("600_sim_", simnames), monochrome = TRUE)
# print_totals(trial600noshrink$metric_totals, name = paste0("600_sim_", simnames), monochrome = FALSE)
# 
# print_totals(trial_text$metric_totals, name = paste0("text_", corpus_names), monochrome = TRUE)
# print_totals(trial_text$metric_totals, name = paste0("text_", corpus_names), monochrome = FALSE)
# 
# print_totals(trial_text_small$metric_totals, name = paste0("text_1third_", corpus_names), monochrome = TRUE)
# print_totals(trial_text_small$metric_totals, name = paste0("text_1third_", corpus_names), monochrome = FALSE)


# Already commented out.
# combined_text <- readRDS("eval_1third_48_all_2018-09-19.RDS")
# separate_text1 <- readRDS("eval_1third_48_1_2018-09-19.RDS")
# separate_text2 <- readRDS("eval_1third_48_2_2018-09-19.RDS")
# separate_text3 <- readRDS("eval_1third_48_3_2018-09-19.RDS")

# str(combined_text)


# k = 455
model_means <- readRDS(paste0("sim_60_dataset_trials_means", k, ".RDS"))
model_means <- rbind(model_means, readRDS(paste0("sim_60_dataset_trials_means", k, ".RDS")))
model_means <- lapply(model_means, function(mat) {
   return(mat[1:400, ])
 })
 saveRDS(model_means, paste0("sim_60_dataset_trials_means", k ,".RDS"))
  model_means <- readRDS("text_2thirds_dataset_trials1000_finished2.RDS")

 for (k in c(62, 122, 227)) {
   tmp <- readRDS(paste0("sim_60_dataset_trials_means", k, ".RDS"))
   for (i in 1:6) {
     model_means[[i]] <- rbind(model_means[[i]], tmp[[i]])
   }
 }
 saveRDS(model_means, paste0("sim_60_dataset_trials_means", 455 ,".RDS"))

 model_means <- readRDS("text_2thirds_dataset_trials_means500_finished.RDS")

 metric_names <- strsplit("F1_train,
                          F1_test,
                          jaccard_train,
                          jaccard_test,
                          youden_train,
                          youden_test,
                          clharm_train,
                          clharm_test,
                          kendall_train,
                          kendall_test,
                          spearman_train,
                          spearman_test,
                          rharm_train,
                          rharm_test,
                          harm_train,
                          harm_test", ",
                          ")[[1]]

 print_means <- function(model_totals, name, sleeptime = 1.0, monochrome = FALSE) {
   if (monochrome) {
     folder <- "images_mc"
   } else {
     folder <- "images"
   }
   for (i in c(1,2,4,6)) {
     print(metric_names[i])
     plot_stats_gg(model_totals[[1]][, i],
                   model_totals[[4]][, i],
                   bw = 0.05,
                   filename = paste0(folder,"/",metric_names[i],"_means_",name[1],".pdf"),
                 custom.palette = c("#000000","#9ad0f3"),
                 monochrome = monochrome)
     Sys.sleep(sleeptime)
     plot_stats_gg(model_totals[[2]][, i],
                   model_totals[[5]][, i],
                   bw = 0.05,
                   filename = paste0(folder,"/",metric_names[i],"_means_",name[2],".pdf"),
                   model.type = "Ordered",
                 custom.palette = c("#000000", "#D55E00"),
                 monochrome = monochrome)
     Sys.sleep(sleeptime)
     plot_stats_gg(model_totals[[3]][, i],
                   model_totals[[6]][, i],
                   bw = 0.05,
                   filename = paste0(folder,"/",metric_names[i],"_means_",name[3],".pdf"),
                   model.type = "Ordered",
                 custom.palette = c("#000000", "#D55E00"),
                 monochrome = monochrome)
     Sys.sleep(sleeptime)
   }
 }

 print_means(trial60$metric_means, name = paste0("60_sim_", simnames), monochrome = TRUE)
 print_means(trial60$metric_means, name = paste0("60_sim_", simnames), monochrome = FALSE)

 print_means(trial600noshrink$metric_means, name = paste0("600_sim_", simnames), monochrome = TRUE)
 print_means(trial600noshrink$metric_means, name = paste0("600_sim_", simnames), monochrome = FALSE)

 print_means(trial_text$metric_means, name = paste0("text_", corpus_names), monochrome = TRUE)
 print_means(trial_text$metric_means, name = paste0("text_", corpus_names), monochrome = FALSE)
























# Compute difference between two sets of values.
#
# Arguments:
# set1: First set of values.
# set2: Second set of values.
difference <- function(set1, set2) {
  d <- set1 - set2
  xmax <- max(abs(d))
  return(list(d,
              mean(d < 0),
              mean(d > 0),
              cutpoint = 0,
              xlim = c(-xmax, xmax)))
}

# Compute ratio between two sets of values.
# 
# Arguments:
# set1: First set of values (numerator).
# set2: Second set of values (denominator).
ratio <- function(set1, set2) {
  r <- set1 / set2
  xmax <- max(abs(r - 1))
  return(list(r,
              mean(r < 1),
              mean(r > 1),
              cutpoint = 1,
              xlim = c(-xmax, xmax) + 1))
}

# Plot the distribution of the difference or ratio between two sets of values.
#
# Arguments:
# separate: Values from the single scale model(s).
# combined: Values from the multi-scale model.
# bw: Bandwidth for the density estimate (default: NULL, which uses the default R method).
# filename: Optional filename to save the plot as a PDF (default: NULL, which does not save).
# model.type: Type of the single scale model ("Binary" or "Ordered", default: "Binary").
# compare: Function to compute the comparison (default: difference).
# custom.palette: Optional custom color palette (default: NULL, which uses default colors).
# monochrome: Whether to use monochrome colors (default: FALSE).
plot_comparison_gg <- function(separate, combined,
                               bw = NULL,
                               filename = NULL,
                               model.type = "Binary",
                               compare = difference,
                               custom.palette = NULL,
                               monochrome = FALSE) {
  # Set color palette based on model type and options
  if (is.null(custom.palette)) {
    if (model.type == "Binary") {
      active.palette <- cbbPalette2[c(2, 1)]
    } else {
      active.palette <- cbbPalette1[c(2, 1)]
    }
  } else {
    active.palette <- custom.palette
  }
  # Override palette for monochrome plots
  if (monochrome) {
    active.palette <- c("#CCCCCC", "#000000")
  }

  # Prepare and clean input data
  df <- data.frame(separate, combined)
  df <- df[complete.cases(df), ]
  separate <- df$separate
  combined <- df$combined

  # Compute comparison statistics (difference or ratio)
  comparison <- compare(combined, separate)

  # Estimate density for the comparison metric
  if (is.null(bw)) {
    dat <- with(density(comparison[[1]]), data.frame(x, y))
  } else {
    dat <- with(density(comparison[[1]], bw = bw), data.frame(x, y))
  }

  # Adjust palette if all x values are above cutpoint
  if (all(dat$x > comparison[[4]])) {
    active.palette <- active.palette[c(2, 1)]
  }

  # Annotate which model is better for each region
  dat <- cbind(dat, Better = factor(as.numeric(dat$x>comparison[[4]]), levels = c(0, 1),
                                    labels = c(paste0(model.type, " (", round(100*comparison[[2]], 2), " %)"),
                                               paste0("Multi-Scale (", round(100*comparison[[3]], 2), " %)"))))

  # Create base ggplot object with density line and cutpoint
  gg <- ggplot(data = dat, mapping = aes(x = x, y = y)) +
    geom_line() +
    geom_vline(data = data.frame(cutpoints = c(comparison[[4]])), aes(xintercept=cutpoints)) + 
    ylim(c(0, max(dat$y)))

  # Add colored area under the density curve to highlight regions
  gg <- gg +
    geom_area(mapping = aes(fill = Better, col = Better), alpha = 0.3)

  # Set plot theme and axis formatting
  gg <- gg + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(fill="white", color="white"),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.text.x = element_text(size=30))

  # Configure legend appearance
  gg <- gg  + theme(legend.title = element_blank(),
                    legend.text=element_text(size=30),
                    legend.justification=c(0,1), 
                    legend.position.inside=c(0, 1),
                    legend.background = element_blank(),
                    legend.key.size = unit(2, 'lines'))

  # Set color and x-axis limits
  gg <- gg + scale_fill_manual(values=active.palette)
  gg <- gg + scale_color_manual(values=active.palette)
  gg <- gg + scale_x_continuous(limits = comparison$xlim)

  # Display the plot
  plot(gg)

  # Optionally save plot as PDF
  if (!is.null(filename)) {
    dev.copy2pdf(file=filename, out.type = "pdf")
  }
}


print_means_diff <- function(model_means, name, sleeptime = 1.0, monochrome = FALSE) {
    if (monochrome) {
      folder <- "images_mc"
    } else {
      folder <- "images"
    }
  for (i in c(1, 2, 4, 6)) {
    print(metric_names[i])
    plot_comparison_gg(model_means[[1]][, i],
                       model_means[[4]][, i],
                       bw = NULL, 
                       filename = paste0(folder,"/",metric_names[i],"_means_diff_",name[1],".pdf"),
                       compare = difference, # difference ?
                       monochrome = monochrome)
    Sys.sleep(sleeptime)
    plot_comparison_gg(model_means[[2]][, i],
                       model_means[[5]][, i],
                       bw = NULL, 
                       filename = paste0(folder,"/",metric_names[i],"_means_diff_",name[2],".pdf"),
                       model.type = "Ordered",
                       compare = difference,
                       monochrome = monochrome)
    
    Sys.sleep(sleeptime)
    plot_comparison_gg(model_means[[3]][, i],
                       model_means[[6]][, i],
                       bw = NULL, 
                       filename = paste0(folder,"/",metric_names[i],"_means_diff_",name[3],".pdf"),
                       model.type = "Ordered",
                       compare = difference,
                       monochrome = monochrome)
    Sys.sleep(sleeptime)
  }  
}

print_means_diff(trial60$metric_means, name = paste0("60_sim_", simnames), monochrome = TRUE)
print_means_diff(trial60$metric_means, name = paste0("60_sim_", simnames), monochrome = FALSE)

print_means_diff(trial600noshrink$metric_means, name = paste0("600_sim_", simnames), monochrome = TRUE)
print_means_diff(trial600noshrink$metric_means, name = paste0("600_sim_", simnames), monochrome = FALSE)

print(trial_text$metric_means[[1]][, 6])

print_means_diff(trial_text$metric_means, name = paste0("text_", corpus_names), monochrome = TRUE)
print_means_diff(trial_text$metric_means, name = paste0("text_", corpus_names), monochrome = FALSE)

print_means_diff(trial_text_small$metric_means, name = paste0("text_1third", corpus_names), monochrome = TRUE)
print_means_diff(trial_text_small$metric_means, name = paste0("text_1third", corpus_names), monochrome = FALSE)

























plot_rmse_gg <- function(separate, combined, known,
                         xlim=NULL, ylim=c(0,4),
                         bw = "nrd0", set_names = NULL,
                         filename = NULL, model.type = "Binary",
                         custom.palette = NULL,
                         monochrome = FALSE) {
  # Set colors
  if (is.null(custom.palette)) {
    if (model.type == "Binary") {
      active.palette <- cbbPalette2
    } else {
      active.palette <- cbbPalette1
    }
  } else {
    active.palette <- custom.palette
  }
  separate <- cpp_rmse_dist(separate, known, 500)
  combined <- cpp_rmse_dist(combined, known, 500)
  if (is.null(xlim)) {
    xlim <- c(0, max(max(separate), max(combined)))
  }
  dat <- data.frame(combined, separate)
  colnames(dat) <- c("Multi-Scale Probit", paste(model.type, "Probit"))
  dat <- melt(dat)
  colnames(dat) <- c("model", "draw")
  means <- ddply(dat, "model", summarise, grp.mean=mean(draw, na.rm = TRUE))
  # Base plot with densities
  if (monochrome) {
    gg <- ggplot(dat, aes(x=draw, linetype=model)) +
      scale_linetype_manual(values = c("solid", "dashed", "dotted"))
    
    # Add mean lines
    gg <- gg + geom_vline(data=means, aes(xintercept=grp.mean, linetype=model))
  } else {
    gg <- ggplot(dat, aes(x=draw, fill=model))
    # Add mean lines
    gg <- gg + geom_vline(data=means, aes(xintercept=grp.mean, color=model))
    gg <- gg + scale_fill_manual(values=active.palette) # color settings
    gg <- gg + scale_color_manual(values=active.palette) # color settings
  }
  gg <- gg +
    coord_cartesian(xlim=xlim) + 
    geom_density(bw=bw, alpha=.5) +
    xlim(xlim)

  # Set grids, axis, and axis texts
  gg <- gg + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(fill="white", color="white"),
                   #panel.border = element_blank(),
                   axis.title = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.text.x = element_text(size=30))
  # Add legend
  gg <- gg  + theme(legend.title = element_blank(),#element_text(size=40, face="bold"),
                    legend.text=element_text(size=30),
                    legend.justification=c(0,1),
                    legend.position=c(0, 1),
                    legend.background = element_blank(),
                    legend.key.size = unit(2, 'lines')) #+
  #guides(shape = guide_legend(override.aes = list(size = 5)))
  plot(gg)
  if (!is.null(filename)) {
    dev.copy2pdf(file=filename, out.type = "pdf")
  }
}
# 
# print(system.time(plot_rmse_gg(trial60$beta_totals[[1]],
#                                trial60$beta_totals[[4]],
#                                trial60$known_beta[1, ])))


print_beta_rmse <- function(beta_totals, known, name, sleeptime = 1.0, monochrome = FALSE) {
  if (monochrome) {
    folder <- "images_mc"
  } else {
    folder <- "images"
  }
    plot_rmse_gg(beta_totals[[1]],
                 beta_totals[[4]],
                 known,
                 filename = paste0(folder,"/rmse_",name[1],".pdf"),
                 monochrome = monochrome)
    Sys.sleep(sleeptime)
    plot_rmse_gg(beta_totals[[2]],
                 beta_totals[[4]],
                 known,
                 filename = paste0(folder,"/rmse_",name[2],".pdf"),
                 model.type = "Ordered",
                 monochrome = monochrome)
    Sys.sleep(sleeptime)
    plot_rmse_gg(beta_totals[[3]],
                 beta_totals[[4]],
                 known,
                 filename = paste0(folder,"/rmse_",name[3],".pdf"),
                 model.type = "Ordered",
                 monochrome = monochrome)
    Sys.sleep(sleeptime)
}


print_beta_rmse(trial60$beta_totals, trial60$known_beta, name = paste0("beta_60_sim_", simnames), monochrome = TRUE)
print_beta_rmse(trial60$beta_totals, trial60$known_beta, name = paste0("beta_60_sim_", simnames), monochrome = FALSE)

print_beta_rmse(trial600noshrink$beta_totals, trial600noshrink$known_beta, name = paste0("beta_600_sim_", simnames), monochrome = TRUE)
print_beta_rmse(trial600noshrink$beta_totals, trial600noshrink$known_beta, name = paste0("beta_600_sim_", simnames), monochrome = FALSE)



print_gamma_rmse <- function(gamma_totals, known, name, sleeptime = 1.0, monochrome = FALSE) {
  if (monochrome) {
    folder <- "images_mc"
  } else {
    folder <- "images"
  }
  plot_rmse_gg(gamma_totals[[1]],
               gamma_totals[[4]],
               known[[1]],
               filename = paste0(folder,"/rmse_",name[1],".pdf"), 
               monochrome = monochrome)
  Sys.sleep(sleeptime)
  plot_rmse_gg(gamma_totals[[2]],
               gamma_totals[[5]],
               known[[2]],
               filename = paste0(folder,"/rmse_",name[2],".pdf"),
               model.type = "Ordered", 
               monochrome = monochrome)
  Sys.sleep(sleeptime)
  plot_rmse_gg(gamma_totals[[3]],
               gamma_totals[[6]],
               known[[3]],
               filename = paste0(folder,"/rmse_",name[3],".pdf"),
               model.type = "Ordered", 
               monochrome = monochrome)
  Sys.sleep(sleeptime)
}


print_gamma_rmse(trial60$gamma_totals, trial60$known_gamma, name = paste0("gamma_60_sim_", simnames), monochrome = TRUE)
print_gamma_rmse(trial60$gamma_totals, trial60$known_gamma, name = paste0("gamma_60_sim_", simnames), monochrome = FALSE)

print_gamma_rmse(trial600noshrink$gamma_totals, trial600noshrink$known_gamma, name = paste0("gamma_600_sim_", simnames), monochrome = TRUE)
print_gamma_rmse(trial600noshrink$gamma_totals, trial600noshrink$known_gamma, name = paste0("gamma_600_sim_", simnames), monochrome = FALSE)














array_split <- function(data, number_of_chunks) {
  rowIdx <- seq_len(nrow(data))    
  t(sapply(split(rowIdx, cut(rowIdx, pretty(rowIdx, number_of_chunks))), function(x) colMeans(data[x, ])))
}

beta_means <- list()
beta_mean_error <- list()
print_beta_rmse_diff <- function(beta_totals, known, name, sleeptime = 1.0, monochrome = FALSE) {
  if (monochrome) {
    folder <- "images_mc"
  } else {
    folder <- "images"
  }
  for (i in 1:4) {
    beta_means[[i]] <- array_split(beta_totals[[i]], 500)
    beta_mean_error[[i]] <- sapply(1:500, function(j) {
      return(sqrt(mean((beta_means[[i]][j, ] - known[j, ])^2)))
    })
  }
  plot_comparison_gg(beta_mean_error[[1]],
                     beta_mean_error[[4]],
                     bw = NULL,
                     compare = difference,
                     filename = paste0(folder,"/rmse_beta_diff_",name[1],".pdf"),
                     monochrome = monochrome)
  Sys.sleep(sleeptime)
  plot_comparison_gg(beta_mean_error[[2]],
                     beta_mean_error[[4]],
                     bw = NULL,
                     compare = difference,
                     model.type = "Ordered",
                     filename = paste0(folder,"/rmse_beta_diff_",name[2],".pdf"),
                     monochrome = monochrome)
  Sys.sleep(sleeptime)
  plot_comparison_gg(beta_mean_error[[3]],
                     beta_mean_error[[4]],
                     bw = NULL,
                     compare = difference,
                     model.type = "Ordered",
                     filename = paste0(folder,"/rmse_beta_diff_",name[3],".pdf"),
                     monochrome = monochrome)
  Sys.sleep(sleeptime)
}

print_beta_rmse_diff(trial60$beta_totals, trial60$known_beta, name = paste0("beta_60_sim_", simnames), monochrome = TRUE)
print_beta_rmse_diff(trial60$beta_totals, trial60$known_beta, name = paste0("beta_60_sim_", simnames), monochrome = FALSE)




print_beta_rmse_ratio <- function(beta_totals, known, name, sleeptime = 1.0, monochrome = FALSE) {
  if (monochrome) {
    folder <- "images_mc"
  } else {
    folder <- "images"
  }
  for (i in 1:4) {
    beta_means[[i]] <- array_split(beta_totals[[i]], 500)
    beta_mean_error[[i]] <- sapply(1:500, function(j) {
      return(sqrt(mean((beta_means[[i]][j, ] - known[j, ])^2)))
    })
  }
  plot_comparison_gg(beta_mean_error[[1]],
                     beta_mean_error[[4]],
                     bw = NULL,
                     compare = ratio,
                     filename = paste0(folder,"/rmse_beta_ratio_",name[1],".pdf"),
                     monochrome = monochrome)
  Sys.sleep(sleeptime)
  plot_comparison_gg(beta_mean_error[[2]],
                     beta_mean_error[[4]],
                     bw = NULL,
                     compare = ratio,
                     model.type = "Ordered",
                     filename = paste0(folder,"/rmse_beta_ratio_",name[2],".pdf"),
                     monochrome = monochrome)
  Sys.sleep(sleeptime)
  plot_comparison_gg(beta_mean_error[[3]],
                     beta_mean_error[[4]],
                     bw = NULL,
                     compare = ratio,
                     model.type = "Ordered",
                     filename = paste0(folder,"/rmse_beta_ratio_",name[3],".pdf"),
                     monochrome = monochrome)
  Sys.sleep(sleeptime)
}

print_beta_rmse_ratio(trial60$beta_totals, trial60$known_beta, name = paste0("beta_60_sim_", simnames), monochrome = TRUE)
print_beta_rmse_ratio(trial60$beta_totals, trial60$known_beta, name = paste0("beta_60_sim_", simnames), monochrome = FALSE)


















plot_oprobit_gg <- function(cutpoints = c(-1, 2),
                            pal = c("black", "white", "black"),
                            xlim = c(-3.5, 3.5), filename = NULL) {
  allpoints <- c(xlim[1], cutpoints, xlim[2])
  label_height <- 0.12
  gg <- ggplot(data = data.frame(x = xlim), aes(x)) 
  for (i in 1:(length(allpoints)-1)) {
    gg <- gg + geom_area(stat = "function", fun = dnorm, args = list(mean = 0, sd = 1),
                         xlim = c(allpoints[i], allpoints[i+1]),
                         color = "black",
                         fill = pal[i],
                         alpha = 0.2) +
      ylab("") +
      scale_y_continuous(breaks = NULL) +
      annotate("text", x=(allpoints[i]+allpoints[i+1])/2, y=label_height, label = paste("y =", i), size = 7)
    
  }
  # gg <- gg + geom_area(stat = "function", fun = dnorm, args = list(mean = 0, sd = 1),
  #                      xlim = c(cutpoints[1], cutpoints[2]),
  #                      color = "black",
  #                      fill = "white",
  #                      alpha = 0.2)
  # gg <- gg + geom_area(stat = "function", fun = dnorm, args = list(mean = 0, sd = 1),
  #                      xlim = c(cutpoints[2], xlim[2]),
  #                      color = "black",
  #                      fill = "black",
  #                      alpha = 0.2)
  labels <- c()
  cutpoints <- sort(c(cutpoints, 0))
  print(length(cutpoints))
  cp <- 1
  for (i in 1:length(cutpoints)) {
    if (cutpoints[i] == 0) {
      exp <- bquote(expression(paste('x'^'T',beta)))
    } else {
      exp <- bquote(expression(paste(gamma,''[.(cp)])))
      cp = cp+1
    }
    labels <- c(labels, eval(exp))
  }
  
  cutpoints <- data.frame(cutpoints=cutpoints, type = rep(1, length(cutpoints)))
  cutpoints$type <- as.factor(cutpoints$type)
  print(str(cutpoints))
  meanpoint <- data.frame(position=c(0))
  gg <- gg + geom_vline(data=meanpoint, aes(xintercept=position), linetype=2)
  # Set grids, axis, and axis texts
  gg <- gg + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_rect(fill="white", color="white"),
                   axis.title = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.ticks.x = element_blank(),
                   legend.position = "none",
                   axis.text.x = element_text(size=30))
  
  gg <- gg + scale_x_continuous(breaks = cutpoints$cutpoints, labels = labels)
  gg <- gg + guides(fill=FALSE)
  # gg <- gg + annotate("text", x=-2.3, y=label_height, label = "y = 1", size = 7) +
  #   annotate("text", x=0.4, y=label_height, label = "y = 2", size = 7) +
  #   annotate("text", x=2.8, y=label_height, label = "y = 3", size = 7)
  plot(gg)
  if (!is.null(filename)) {
    dev.copy2pdf(file=filename, out.type = "pdf")
  }
}
plot_oprobit_gg(filename = "images/example_oprobit.pdf")
plot_oprobit_gg(cutpoints = c(1), pal = c("white", "black"), filename = "images/example_probit2.pdf")
plot_oprobit_gg(filename = "images_mc/example_oprobit.pdf")
plot_oprobit_gg(cutpoints = c(1), pal = c("white", "black"), filename = "images_mc/example_probit2.pdf")


plot_probit_gg <- function(xlim = c(-3.5, 3.5), filename = NULL) {
  
  gg <- ggplot(data = data.frame(x = xlim), aes(x)) +
    # geom_area(stat = "function", fun = dnorm, args = list(mean = 0, sd = 1),
    #           xlim = c(xlim[1], 0),
    #           color = "black",
    #           fill = "red",
    #           alpha = 0.2) +
    # geom_area(stat = "function", fun = dnorm, args = list(mean = 0, sd = 1),
    #           xlim = c(0, xlim[2]),
    #           color = "black",
    #           fill = "darkgreen",
    #           alpha = 0.2) +
    # geom_area(stat = "function", fun = dnorm, args = list(mean = 0, sd = 1),
    #           xlim = xlim,
    #           color = "black",
    #           fill = "white") +
    labs(x = expression(paste("x"^"T", beta)), y = expression(paste(Phi, "(x"^"T", beta, ")")))
  gg <- gg + stat_function(fun = function(x) {return(pnorm(x, mean = 0, sd = 1)/2)}, n = 101)

#  gg <- gg + geom_vline(data=cutpoints, aes(xintercept=cutpoints, linetype="dotted"))
  limits <- c(0, 0.5)
  limits <- as.data.frame(limits)
  gg <- gg + scale_y_continuous(breaks = limits$limits,
                                labels = paste("y = ", limits$limits*2),
                                position = "right")
  # Set grids, axis, and axis texts
  gg <- gg + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   #panel.border = element_blank(),
                   panel.background = element_rect(fill="white", color="white"),
                   axis.title.x = element_text(size = 24),
                   axis.title.y = element_blank(), #element_text(angle = 180, size = 24),
                   # axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   # axis.ticks.x = element_blank(),
                   legend.position = "none",
                   axis.text.x = element_text(size=24),
                   axis.text.y = element_text(size=24))
  

  # gg <- gg + guides(fill=FALSE)
  gg <- gg + annotate("text", x=-3, y=0.45,
                      label = expression(paste(Phi, "(x"^"T", beta, ")")),
                      size = 10)
  #   annotate("text", x=0.4, y=label_height, label = "y = 2", size = 7) +
  #   annotate("text", x=2.8, y=label_height, label = "y = 3", size = 7)
  plot(gg)
  if (!is.null(filename)) {
    dev.copy2pdf(file=filename, out.type = "pdf")
  }
}
plot_probit_gg(filename = "images/example_probit.pdf")
plot_probit_gg(filename = "images_mc/example_probit.pdf")













all_data <- readRDS("8_chains_fitted_on_all_text_data_prec10.RDS") #originally "8_chains_fitted_on_all_text_data.RDS"
all_beta_draws <- all_data[[1]]$simulations[[4]]$beta
for (i in 2:8) {
  all_beta_draws <- rbind(all_beta_draws, all_data[[i]]$simulations[[4]]$beta)
}

gamma_draws = list()
gamma_draws[[1]] <- all_data[[1]]$simulations[[4]]$gammas[[1]]
for (i in 2:8) {
  gamma_draws[[1]] <- rbind(gamma_draws[[1]], all_data[[i]]$simulations[[4]]$gammas[[1]])
}
colnames(gamma_draws[[1]]) <- c("gamma1")
gamma_draws[[2]] <- all_data[[1]]$simulations[[4]]$gammas[[2]]
for (i in 2:8) {
  gamma_draws[[2]] <- rbind(gamma_draws[[2]], all_data[[i]]$simulations[[4]]$gammas[[2]])
}
colnames(gamma_draws[[2]]) <- paste0("gamma", c(1, 2, 3))
gamma_draws[[3]] <- all_data[[1]]$simulations[[4]]$gammas[[3]]
for (i in 2:8) {
  gamma_draws[[3]] <- rbind(gamma_draws[[3]], all_data[[i]]$simulations[[4]]$gammas[[3]])
}
colnames(gamma_draws[[3]]) <- paste0("gamma", c(1, 2, 3))



separate_gamma_draws = list()
separate_gamma_draws[[1]] <- all_data[[1]]$simulations[[1]]$gammas[[1]]
for (i in 2:8) {
  separate_gamma_draws[[1]] <- rbind(separate_gamma_draws[[1]], all_data[[i]]$simulations[[1]]$gammas[[1]])
}
colnames(separate_gamma_draws[[1]]) <- c("gamma1")
separate_gamma_draws[[2]] <- all_data[[1]]$simulations[[2]]$gammas[[1]]
for (i in 2:8) {
  separate_gamma_draws[[2]] <- rbind(separate_gamma_draws[[2]], all_data[[i]]$simulations[[2]]$gammas[[1]])
}
colnames(separate_gamma_draws[[2]]) <- paste0("gamma", c(1, 2, 3))
separate_gamma_draws[[3]] <- all_data[[1]]$simulations[[3]]$gammas[[1]]
for (i in 2:8) {
  separate_gamma_draws[[3]] <- rbind(separate_gamma_draws[[3]], all_data[[i]]$simulations[[3]]$gammas[[1]])
}
colnames(separate_gamma_draws[[3]]) <- paste0("gamma", c(1, 2, 3))



# Plot posterior estimates for beta coefficients.
#
# Arguments:
# values: A matrix or data frame where each column contains posterior draws for a 
#         coefficient.
# bw: Bandwidth for the density estimation.
# color_gradient: A vector of two colors for the gradient (negative to positive).
# filename: If provided, the plots will be saved to files with this prefix.
# interactive: If TRUE, the function will wait for user input before proceeding to the next plot
# monochrome: If TRUE, the plots will be in black and white.
plot_posterior_estimates <- function(values, 
                                     bw = 0.15,
                                     color_gradient = c("red3", "green3"),
                                     filename = NULL,
                                     interactive = FALSE,
                                     monochrome = FALSE) {
  # If monochrome is requested, set the color gradient to black and white
  if (monochrome) {
    color_gradient = c("white", "black")
  }
  # Create a color palette for the density fill
  pal <- colorRampPalette(color_gradient)(100)
  print(pal) # Print the palette for debugging
  # Convert the input values to a data frame for plotting
  dat = data.frame(values)
  val <- 0
  # Loop over each parameter (column) in the values matrix/data frame
  for (i in 1:ncol(values)) {
    # Calculate the proportion of samples greater than zero for the parameter
    posneg <- mean(values[, i] > 0)
    print(posneg) # Print the proportion for debugging
    name = "test"
    fillcolor = "pink"
    # Determine the fill color and label based on the sign certainty
    if (posneg > 0.5) {
      # More samples are positive: use a color from the positive end of the palette
      val <- floor(200*(posneg-0.5))
      fillcolor <- pal[max(1, val)]
      name <- paste("Positive with", val, "% certainty")
    } else if (posneg < 0.5) {
      # More samples are negative: use a color from the negative end of the palette
      val <- floor(200*(0.5-posneg))
      fillcolor <- pal[max(1, val)]
      name <- paste("Negative with", val, "% certainty")
    } else {
      # Uncertain sign: use the first color in the palette
      print("HERERERERERERER!!!")
      fillcolor <- pal[1]
      name <- "Uncertain"
    }
    print(str(dat[, i])) # Print the structure of the parameter samples
    print(name) # Print the label
    print(fillcolor) # Print the fill color
    # Create the density plot for the parameter
    gg <- ggplot(dat, aes(x=!!dat[, i])) +
      geom_density(bw=bw, fill=fillcolor, alpha=.6) + # Density plot with fill
      # labs(x = paste(colnames(values)[i], name)) +
      labs(x = name) + # X-axis label with certainty
      geom_vline(xintercept = 0) + # Vertical line at zero
      theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_rect(fill="white", color="white"),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_blank(), #element_text(angle = 180, size = 24),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.ticks.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size=30))
        #axis.text.y = element_text(size=24))
    # Optionally annotate with the parameter name (commented out)
    # gg <- gg + annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, label = colnames(values)[i], size = 10)
    plot(gg) # Display the plot
    # If a filename is provided, save the plot to a PDF file
    if (!is.null(filename)) {
      dev.copy2pdf(file=paste0(filename, "_", val, "_", gsub("\\.", "dot", colnames(values)[i]), ".pdf"), out.type = "pdf")
    }
    # If interactive mode is enabled, wait for user input before next plot
    if (!i==ncol(values) & interactive){
      invisible(readline(prompt="Press [enter] to continue"))
    }
  }
}

# colnames(trial_text$beta_totals[[4]]) <- colnames(all_beta_draws)
# colnames(trial_text$beta_means[[4]]) <- colnames(all_beta_draws)


plot_posterior_estimates(all_beta_draws, filename = "images/beta_posterior", interactive = FALSE)
plot_posterior_estimates(all_beta_draws, filename = "images_mc/beta_posterior", interactive = FALSE, monochrome = TRUE)

#plot_posterior_estimates(trial_text$beta_totals[[4]], filename = "images/beta_posterior_2thirds", interactive = FALSE)


# Plot posterior estimates for gamma parameters.
# 
# Arguments:
# scale1: A matrix or data frame where each column contains posterior draws for the first 
#         set of gamma parameters.
# scale2: A matrix or data frame where each column contains posterior draws for the second 
#         set of gamma parameters.
# scale3: A matrix or data frame where each column contains posterior draws for the third
#         set of gamma parameters.
# colors: A vector of colors for the different scales.
# names: A vector of names for the different scales.
# xlim: A vector specifying the x-axis limits for the plot.
# filename: If provided, the plot will be saved to a file with this name.
# monochrome: If TRUE, the plot will be in black and white.
plot_posterior_gammas <- function(scale1, scale2, scale3,
                                  colors = cbbPalette2[c(6,2,8)],
                                  names = c("1", "2"),
                                  xlim = NULL,
                                  filename = NULL,
                                  monochrome = FALSE) {
  # Melt and label the first scale's gamma samples for plotting
  scale1 = melt(scale1)
  scale1 = cbind(scale1[, -1], scale=names[1])
  print(str(scale1))
  # Melt and label the second scale's gamma samples for plotting
  scale2 = melt(scale2)
  scale2 = cbind(scale2[, -1], scale=names[2])
  print(str(scale2))
  # Melt and label the third scale's gamma samples for plotting
  scale3 = melt(scale3)
  scale3 = cbind(scale3[, -1], scale=names[3])
  print(str(scale2))
  # Combine all scales into one data frame for plotting
  dat <- rbind(rbind(scale1, scale2), scale3)
  # Convert the scale column to a factor for ggplot
  dat[, 3] <- factor(dat[, 3])#, labels = c(1, 2), levels(c("Scale 1", "Scale 2")))
  colnames(dat)[1] = "variable"
  print(str(dat))
  # Choose plot style: monochrome (linetype) or color (fill/color)
  if (monochrome) {
    gg <- ggplot(dat, aes(x=value, linetype=scale)) +
      scale_linetype_manual(values = c("solid", "dashed", "dotted"))
  } else {
    gg <- ggplot(dat, aes(x=value, fill=scale, color=scale))
    gg <- gg + scale_fill_manual(values=colors) # color settings
    gg <- gg + scale_color_manual(values=colors) # color settings
  }
  # Plot the density of each gamma parameter for each scale
  gg <- gg +
    geom_density(data=dat, mapping=aes(group=paste0(scale,variable)), alpha=0.3, bw = 1) +
    xlab(expression(paste('x'^'T',beta)))
  # Set plot theme and axis formatting
  gg <- gg + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   #panel.border = element_blank(),
                   panel.background = element_rect(fill="white", color="white"),
                   axis.title.x = element_text(size = 24),
                   axis.title.y = element_blank(), #element_text(angle = 180, size = 24),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   # axis.ticks.x = element_blank(),
                   axis.text.x = element_text(size=24))
   gg <- gg  + theme(legend.title = element_blank(),
                     legend.text=element_text(size=30),
                     legend.justification=c(0, 1), 
                     legend.position=c(0, 1),  
                     legend.background = element_blank(),
                     legend.key.size = unit(2, 'lines')) #+
  # Optionally set the x-axis limits
  if (!is.null(xlim)) {
    gg <- gg + coord_cartesian(xlim = xlim)
  }
  #axis.text.y = element_text(size=24))
  plot(gg) # Display the plot
  # Optionally save the plot to a PDF file
  if (!is.null(filename)) {
    dev.copy2pdf(file=paste0(filename, ".pdf"), out.type = "pdf")
  }
}

xlim <- c(Reduce(min, gamma_draws), Reduce(max, gamma_draws))
xlim_sep <- c(Reduce(min, separate_gamma_draws), Reduce(max, separate_gamma_draws))


plot_posterior_gammas(gamma_draws[[1]], gamma_draws[[2]], gamma_draws[[3]],
                      names = c("LL", "Legimus", "Hegas"),
                      colors = cbbPalette2[c(1,4,3)],
                      xlim = xlim,
                      filename = "images/posterior_gamma_combined")
plot_posterior_gammas(gamma_draws[[1]], gamma_draws[[2]], gamma_draws[[3]],
                      names = c("LL", "Legimus", "Hegas"),
                      colors = cbbPalette2[c(1,4,3)],
                      xlim = xlim,
                      filename = "images_mc/posterior_gamma_combined",
                      monochrome = TRUE)

plot_posterior_gammas(separate_gamma_draws[[1]], separate_gamma_draws[[2]], separate_gamma_draws[[3]],
                      names = c("LL", "Legimus", "Hegas"),
                      colors = cbbPalette2[c(1,4,3)], 
                      xlim = xlim,
                      filename = "images/posterior_gamma_separate")
plot_posterior_gammas(separate_gamma_draws[[1]], separate_gamma_draws[[2]], separate_gamma_draws[[3]],
                      names = c("LL", "Legimus", "Hegas"),
                      colors = cbbPalette2[c(1,4,3)], 
                      xlim = xlim,
                      filename = "images_mc/posterior_gamma_separate",
                      monochrome = TRUE)


# 
# set_palette <- cbbPalette2[c(6,2,8)]
# 
# print_posterior_gammas <- function(gamma_draws, palette = set_palette, names = corpus_names, sleeptime=1.0) {
#   xlim <- c(Reduce(min, gamma_draws), Reduce(max, gamma_draws))
#   for (i in 1:3) {
#     j <- (i+1)
#     if (j == 4) {
#       j <- 1
#     }
#     plot_posterior_gammas(gamma_draws[[i]], gamma_draws[[j]],
#                           names = c(names[i], names[j]),
#                           colors = palette[c(i, j)],
#                           xlim = xlim)
#     Sys.sleep(sleeptime)
#   }
# }
# 
# print_posterior_gammas(gamma_draws)
# 
# 
# plot_posterior_gammas(gamma_draws[[1]], gamma_draws[[2]], gamma_draws[[3]],
#                       names = corpus_names,
#                       colors = set_palette[c(1,2,3)])
# plot_posterior_gammas(gamma_draws[[2]], gamma_draws[[3]])
# plot_posterior_gammas(gamma_draws[[3]], gamma_draws[[1]])
# 
# 
# #Testing color gradients
# colfunc <- colorRampPalette(c("red3", "green3"))
# plot(rep(1,100),col=colfunc(100),pch=19,cex=3)




for (i in seq(1,16,1)) {
  print(paste0(features[i], " & ", features[i+16], " & ", features[i+32], " \\"))
}

