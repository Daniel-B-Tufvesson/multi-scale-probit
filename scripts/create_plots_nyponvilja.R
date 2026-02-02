# MOVE FILE TO scripts/create_plots_nyponvilja.R

library(ggplot2)
#install.packages("tikzDevice")
library(grDevices)
library(reshape2)
library(plyr)
library(coda)

source("util_nyponvilja.R")

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(1,2,4,5,6,7,8,3)]

# The palette with black:
cbbPalette1 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#cbbPalette1 <- c("#bb11ee", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbbPalette2 <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
                "#CC79A7", "#F0E442")

plot_stats_gg <- function(separate, combined,
                       xlim=c(0, 1), ylim=c(0,4),
                       bw = "nrd0", set_names = NULL,
                       filename = NULL, model.type = "Ordered",
                       custom.palette = NULL,
                       monochrome = TRUE,
                       i_value = NULL,
                       baseline_filepath = NULL,
                       dataset = NULL) {
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
  # Organize data
  #print("1")
  dat <- data.frame(combined, separate)
  names(dat) <- c("Multi-Scale", "Single Scale")
  #print("1.5")
  print(str(dat))
  dat <- melt(dat)
  colnames(dat) <- c("model", "draw")
  #print("3")
  means <- ddply(dat[which(dat[, 2] >= 0 & dat[, 2] <= 1), ], "model", summarise, grp.mean=mean(draw, na.rm = TRUE))
  #print(str(dat))
  # Base plot with densities
  if (monochrome) {
    gg <- ggplot(dat, aes(x=draw, linetype=model)) +
      scale_linetype_manual(values = c("solid", "dashed", "dotted"))
    # Add mean lines
    gg <- gg + geom_vline(data=means, aes(xintercept=grp.mean, linetype=model), size = 1.5)
  } else {
    gg <- ggplot(dat, aes(x=draw, fill=model))
    # Add mean lines
    gg <- gg + geom_vline(data=means, aes(xintercept=grp.mean, color=model), size = 1.5)
    gg <- gg + scale_fill_manual(values=active.palette) # color settings
    gg <- gg + scale_color_manual(values=active.palette) # color settings
  }
  if (i_value == 2 || i_value == 6) {
    xlim = c(0,1)
  } else {
    xlim = c(0,1)
  }
  gg <- gg + labs(x = paste(dataset))
  gg <- gg + theme(panel.background = element_rect(fill="white", color="white"))
  gg <- gg + coord_cartesian(xlim=xlim) + 
    geom_density(bw=bw, alpha=.5) +
    xlim(xlim)
  #print("5")
  # Set grids, axis, and axis texts
  gg <- gg + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   #panel.border = element_blank(),
                   axis.title.y = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.text.x = element_text(size=30),
                   axis.title.x = element_text(size = 40, face = "italic"))
                   #plot.margin = margin(t = 0, r = 4, b = 0, l = 1, unit = "cm"))  # Add margin to the right side)
  # Add legend
  gg <- gg  + theme(legend.title = element_blank(),
                    legend.text=element_text(size=30),
                    legend.justification=c(0, 1),
                    legend.position=c(0, 1),
                    legend.background = element_blank(),
                    legend.key.size = unit(2, 'lines'))
    #guides(shape = guide_legend(override.aes = list(size = 5)))
  # Load the .RDS file containing the vector
  F1_baseline <- readRDS(baseline_filepath)
  #print(F1_baseline)
  # Calculate the mean value of the vector
  mean_F1_baseline <- mean(F1_baseline)
  print(mean_F1_baseline)
  #print(mean_F1_baseline)
  # Check if i is 1 or 2, meaning F1 train or F1 test
  if (i_value == 1 || i_value == 2) {
    # Add vertical line for mean baseline F1
    gg <- gg + geom_vline(xintercept = mean_F1_baseline, color = "grey", size = 1.5)
  }
  plot(gg)
  if (!is.null(filename)) {
    dev.copy2pdf(file=filename, out.type = "pdf", width = 10, height = 5)
  }
}


print_totals <- function(model_totals, name, sleeptime = 1.0, monochrome = TRUE) {
  if (monochrome) {
    folder <- "images_mc_22cov_withsuc/metric_totals"
  } else {
    folder <- "images22cov/images500withsuc/metric_totals"
  }
  for (i in c(1,2,4,6)) {
    print(metric_names[i])
    plot_stats_gg(model_totals[[1]][, i], #Change the 1,3 and 2,4 depending on how many datasets.
                  model_totals[[3]][, i], #If 3 datasets, it should be 1,4 + 2,5 + 3,6 (add another plot code chunk)
                  filename = paste0(folder,"/",metric_names[i],"_",name[1],".pdf"),
                  model.type = "Ordered",
                  monochrome = monochrome,
                  i_value = i,
                  baseline_filepath = "baseline/nypon_F1_withsuc.RDS",
                  dataset = "Nypon")
    Sys.sleep(sleeptime)
    plot_stats_gg(model_totals[[2]][, i],
                  model_totals[[4]][, i],
                  filename = paste0(folder,"/",metric_names[i],"_",name[2],".pdf"),
                  model.type = "Ordered",
                  monochrome = monochrome,
                  i_value = i,
                  baseline_filepath = "baseline/vilja_F1_withsuc.RDS",
                  dataset = "Vilja")
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
simnames <- c("ordered1", "ordered2")
corpus_names <- c("Nypon", "Vilja") #Nypon is first, Vilja second

trial_text <- readRDS("text_2thirds_dataset_all500_22cov_withsuc_reorganized.RDS")
trial_text_small <- readRDS("text_1thirds_dataset_all500_22cov_withsuc_reorganized.RDS")
#summary(trial_text$metric_totals[[4]], na.rm = TRUE)
#print(trial_text$metric_totals)
print_totals(trial_text$metric_totals, name = paste0("text_", corpus_names), monochrome = TRUE)
#print_totals(trial_text$metric_totals, name = paste0("text_2third", corpus_names), monochrome = FALSE)

print_totals(trial_text_small$metric_totals, name = paste0("text_1third_", corpus_names), monochrome = TRUE)
#print_totals(trial_text_small$metric_totals, name = paste0("text_1third_", corpus_names), monochrome = FALSE)



#----------------------METRIC MEANS-------------------------------------


# 
#  metric_names <- strsplit("F1_train,
#                           F1_test,
#                           jaccard_train,
#                           jaccard_test,
#                           youden_train,
#                           youden_test,
#                           clharm_train,
#                           clharm_test,
#                           kendall_train,
#                           kendall_test,
#                           spearman_train,
#                           spearman_test,
#                           rharm_train,
#                           rharm_test,
#                           harm_train,
#                           harm_test", ",
#                           ")[[1]]
# 
#  print_means <- function(model_totals, name, sleeptime = 1.0, monochrome = FALSE) {
#    if (monochrome) {
#      folder <- "images_mc"
#    } else {
#      folder <- "images500withsuc/metric_means"
#    }
#    for (i in c(1,2,4,6)) {
#      print(metric_names[i])
#      plot_stats_gg(model_totals[[1]][, i],
#                    model_totals[[3]][, i],
#                    bw = 0.05,
#                    filename = paste0(folder,"/",metric_names[i],"_means_",name[2],".pdf"),
#                    model.type = "Ordered",
#                  custom.palette = c("#000000", "#D55E00"),
#                  monochrome = monochrome,
#                  i_value = i,
#                  baseline_filepath = "baseline/nypon_F1_withsuc.RDS",
#                  dataset = "Nypon")
#      Sys.sleep(sleeptime)
#      plot_stats_gg(model_totals[[2]][, i],
#                    model_totals[[4]][, i],
#                    bw = 0.05,
#                    filename = paste0(folder,"/",metric_names[i],"_means_",name[3],".pdf"),
#                    model.type = "Ordered",
#                  custom.palette = c("#000000", "#D55E00"),
#                  monochrome = monochrome,
#                  i_value = i,
#                  baseline_filepath = "baseline/vilja_F1_withsuc.RDS",
#                  dataset = "Vilja")
#      Sys.sleep(sleeptime)
#    }
#  }
# 
#  trial_text <- readRDS("text_2thirds_dataset_all500_withsuc_reorganized.RDS")
#  trial_text_small <- readRDS("text_1thirds_dataset_all500_withsuc_reorganized.RDS")
#  #print_means(trial_text$metric_means, name = paste0("text_", corpus_names), monochrome = TRUE)
#  print_means(trial_text$metric_means, name = paste0("text_", corpus_names), monochrome = FALSE)


#----------------------MEANS DIFF-----------------------------




difference <- function(set1, set2) {
  d <- set1 - set2
  xmax <- max(abs(d))
  return(list(d,
              mean(d < 0),
              mean(d > 0),
              cutpoint = 0,
              xlim = c(-xmax, xmax)))
}

ratio <- function(set1, set2) {
  r <- set1 / set2
  xmax <- max(abs(r - 1))
  return(list(r,
              mean(r < 1),
              mean(r > 1),
              cutpoint = 1,
              xlim = c(-xmax, xmax) + 1))
}

plot_comparison_gg <- function(separate, combined,
                               bw = NULL,
                               filename = NULL,
                               model.type = "Binary",
                               compare = difference,
                               custom.palette = NULL,
                               monochrome = TRUE,
                               dataset = NULL) {
  if (is.null(custom.palette)) {
    if (model.type == "Binary") {
      active.palette <- cbbPalette2[c(2, 1)]
    } else {
      active.palette <- cbbPalette1[c(2, 1)]
    }
  } else {
    active.palette <- custom.palette
  }
  if (monochrome) {
    active.palette <- c("#CCCCCC", "#000000")
  }
  df <- data.frame(separate, combined)
  df <- df[complete.cases(df), ]
  separate <- df$separate
  combined <- df$combined
  comparison <- compare(combined, separate)
  if (is.null(bw)) {
    dat <- with(density(comparison[[1]]), data.frame(x, y))
  } else {
    dat <- with(density(comparison[[1]], bw = bw), data.frame(x, y))
  }
  if (all(dat$x > comparison[[4]])) {
    active.palette <- active.palette[c(2, 1)]
  }
  dat <- cbind(dat, Better = factor(as.numeric(dat$x>comparison[[4]]), levels = c(0, 1),
                                    labels = c(paste0(model.type, " (", round(100*comparison[[2]], 2), " %)"),
                                               paste0("Multi-Scale (", round(100*comparison[[3]], 2), " %)"))))
  gg <- ggplot(data = dat, mapping = aes(x = x, y = y)) +
    geom_line()+
    #geom_area(mapping = aes(x = ifelse(x <= 0 , x, 0)), fill = active.palette[[1]], alpha = 0.3) +
    #geom_area(mapping = aes(x = ifelse(x > 0 , x, 0)), fill = active.palette[[2]], alpha = 0.3) +
    geom_vline(data = data.frame(cutpoints = c(comparison[[4]])), aes(xintercept=cutpoints)) + 
    ylim(c(0, max(dat$y)))
  # if (monochrome) {
  #   gg <- gg + geom_area(mapping = aes(linetype = Better), alpha = 0.3) +
  #     scale_linetype_manual(values = c("solid", "dashed", "dotted"))
  # } else {
    gg <- gg +
      geom_area(mapping = aes(fill = Better, col = Better), alpha = 0.3)
  # }
    gg <- gg + labs(x = paste(dataset))
  gg <- gg + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   #panel.border = element_blank(),
                   panel.background = element_rect(fill="white", color="white"),
                   #axis.title.x = element_blank(),
                   axis.title.y = element_blank(), #element_text(angle = 180, size = 24),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   # axis.ticks.x = element_blank(),
                   axis.text.x = element_text(size=30),
                   axis.title.x = element_text(size = 40, face = "italic"))
  gg <- gg  + theme(legend.title = element_blank(),
                    legend.text=element_text(size=30),
                    legend.justification=c(1,1), 
                    legend.position=c(1, 1),
                    legend.background = element_blank(),
                    legend.key.size = unit(2, 'lines'))
  gg <- gg + scale_fill_manual(values=active.palette) # color settings
  gg <- gg + scale_color_manual(values=active.palette)
  gg <- gg + scale_x_continuous(limits = comparison$xlim)
  plot(gg)
  if (!is.null(filename)) {
    dev.copy2pdf(file=filename, out.type = "pdf",width = 10, height = 5)
  }
}


print_means_diff <- function(model_means, name, sleeptime = 1.0, monochrome = TRUE) {
    if (monochrome) {
      folder <- "images_mc_22cov_withsuc/means_diff"
    } else {
      folder <- "images22cov/images500withsuc/means_diff"
    }
  for (i in c(1, 2, 4, 6)) {
    print(metric_names[i])
    plot_comparison_gg(model_means[[1]][, i],
                       model_means[[3]][, i],
                       bw = NULL, 
                       filename = paste0(folder,"/",metric_names[i],"_means_diff_",name[1],".pdf"),
                       model.type = "Single Scale",
                       compare = difference,
                       monochrome = monochrome,
                       dataset = "Nypon")
    
    Sys.sleep(sleeptime)
    plot_comparison_gg(model_means[[2]][, i],
                       model_means[[4]][, i],
                       bw = NULL, 
                       filename = paste0(folder,"/",metric_names[i],"_means_diff_",name[2],".pdf"),
                       model.type = "Single Scale",
                       compare = difference,
                       monochrome = monochrome,
                       dataset = "Vilja")
    Sys.sleep(sleeptime)
  }  
}

trial_text <- readRDS("text_2thirds_dataset_all500_22cov_withsuc_reorganized.RDS")
trial_text_small <- readRDS("text_1thirds_dataset_all500_22cov_withsuc_reorganized.RDS")

#print(trial_text$metric_means[[4]][, 1])

print_means_diff(trial_text$metric_means, name = paste0("text_", corpus_names), monochrome = TRUE)
#print_means_diff(trial_text$metric_means, name = paste0("text_", corpus_names), monochrome = FALSE)

print_means_diff(trial_text_small$metric_means, name = paste0("text_1third", corpus_names), monochrome = TRUE)
#print_means_diff(trial_text_small$metric_means, name = paste0("text_1third", corpus_names), monochrome = FALSE)





#---------------------EXAMPLE PROBIT----------






# 
# 
# plot_oprobit_gg <- function(cutpoints = c(-1, 2),
#                             pal = c("black", "white", "black"),
#                             xlim = c(-3.5, 3.5), filename = NULL) {
#   allpoints <- c(xlim[1], cutpoints, xlim[2])
#   label_height <- 0.12
#   gg <- ggplot(data = data.frame(x = xlim), aes(x)) 
#   for (i in 1:(length(allpoints)-1)) {
#     gg <- gg + geom_area(stat = "function", fun = dnorm, args = list(mean = 0, sd = 1),
#                          xlim = c(allpoints[i], allpoints[i+1]),
#                          color = "black",
#                          fill = pal[i],
#                          alpha = 0.2) +
#       ylab("") +
#       scale_y_continuous(breaks = NULL) +
#       annotate("text", x=(allpoints[i]+allpoints[i+1])/2, y=label_height, label = paste("y =", i), size = 7)
#     
#   }
#   # gg <- gg + geom_area(stat = "function", fun = dnorm, args = list(mean = 0, sd = 1),
#   #                      xlim = c(cutpoints[1], cutpoints[2]),
#   #                      color = "black",
#   #                      fill = "white",
#   #                      alpha = 0.2)
#   # gg <- gg + geom_area(stat = "function", fun = dnorm, args = list(mean = 0, sd = 1),
#   #                      xlim = c(cutpoints[2], xlim[2]),
#   #                      color = "black",
#   #                      fill = "black",
#   #                      alpha = 0.2)
#   labels <- c()
#   cutpoints <- sort(c(cutpoints, 0))
#   print(length(cutpoints))
#   cp <- 1
#   for (i in 1:length(cutpoints)) {
#     if (cutpoints[i] == 0) {
#       exp <- bquote(expression(paste('x'^'T',beta)))
#     } else {
#       exp <- bquote(expression(paste(gamma,''[.(cp)])))
#       cp = cp+1
#     }
#     labels <- c(labels, eval(exp))
#   }
#   
#   cutpoints <- data.frame(cutpoints=cutpoints, type = rep(1, length(cutpoints)))
#   cutpoints$type <- as.factor(cutpoints$type)
#   print(str(cutpoints))
#   meanpoint <- data.frame(position=c(0))
#   gg <- gg + geom_vline(data=meanpoint, aes(xintercept=position), linetype=2, size = 1.5)
#   # Set grids, axis, and axis texts
#   gg <- gg + theme(panel.grid.major = element_blank(),
#                    panel.grid.minor = element_blank(),
#                    panel.border = element_blank(),
#                    panel.background = element_rect(fill="white", color="white"),
#                    axis.title = element_blank(),
#                    axis.text.y = element_blank(),
#                    axis.ticks.y = element_blank(),
#                    axis.ticks.x = element_blank(),
#                    legend.position = "none",
#                    axis.text.x = element_text(size=30))
#   
#   gg <- gg + scale_x_continuous(breaks = cutpoints$cutpoints, labels = labels)
#   gg <- gg + guides(fill=FALSE)
#   # gg <- gg + annotate("text", x=-2.3, y=label_height, label = "y = 1", size = 7) +
#   #   annotate("text", x=0.4, y=label_height, label = "y = 2", size = 7) +
#   #   annotate("text", x=2.8, y=label_height, label = "y = 3", size = 7)
#   plot(gg)
#   if (!is.null(filename)) {
#     dev.copy2pdf(file=filename, out.type = "pdf", width = 10, height = 5)
#   }
# }
# plot_oprobit_gg(filename = "example_oprobit.pdf")
# plot_oprobit_gg(cutpoints = c(1), pal = c("white", "black"), filename = "example_probit2.pdf")
# plot_oprobit_gg(filename = "example_oprobit.pdf")
# plot_oprobit_gg(cutpoints = c(1), pal = c("white", "black"), filename = "example_probit2.pdf")
# 
# 
# plot_probit_gg <- function(xlim = c(-3.5, 3.5), filename = NULL) {
#   
#   gg <- ggplot(data = data.frame(x = xlim), aes(x)) +
#     # geom_area(stat = "function", fun = dnorm, args = list(mean = 0, sd = 1),
#     #           xlim = c(xlim[1], 0),
#     #           color = "black",
#     #           fill = "red",
#     #           alpha = 0.2) +
#     # geom_area(stat = "function", fun = dnorm, args = list(mean = 0, sd = 1),
#     #           xlim = c(0, xlim[2]),
#     #           color = "black",
#     #           fill = "darkgreen",
#     #           alpha = 0.2) +
#     # geom_area(stat = "function", fun = dnorm, args = list(mean = 0, sd = 1),
#     #           xlim = xlim,
#     #           color = "black",
#     #           fill = "white") +
#     labs(x = expression(paste("x"^"T", beta)), y = expression(paste(Phi, "(x"^"T", beta, ")")))
#   gg <- gg + stat_function(fun = function(x) {return(pnorm(x, mean = 0, sd = 1)/2)}, n = 101)
# 
# #  gg <- gg + geom_vline(data=cutpoints, aes(xintercept=cutpoints, linetype="dotted"))
#   limits <- c(0, 0.5)
#   limits <- as.data.frame(limits)
#   gg <- gg + scale_y_continuous(breaks = limits$limits,
#                                 labels = paste("y = ", limits$limits*2),
#                                 position = "right")
#   # Set grids, axis, and axis texts
#   gg <- gg + theme(panel.grid.major = element_blank(),
#                    panel.grid.minor = element_blank(),
#                    #panel.border = element_blank(),
#                    panel.background = element_rect(fill="white", color="white"),
#                    axis.title.x = element_text(size = 24),
#                    axis.title.y = element_blank(), #element_text(angle = 180, size = 24),
#                    # axis.text.y = element_blank(),
#                    axis.ticks.y = element_blank(),
#                    # axis.ticks.x = element_blank(),
#                    legend.position = "none",
#                    axis.text.x = element_text(size=24),
#                    axis.text.y = element_text(size=24))
#   
# 
#   # gg <- gg + guides(fill=FALSE)
#   gg <- gg + annotate("text", x=-3, y=0.45,
#                       label = expression(paste(Phi, "(x"^"T", beta, ")")),
#                       size = 10)
#   #   annotate("text", x=0.4, y=label_height, label = "y = 2", size = 7) +
#   #   annotate("text", x=2.8, y=label_height, label = "y = 3", size = 7)
#   plot(gg)
#   if (!is.null(filename)) {
#     dev.copy2pdf(file=filename, out.type = "pdf")
#   }
# }
# plot_probit_gg(filename = "example_probit.pdf")
# plot_probit_gg(filename = "example_probit.pdf")
# 
# 



#---------------------------POSTERIOR ANALYSIS--------------







all_data <- readRDS("20_chains_fitted_on_all_text_data_22cov_withsuc.RDS") #originally "8_chains_fitted_on_all_text_data.RDS"
all_beta_draws <- all_data[[1]]$simulations[[3]]$beta
for (i in 2:20) {
  all_beta_draws <- rbind(all_beta_draws, all_data[[i]]$simulations[[3]]$beta)
}

gamma_draws = list()
#Nypon
gamma_draws[[1]] <- all_data[[1]]$simulations[[3]]$gammas[[1]]
for (i in 2:20) {
  gamma_draws[[1]] <- rbind(gamma_draws[[1]], all_data[[i]]$simulations[[3]]$gammas[[1]])
}
#colnames(gamma_draws[[1]]) <- paste0("gamma", c(1, 2, 3, 4)) #without suc
colnames(gamma_draws[[1]]) <- paste0("gamma", c(1, 2, 3, 4, 5)) #with suc
#Vilja
gamma_draws[[2]] <- all_data[[1]]$simulations[[3]]$gammas[[2]]
for (i in 2:20) {
  gamma_draws[[2]] <- rbind(gamma_draws[[2]], all_data[[i]]$simulations[[3]]$gammas[[2]])
}
#colnames(gamma_draws[[2]]) <- paste0("gamma", c(1, 2, 3, 4, 5)) #without suc
colnames(gamma_draws[[2]]) <- paste0("gamma", c(1, 2, 3, 4, 5, 6)) #with suc


separate_gamma_draws = list()
separate_gamma_draws[[1]] <- all_data[[1]]$simulations[[1]]$gammas[[1]]
for (i in 2:20) {
  separate_gamma_draws[[1]] <- rbind(separate_gamma_draws[[1]], all_data[[i]]$simulations[[1]]$gammas[[1]])
}
#colnames(separate_gamma_draws[[1]]) <- paste0("gamma", c(1, 2, 3, 4)) #without suc
colnames(separate_gamma_draws[[1]]) <- paste0("gamma", c(1, 2, 3, 4, 5)) #with suc
separate_gamma_draws[[2]] <- all_data[[1]]$simulations[[2]]$gammas[[1]]
for (i in 2:20) {
  separate_gamma_draws[[2]] <- rbind(separate_gamma_draws[[2]], all_data[[i]]$simulations[[2]]$gammas[[1]])
}
#colnames(separate_gamma_draws[[2]]) <- paste0("gamma", c(1, 2, 3, 4, 5)) #without suc
colnames(separate_gamma_draws[[2]]) <- paste0("gamma", c(1, 2, 3, 4, 5, 6)) #with suc
# print(separate_gamma_draws[[2]])
# 
# melted_df <- reshape2::melt(separate_gamma_draws[[2]])
# 
# # Create a distribution plot
# ggplot(melted_df, aes(x=value)) +
#   geom_density() +
#   labs(title="Distribution Plot", x="Values", y="Density")


plot_posterior_estimates <- function(values, bw = 0.05,
                                     color_gradient = c("red3", "green3"),
                                     filename = NULL,
                                     interactive = FALSE,
                                     monochrome = TRUE) {
  if (monochrome) {
    color_gradient = c("white", "black")
  }
  pal <- colorRampPalette(color_gradient)(100)
  print(pal)
  dat = data.frame(values)
  val <- 0
  for (i in 1:ncol(values)) {
    posneg <- mean(values[, i] > 0)
    print(posneg)
    name = "test"
    fillcolor = "pink"
    if (posneg > 0.5) {
      val <- floor(200*(posneg-0.5))
      fillcolor <- pal[max(1, val)]
      name <- paste("Positive with", val, "% certainty")
    } else if (posneg < 0.5) {
      val <- floor(200*(0.5-posneg))
      fillcolor <- pal[max(1, val)]
      name <- paste("Negative with", val, "% certainty")
    } else {
      print("HERERERERERERER!!!")
      fillcolor <- pal[1]
      name <- "Uncertain"
    }
    #print(str(dat[, i]))
    print(name)
    print(fillcolor)
    
    gg <- ggplot(dat, aes(x=!!dat[, i])) +
      coord_cartesian(clip = "off") +
      geom_density(bw=bw, fill=fillcolor, alpha=.6) +
      # labs(x = paste(colnames(values)[i], name)) +
      labs(x = name) +
      geom_vline(xintercept = 0) +
      theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_rect(fill="white", color="white"),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_blank(), #element_text(angle = 180, size = 24),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.ticks.x = element_blank(),
        # legend.title = element_blank(),
        # legend.text=element_text(size=30),
        # legend.justification=c(1, 1),
        # legend.position=c(1, 1),
        # legend.background = element_blank(),
        # legend.key.size = unit(2, 'lines'),
        axis.text.x = element_text(size=30),
        plot.title = element_text(hjust = 0.5, vjust = 1, size = 30, face = "italic"))
        #axis.text.y = element_text(size=24))
    #gg <- gg + labs(x = paste(colnames(values)[i]))
    #gg <- gg + annotate("text", x = Inf, y = -Inf, hjust = -1, vjust = 0, label = colnames(values)[i], size = 10)
    gg <- gg + ggtitle(colnames(values)[i])
    plot(gg)
    if (!is.null(filename)) {
      dev.copy2pdf(file=paste0(filename, "_", val, "_", gsub("\\.", "dot", colnames(values)[i]), ".pdf"), out.type = "pdf")
    }
    if (!i==ncol(values) & interactive){
      invisible(readline(prompt="Press [enter] to continue"))
    }
  }
}

# colnames(trial_text$beta_totals[[4]]) <- colnames(all_beta_draws)
# colnames(trial_text$beta_means[[4]]) <- colnames(all_beta_draws)

#plot_posterior_estimates(all_beta_draws, filename = "images_mc_22cov_withsuc/20chains/features_posteriors/beta_posterior", interactive = FALSE)
plot_posterior_estimates(all_beta_draws, filename = "images_mc_22cov_withsuc/20chains/features_posteriors/beta_posterior", interactive = FALSE, monochrome = TRUE)

#plot_posterior_estimates(trial_text$beta_totals[[4]], filename = "images/beta_posterior_2thirds", interactive = FALSE)



plot_posterior_gammas <- function(scale1, scale2,
                                  colors = cbbPalette2[c(6,2)],
                                  names = c("1", "2"),
                                  xlim = NULL,
                                  filename = NULL,
                                  monochrome = TRUE) {
  scale1 = melt(scale1)
  scale1 = cbind(scale1[, -1], scale=names[1])
  #print("scale 1")
  #print(str(scale1)$scale)
  scale2 = melt(scale2)
  scale2 = cbind(scale2[, -1], scale=names[2])
  #print("scale 2")
  #print(str(scale2)$scale)

  dat <- rbind(scale1, scale2)
  #print(dat[, 2])
  #print("new line")
  #print(print(dat[, 3]))
  dat[, 3] <- factor(dat[, 3])#, labels = c(1, 2), levels(c("Scale 1", "Scale 2")))
  colnames(dat)[1] = "variable"
  #print(str(dat))
  if (monochrome) {
    gg <- ggplot(dat, aes(x=value, linetype=scale)) +
      scale_linetype_manual(values = c("solid", "dashed", "dotted"))
  } else {
    gg <- ggplot(dat, aes(x=value, fill=scale, color=scale))
    gg <- gg + scale_fill_manual(values=colors) # color settings
    gg <- gg + scale_color_manual(values=colors) # color settings
  }
  gg <- gg +
    geom_density(data=dat, mapping=aes(group=paste0(scale,variable)), alpha=0.3) +
    xlab(expression(paste('x'^'T',beta)))
  gg <- gg + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   #panel.border = element_blank(),
                   panel.background = element_rect(fill="white", color="white"),
                   axis.title.x = element_text(size = 40),
                   axis.title.y = element_blank(), #element_text(angle = 180, size = 24),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   # axis.ticks.x = element_blank(),
                   axis.text.x = element_text(size=30))
   gg <- gg  + theme(legend.title = element_blank(),
                     legend.text=element_text(size=40),
                     legend.justification=c(0, 1), 
                     legend.position=c(0, 1),  
                     legend.background = element_blank(),
                     legend.key.size = unit(2, 'lines'),
                     plot.margin = margin(t = 0, r = 0, b = 0, l = 1, unit = "cm")) #+
  if (!is.null(xlim)) {
    gg <- gg + coord_cartesian(xlim = xlim)
  }
  #axis.text.y = element_text(size=24))
  plot(gg)
  if (!is.null(filename)) {
    dev.copy2pdf(file=paste0(filename, ".pdf"), out.type = "pdf", width = 10)
  }
}

#print(gamma_draws[[2]][, "gamma1", drop = FALSE])
xlim <- c(Reduce(min, gamma_draws), Reduce(max, gamma_draws))
xlim_sep <- c(Reduce(min, separate_gamma_draws), Reduce(max, separate_gamma_draws))
print(xlim)

plot_posterior_gammas(gamma_draws[[1]], gamma_draws[[2]],
                      names = c("Nypon", "Vilja"),
                      colors = cbbPalette2[c(5,7)],
                      xlim = xlim,
                      filename = "images_mc_22cov_withsuc/20chains/gamma_posteriors/posterior_gamma_combined")

plot_posterior_gammas(separate_gamma_draws[[1]], separate_gamma_draws[[2]],
                      names = c("Nypon", "Vilja"),
                      colors = cbbPalette2[c(5,7)], 
                      xlim = xlim_sep,
                      filename = "images_mc_22cov_withsuc/20chains/gamma_posteriors/posterior_gamma_separate")






#------------------------------------LATENT VALUE SCATTER----------------------



# Load required libraries
library(ggplot2)
library(dplyr)

# Read the data
data <- readRDS("latent_values_all20.RDS")

# Define the order of levels for nypon and vilja datasets
nypon_levels <- c("1", "2", "3", "4", "5")
vilja_levels <- c("x-small", "small", "medium", "large", "x-large", "xx-large")
print(vilja_data)

# Plot for nypon dataset
nypon_data <- data %>%
  filter(dataset == "nypon") %>%
  mutate(level = factor(level, levels = nypon_levels))

# Plot for vilja dataset
vilja_data <- data %>%
  filter(dataset == "vilja") %>%
  mutate(level = factor(level, levels = vilja_levels))

# Rename levels for vilja dataset
vilja_data$level <- factor(vilja_data$level,
                           levels = vilja_levels,
                           labels = c("X-S", "S", "M", "L", "X-L", "XX-L"))

head(nypon_data)
str(vilja_data$level)

# Calculate mean and standard deviation manually
nypon_summary <- nypon_data %>%
  group_by(level) %>%
  summarise(
    mean_latent_variable = mean(latent_variable),
    sd_latent_variable = sd(latent_variable)
  )

vilja_summary <- vilja_data %>%
  group_by(level) %>%
  summarise(
    mean_latent_variable = mean(latent_variable),
    sd_latent_variable = sd(latent_variable)
  )

# Define function to create the plot
latent_plot <- function(data, summary, title) {
  plot <- ggplot() +
    geom_point(data = data, aes(x = level, y = latent_variable, color = level), size = 2, color = "#999999", show.legend = FALSE) +
    geom_errorbar(data = summary, aes(x = level, ymin = mean_latent_variable - sd_latent_variable, ymax = mean_latent_variable + sd_latent_variable), width = 0.4, color = "#e79f00", alpha = 0.5, size = 1.5) +
    geom_point(data = summary, aes(x = level, y = mean_latent_variable), color = "#e79f00", size = 6, shape = 18) +
    geom_smooth(data = data, aes(x = level, y = latent_variable, group = 1), method = "lm", color = alpha("#0072B2", 0.6), linetype = "solid", se = FALSE) +
    #scale_color_manual(values = common_colors) +
    labs(title = title, x = "Level", y = "Latent Variable") +
    theme(
      plot.title = element_text(size = 45, hjust = 0.5),
      axis.title.x = element_text(size = 35),
      axis.text = element_text(size = 30),
      axis.title.y = element_text(size = 35),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "lightgray"),
      panel.grid.minor = element_blank()
    )
  return(plot)
}

# Plot for nypon dataset
nypon_plot <- latent_plot(nypon_data, nypon_summary, "Nypon")

# Plot for vilja dataset
vilja_plot <- latent_plot(vilja_data, vilja_summary, "Vilja")

# Plot each graph separately
print(nypon_plot)

print(vilja_plot)

ggsave("images_mc_22cov_withsuc/scatter/nypon_latent.pdf", plot = nypon_plot, device = "pdf")
ggsave("images_mc_22cov_withsuc/scatter/vilja_latent.pdf", plot = vilja_plot, device = "pdf")


#-------------------------------LIX SCATTER-------------------------


# Load required libraries
library(ggplot2)
library(dplyr)

# Read the data
data <- readRDS("stripped_data_lix.RDS")

# Define the order of levels for nypon and vilja datasets
nypon_levels <- c("1", "2", "3", "4", "5")
vilja_levels <- c("x-small", "small", "medium", "large", "x-large", "xx-large")

# Plot for nypon dataset
nypon_data <- data %>%
  filter(dataset == "nypon") %>%
  mutate(level = factor(level, levels = nypon_levels))

# Plot for vilja dataset
vilja_data <- data %>%
  filter(dataset == "vilja") %>%
  mutate(level = factor(level, levels = vilja_levels))

# Rename levels for vilja dataset
vilja_data$level <- factor(vilja_data$level,
                           levels = vilja_levels,
                           labels = c("X-S", "S", "M", "L", "X-L", "XX-L"))


# Calculate mean and standard deviation manually
nypon_summary <- nypon_data %>%
  group_by(level) %>%
  summarise(
    mean_lix_variable = mean(lixValue),
    sd_lix_variable = sd(lixValue)
  )

vilja_summary <- vilja_data %>%
  group_by(level) %>%
  summarise(
    mean_lix_variable = mean(lixValue),
    sd_lix_variable = sd(lixValue)
  )

lix_plot <- function(data, summary, title) {
  plot <- ggplot() +
    geom_point(data = data, aes(x = level, y = lixValue, color = level), size = 2, color = "#999999", show.legend = FALSE) +
    geom_errorbar(data = summary, aes(x = level, ymin = mean_lix_variable - sd_lix_variable, ymax = mean_lix_variable + sd_lix_variable), width = 0.4, color = "#e79f00", alpha = 0.5, size = 1.5) +
    geom_point(data = summary, aes(x = level, y = mean_lix_variable), color = "#e79f00", size = 6, shape = 18) +
    geom_smooth(data = data, aes(x = level, y = lixValue, group = 1), method = "lm", color = alpha("#0072B2", 0.6), linetype = "solid", se = FALSE) +
    labs(title = title, x = "Level", y = "LIX") +
    theme(
      plot.title = element_text(size = 45, hjust = 0.5),
      axis.title.x = element_text(size = 35),
      axis.text = element_text(size = 30),
      axis.title.y = element_text(size = 35),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "lightgray"),
      panel.grid.minor = element_blank()
    )
  return(plot)
}


# Plot for nypon dataset
nypon_plot <- lix_plot(nypon_data, nypon_summary, "Nypon")

# Plot for vilja dataset
vilja_plot <- lix_plot(vilja_data, vilja_summary, "Vilja")

# Plot each graph separately
print(nypon_plot)
print(vilja_plot)

ggsave("images_mc_22cov_withsuc/scatter/nypon_lix.pdf", plot = nypon_plot, device = "pdf")
ggsave("images_mc_22cov_withsuc/scatter/vilja_lix.pdf", plot = vilja_plot, device = "pdf")




#--------------------LATENT DISTRIBUTIONS---------

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Load data
predicted_data <- readRDS("aggregated_predicted_values.RDS")
latent_data <- readRDS("latent_values_all20.RDS")
str(predicted_data)
str(latent_data)
# Extract predicted values for the two observations
observation_12 <- predicted_data[[12]]
observation_204 <- predicted_data[[204]]

# Extract latent variable values for the two observations
latent_value_12 <- latent_data$latent_variable[12]
latent_value_204 <- latent_data$latent_variable[204]


plot_distribution <- function(observation_index, title=FALSE, predicted_values, latent_value, filename) {
  p <- ggplot(data.frame(value = predicted_values), aes(x = value)) +
    geom_density(fill = "#D55E00", alpha = 0.5) +
    geom_vline(xintercept = latent_value, color = "black", linetype = "solid") +
    labs(title = paste(title), x = "Latent Variable") +
    theme(plot.title = element_text(size = 45, hjust = 0.5),
          axis.text = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "lightgray"),
          panel.grid.minor = element_blank())
  
  p  # Return the plot object
}

plot_12 <- plot_distribution(12, title="Saad är byggarbetare", observation_12, latent_value_12, NULL)
plot_204 <- plot_distribution(204, title="Det går an", observation_204, latent_value_204, NULL)

print(plot_12)
print(plot_204)

ggsave("images22cov/scatter/distribution_plot_12.pdf", plot = plot_12, width = 8, height = 6, dpi = 300)
ggsave("images22cov/scatter/distribution_plot_204.pdf", plot = plot_204, width = 8, height = 6, dpi = 300)






#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#-----------------------------  MONOCHROME VARIATIONS OF SCATTER PLOTS FOLLOWS ------------------------
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------



#------------------------------------LATENT VALUE SCATTER MC ----------------------



# Load required libraries
library(ggplot2)
library(dplyr)

# Read the data
data <- readRDS("latent_values_all20.RDS")

# Define the order of levels for nypon and vilja datasets
nypon_levels <- c("1", "2", "3", "4", "5")
vilja_levels <- c("x-small", "small", "medium", "large", "x-large", "xx-large")
print(vilja_data)

# Plot for nypon dataset
nypon_data <- data %>%
  filter(dataset == "nypon") %>%
  mutate(level = factor(level, levels = nypon_levels))

# Plot for vilja dataset
vilja_data <- data %>%
  filter(dataset == "vilja") %>%
  mutate(level = factor(level, levels = vilja_levels))

# Rename levels for vilja dataset
vilja_data$level <- factor(vilja_data$level,
                           levels = vilja_levels,
                           labels = c("X-S", "S", "M", "L", "X-L", "XX-L"))

head(nypon_data)
str(vilja_data$level)

# Calculate mean and standard deviation manually
nypon_summary <- nypon_data %>%
  group_by(level) %>%
  summarise(
    mean_latent_variable = mean(latent_variable),
    sd_latent_variable = sd(latent_variable)
  )

vilja_summary <- vilja_data %>%
  group_by(level) %>%
  summarise(
    mean_latent_variable = mean(latent_variable),
    sd_latent_variable = sd(latent_variable)
  )

# Define function to create the plot
latent_plot <- function(data, summary, title) {
  plot <- ggplot() +
    geom_point(data = data, aes(x = level, y = latent_variable, color = level), size = 2, shape = 1, color = "#2D2D2D", show.legend = FALSE) +
    geom_errorbar(data = summary, aes(x = level, ymin = mean_latent_variable - sd_latent_variable, ymax = mean_latent_variable + sd_latent_variable), width = 0.4, color = "#8E9090", alpha = 0.5, size = 1.5) +
    geom_point(data = summary, aes(x = level, y = mean_latent_variable), color = "#2D2D2D", size = 6, shape = 18) +
    geom_smooth(data = data, aes(x = level, y = latent_variable, group = 1), method = "lm", color = alpha("black", 0.6), linetype = "solid", se = FALSE) +
    #scale_color_manual(values = common_colors) +
    labs(title = title, x = "Level", y = "Latent Variable") +
    theme(
      plot.title = element_text(size = 45, hjust = 0.5),
      axis.title.x = element_text(size = 35),
      axis.text = element_text(size = 30),
      axis.title.y = element_text(size = 35),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "lightgray"),
      panel.grid.minor = element_blank()
    )
  return(plot)
}

# Plot for nypon dataset
nypon_plot <- latent_plot(nypon_data, nypon_summary, "Nypon")

# Plot for vilja dataset
vilja_plot <- latent_plot(vilja_data, vilja_summary, "Vilja")

# Plot each graph separately
print(nypon_plot)

print(vilja_plot)

ggsave("images_mc_22cov_withsuc/scatter/nypon_latent.pdf", plot = nypon_plot, device = "pdf")
ggsave("images_mc_22cov_withsuc/scatter/vilja_latent.pdf", plot = vilja_plot, device = "pdf")


#-------------------------------LIX SCATTER MC-------------------------


# Load required libraries
library(ggplot2)
library(dplyr)

# Read the data
data <- readRDS("stripped_data_lix.RDS")

# Define the order of levels for nypon and vilja datasets
nypon_levels <- c("1", "2", "3", "4", "5")
vilja_levels <- c("x-small", "small", "medium", "large", "x-large", "xx-large")

# Plot for nypon dataset
nypon_data <- data %>%
  filter(dataset == "nypon") %>%
  mutate(level = factor(level, levels = nypon_levels))

# Plot for vilja dataset
vilja_data <- data %>%
  filter(dataset == "vilja") %>%
  mutate(level = factor(level, levels = vilja_levels))

# Rename levels for vilja dataset
vilja_data$level <- factor(vilja_data$level,
                           levels = vilja_levels,
                           labels = c("X-S", "S", "M", "L", "X-L", "XX-L"))


# Calculate mean and standard deviation manually
nypon_summary <- nypon_data %>%
  group_by(level) %>%
  summarise(
    mean_lix_variable = mean(lixValue),
    sd_lix_variable = sd(lixValue)
  )

vilja_summary <- vilja_data %>%
  group_by(level) %>%
  summarise(
    mean_lix_variable = mean(lixValue),
    sd_lix_variable = sd(lixValue)
  )

lix_plot <- function(data, summary, title) {
  plot <- ggplot() +
    geom_point(data = data, aes(x = level, y = lixValue, color = level), size = 2, shape = 1, color = "#2D2D2D", show.legend = FALSE) +
    geom_errorbar(data = summary, aes(x = level, ymin = mean_lix_variable - sd_lix_variable, ymax = mean_lix_variable + sd_lix_variable), width = 0.4, color = "#8E9090", alpha = 0.5, size = 1.5) +
    geom_point(data = summary, aes(x = level, y = mean_lix_variable), color = "#2D2D2D", size = 6, shape = 18) +
    geom_smooth(data = data, aes(x = level, y = lixValue, group = 1), method = "lm", color = alpha("black", 0.6), linetype = "solid", se = FALSE) +
    labs(title = title, x = "Level", y = "LIX") +
    theme(
      plot.title = element_text(size = 45, hjust = 0.5),
      axis.title.x = element_text(size = 35),
      axis.text = element_text(size = 30),
      axis.title.y = element_text(size = 35),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "lightgray"),
      panel.grid.minor = element_blank()
    )
  return(plot)
}


# Plot for nypon dataset
nypon_plot <- lix_plot(nypon_data, nypon_summary, "Nypon")

# Plot for vilja dataset
vilja_plot <- lix_plot(vilja_data, vilja_summary, "Vilja")

# Plot each graph separately
print(nypon_plot)
print(vilja_plot)

ggsave("images_mc_22cov_withsuc/scatter/nypon_lix.pdf", plot = nypon_plot, device = "pdf")
ggsave("images_mc_22cov_withsuc/scatter/vilja_lix.pdf", plot = vilja_plot, device = "pdf")




#--------------------LATENT DISTRIBUTIONS MC---------

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Load data
predicted_data <- readRDS("aggregated_predicted_values.RDS")
latent_data <- readRDS("latent_values_all20.RDS")
str(predicted_data)
str(latent_data)
# Extract predicted values for the two observations
observation_12 <- predicted_data[[12]]
observation_204 <- predicted_data[[204]]

# Extract latent variable values for the two observations
latent_value_12 <- latent_data$latent_variable[12]
latent_value_204 <- latent_data$latent_variable[204]


plot_distribution <- function(observation_index, title=FALSE, predicted_values, latent_value, filename) {
  p <- ggplot(data.frame(value = predicted_values), aes(x = value)) +
    geom_density(fill = "white", alpha = 0.5) +
    geom_vline(xintercept = latent_value, color = "black", linetype = "solid") +
    labs(title = paste(title), x = "Latent Variable") +
    theme(plot.title = element_text(size = 45, hjust = 0.5),
          axis.text = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "lightgray"),
          panel.grid.minor = element_blank())
  
  p  # Return the plot object
}

plot_12 <- plot_distribution(12, title="Saad är byggarbetare", observation_12, latent_value_12, NULL)
plot_204 <- plot_distribution(204, title="Det går an", observation_204, latent_value_204, NULL)

print(plot_12)
print(plot_204)

ggsave("images_mc_22cov_withsuc/scatter/distribution_plot_12.pdf", plot = plot_12, width = 8, height = 6, dpi = 300)
ggsave("images_mc_22cov_withsuc/scatter/distribution_plot_204.pdf", plot = plot_204, width = 8, height = 6, dpi = 300)

