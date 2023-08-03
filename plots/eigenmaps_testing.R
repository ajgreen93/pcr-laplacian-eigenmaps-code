#--------------------------------------------------#
# Script to generate plots used for Figure 2 of Laplacian Eigenmaps paper.
#--------------------------------------------------#
library(dplyr)
library(ggplot2)
library(gridExtra)
library(viridis)
source("misc.R")
source("plot_methods.R")
source("sample.R")

# User entered information.
data_directory <- "data/testing/eigenfunction_2d_2s"
plot_directory <- "plots/testing/eigenfunction" # Please change this to whichever directory you prefer.
if(!exists(plot_directory)){dir.create(plot_directory,recursive = T)}

# Plotting parameters
cols <- c("red","green")
sd <- T

# Load data
load(file.path(data_directory,"configs.R"))
d <- configs$d; ns <- configs$ns; s <- configs$s; M <- configs$M; make_f0s <- configs$make_f0s
methods <- configs$methods
load(file.path(data_directory,"thetas.R"))
load(file.path(data_directory,"err.R"))
load(file.path(data_directory,"power.R"))

# Subset data to methods you actually want to plot.
plot_methods <- c("laplacian_eigenmaps",
                  "spectral_projection")
err <- lapply(err,FUN = function(m){m[names(methods) %in% plot_methods]})
power <- lapply(power,FUN = function(m){m[names(methods) %in% plot_methods]})
methods <- methods[names(methods) %in% plot_methods]

# Plot critical radius.
plot_name <- paste0("critical_radius_by_sample_size_",d,"d_",s,"s.pdf")
pdf(file.path(plot_directory,plot_name))
par(mar = c(5.1,6,4.1,6),cex.main = 2.5)
plot_testing_critical_radius(methods,err,type_II_error = runif(100,0,1),
                             sd, cols, verbose = T)
dev.off()