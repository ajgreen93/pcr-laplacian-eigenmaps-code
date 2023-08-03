#--------------------------------------------------#
# Script to generate plots analyzing impact of LE tuning parameters.
#--------------------------------------------------#
library(dplyr)
library(ggplot2)
library(gridExtra)
library(viridis)
source("plot_methods.R")
source("misc.R")

# User entered information.
data_directory <- "data/tuning/sobolev_1s_1d"
plot_directory <- "plots/tuning/sobolev" # Please change this to whichever directory you prefer.
if(!exists(plot_directory)){dir.create(plot_directory,recursive = T)}

# Load data
load(file.path(data_directory,"configs.R"))
d <- configs$d; s <- configs$s; ns <- configs$ns; methods <- configs$methods
load(file.path(data_directory,"best_fits_by_method.R"))
load(file.path(data_directory,"thetas.R"))
load(file.path(data_directory,"mse.R"))
if("test_mse.R" %in% list.files(data_directory)) load(file.path(data_directory,"test_mse.R"))
load(file.path(data_directory,"Xs.R"))
load(file.path(data_directory,"f0s.R"))
load(file.path(data_directory,"Ys.R"))

# Subset data to methods you actually want to plot.
plot_n <- 1000
plot_methods <- c("laplacian_eigenmaps",
                  "spectral_projection")
stopifnot(plot_n %in% ns)
stopifnot(all(plot_methods %in% names(methods)))
mse <- mse[ns %in% plot_n][[1]][names(methods) %in% plot_methods]
if(exists("test_mse")) test_mse <- test_mse[ns %in% plot_n][[1]][names(methods) %in% plot_methods]
thetas <- thetas[ns %in% plot_n][[1]][names(methods) %in% plot_methods]
methods <- methods[names(methods) %in% plot_methods]

# Plotting parameters for all plots
# Capitalization function from https://rstudio-pubs-static.s3.amazonaws.com/408658_512da947714740b99253228f084a08a9.html.
CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}
function_names <- c("eigenfunction","sobolev")
function_name <- function_names[which(sapply(function_names,FUN = function(name){grepl(name,data_directory)}))]
if(length(function_name) == 0) function_name = ""
title <- paste0("d = ",d,", s = ",s,". ", CapStr(function_name),".")

## Plot 1: Mean squared error as a function of K.

# Find best parameters
best_parameters <- find_best_parameters(list(mse),list(thetas))[[1]]
alg_indx <- sapply(thetas,FUN = function(theta){"K" %in% names(theta)}) %>% which()
Ks <- sapply(thetas[alg_indx],FUN = function(theta){unique(theta$K)})

# Mse as a function of K, for best other parameters.
plot_mse <- matrix(nrow = nrow(Ks),ncol = ncol(Ks))
colnames(plot_mse) <- names(methods)[alg_indx]
for(ii in alg_indx)
{
  best_parameters_ii <- best_parameters[[ii]]
  mse_ii <- mse[[ii]]
  thetas_ii <- thetas[[ii]]
  if(exists("test_mse")) test_mse_ii <- test_mse[[ii]]
  if(is.null(names(best_parameters_ii))){
    # Spectral projection or least squares
    plot_mse[,ii] <- rowMeans(mse_ii)
  } else if(all(names(best_parameters_ii) == c("r","K"))){
    # Laplacian eigenmaps
    plot_mse[,ii] <- rowMeans(mse_ii[thetas_ii$r == best_parameters_ii$r,])
  } else if(all(names(best_parameters_ii) == c("r","K","h"))){
    # Laplacian eigenmaps plus kernel smoothing
    plot_mse[,ii] <- rowMeans(test_mse_ii[thetas_ii$r == best_parameters_ii$r & thetas_ii$h == best_parameters_ii$h,])
  }
}

# Plotting parameters
xlims <- c(min(Ks),max(Ks))
ylims <- c(min(plot_mse),max(plot_mse))
cols <- c("red","green")
stopifnot(ncol(plot_mse) <= 2)

plot_name <- paste0("mse_by_number_of_eigenvectors_",plot_n,"n_",d,"d_",s,"s.pdf")
pdf(file.path(plot_directory,plot_name))
par(mar = c(5.1,6,4.1,4.1)) # Prevent the left hand side from getting cut off.
# shadow plot
plot(x = Ks[,1], xlim = xlims, ylim = ylims, xlab = "Number of eigenvectors", ylab = "Mean squared error", 
     main = title,
     cex.main = 2.5, cex.lab = 2.5, cex.axis = 2, type = "n")

# add points and lines
for(jj in 1:ncol(plot_mse))
{
  points(x = Ks[,jj], y =  plot_mse[,jj], col = cols[jj], pch = 20)
  lines(x = Ks[,jj], y =  plot_mse[,jj], col = cols[jj], lwd = 1.5)
}
grid(lwd = 2)
dev.off()

## Plot 2: Mean squared error as a function of radius.
alg_indx_r <- sapply(thetas,FUN = function(theta){"r" %in% names(theta)}) %>% which()
rs <- sapply(thetas[alg_indx_r],FUN = function(theta){unique(theta$r)})

# Mse as a function of r, for best other parameters.
plot_mse <- matrix(nrow = nrow(rs),ncol = ncol(rs))
for(jj in 1:length(alg_indx_r))
{
  ii <- alg_indx_r[jj]
  best_parameters_ii <- best_parameters[[ii]]
  mse_ii <- mse[[ii]]
  thetas_ii <- thetas[[ii]]
  if(exists("test_mse")) test_mse_ii <- test_mse[[ii]]
  if(all(names(best_parameters_ii) == c("r","K"))){
    # Laplacian eigenmaps
    plot_mse[,jj] <- unique(rowMeans(mse_ii[thetas_ii$K == best_parameters_ii$K,]))
  } else if(all(names(best_parameters_ii) == c("r","K","h"))){
    # Laplacian eigenmaps plus kernel smoothing
    plot_mse[,jj] <- rowMeans(test_mse_ii[thetas_ii$K == best_parameters_ii$K & thetas_ii$h == best_parameters_ii$h,])
  }
}

# Plotting parameters
xlims <- c(min(rs),max(rs))
ylims <- c(min(plot_mse),max(plot_mse))
cols <- c("red")
stopifnot(ncol(plot_mse) <= 1)

plot_name <- paste0("mse_by_radius_",plot_n,"n_",d,"d_",s,"s.pdf")
pdf(file.path(plot_directory,plot_name))
par(mar = c(5.1,6,4.1,4.1)) # Prevent the left hand side from getting cut off.
# shadow plot
plot(x = rs[,1], xlim = xlims, ylim = ylims, xlab = "Radius", ylab = "Mean Squared Error", 
     main = title,
     cex.main = 2.5, cex.lab = 2.5, cex.axis = 2, type = "n")

# add points and lines
for(jj in 1:ncol(plot_mse))
{
  points(x = rs[,jj], y =  plot_mse[,jj], col = cols[jj], pch = 20)
  lines(x = rs[,jj], y =  plot_mse[,jj], col = cols[jj], lwd = 1.5)
}
grid(lwd = 2)
dev.off()