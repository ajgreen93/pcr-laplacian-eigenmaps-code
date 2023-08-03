library(RANN)
library(dplyr)
library(Matrix)
library(reshape2)
library(RSpectra)
source("precompute_data.R")
source("sample.R")
source("graph.R")
source("misc.R")
source("initialize_parameters.R")
source("tests.R")

# Please change this line to whichever config you wish to run.
source("configs/laplacian_eigenmaps/testing/eigenfunction_1d_1s.R")

# Options for running the pipeline
verbose <- T

#----------------------------------------------------#
# This is the pipeline for running regression estimation experiments.
#----------------------------------------------------#

# For output.
power <- vector(mode = "list",length = length(ns))
err <- vector(mode = "list",length = length(ns))
Xs <- vector(mode = "list", length = length(ns))
f0s <- vector(mode = "list", length = length(ns))
Ys <- vector(mode = "list", length = length(ns))
thetas <- vector(mode = "list", length = length(ns))

for(ii in 1:length(ns))
{
  n <- ns[ii]
  f0s <- make_f0s(d,n)
  power[[ii]] <- vector(mode = "list", length = length(methods))
  err[[ii]] <- vector(mode = "list", length = length(methods))
  
  # Initialize tuning parameters (thetas).
  thetas[[ii]] <- vector(mode = "list", length(methods))
  
  for(jj in 1:length(methods))
  {
    thetas[[ii]][[jj]] <- initialize_thetas[[jj]](sample_X, n)
    err[[ii]][[jj]] <- array(dim = c(nrow(thetas[[ii]][[jj]]),length(f0s),iters))
    power[[ii]][[jj]] <- matrix(nrow = nrow(thetas[[ii]][[jj]]), ncol = length(f0s))
  }
    
  for(iter in 1:iters)
  {
    ### Data. ###
    X <- sample_X(n)
    Y_list <- vector(mode = "list",length = length(f0s))
    
    for(hh in 1:length(f0s))
    {
      f0 <- f0s[[hh]]
      
      ### Response data. ###
      f0_evaluations <- apply(X,1,f0)
      Y_list[[hh]] <- f0_evaluations + rnorm(n,0,1)
    }
    
    ### Precomputed data. ###
    precomputed_data_list <- list() 
    threshold_list <- list()
    for(jj in 1:length(methods))
    {
      method <- methods[[jj]]
      
      # Do any computations here that we would have to repeat for each hyperparameter (theta)
      # E.g. compute eigenvectors.
      precompute_data <- attributes(method)$precompute_data
      if(!is.null(precompute_data)){
        precomputed_data_list[[jj]] <- precompute_data(Y_list,X,thetas[[ii]][[jj]])
      } else{
        precomputed_data_list[[jj]] <- NULL
      }
      
      # Compute thresholds
      compute_threshold <- attributes(method)$compute_threshold
      threshold_list[[jj]] <- compute_threshold(X,thetas[[ii]][[jj]],precomputed_data_list[[jj]])
    }
    
    ### Cycle through f0s, perform tests. ###
    for(hh in 1:length(f0s))
    {
      Y <- Y_list[[hh]]
      
      ### Analysis. ###
      for(jj in 1:length(methods))
      {
        method <- methods[[jj]]
        precomputed_data <- precomputed_data_list[[jj]]
        threshold <- threshold_list[[jj]]
        
        for(kk in 1:nrow(thetas[[ii]][[jj]])){
          theta <- slice(thetas[[ii]][[jj]],kk)
          test <- method(theta)
          environment(test)$precomputed_data <- precomputed_data[["train"]]
          environment(test)$threshold <- threshold[kk]
          err[[ii]][[jj]][kk,hh,iter] <- test(Y,X)
        }
        if(verbose){
          logger::log_info("n: ", n, ".",
                           "Iter: ", iter," out of ", iters, ".",
                           "Method: ", jj, " out of ", length(methods), ".",
                           "Function: ",hh, " out of ", length(f0s),".")
        }
      }
    }
    logger::log_info("n: ", n, ".",
                     "Iter: ", iter," out of ", iters, ".")
  }
  
  ### Compute power. ###
  for(jj in 1:length(methods))
  {
    power[[ii]][[jj]] <- rowMeans(err[[ii]][[jj]], dims = 2)
  }
  
  # Taking the last iteration X, it doesn't matter. 
  Xs[[ii]] <- X
  f0s[[ii]] <- f0
  Ys[[ii]] <- Y
}

### Save. ###
save_directory <- file.path("data",gsub("[^[:alnum:]]", "", Sys.time()))
for(directory in c(save_directory))
{
  dir.create(directory)
}
configs <- list(d = d,
                s = s,
                ns = ns,
                M = M,
                methods = methods,
                sample_X = sample_X,
                make_f0s = make_f0s)
save(configs, file = paste0(save_directory, '/configs.R'))
save(thetas, file = paste0(save_directory, '/thetas.R'))
save(err, file = paste0(save_directory, '/err.R'))
save(power, file = paste0(save_directory, '/power.R'))