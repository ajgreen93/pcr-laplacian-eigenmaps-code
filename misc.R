#------------------------------------------#
# Helper functions
#------------------------------------------#

get_trigonometric_basis <- function(d,K){
  # Get the indices for all tensor product polynomials.
  max_indx <- ceiling( sqrt(d*K^(2/d)) )
  indx <- expand.grid(rep(list(0:(max_indx - 1)), d)) %>% as.matrix()
  eigenvalues <- apply(indx,1,FUN = function(k){sum(k^2)})
  indx <- indx[order(eigenvalues),,drop = FALSE][1:K,,drop = FALSE]
  
  # tensor produt trigonometric polynomials
  trig_basis <- function(x){
    basis_by_dim <- sqrt(2) * t(cos(t(indx)*pi*(x + 1)/2)) * (1/sqrt(2))^(indx == 0)
    basis <- rep(1,K)
    for(ii in 1:d){
      basis <- basis * basis_by_dim[,ii]
    }
    return(basis)
  }
  return(trig_basis)
}

get_eigenfunction <- function(d,K){
  max_indx <- ceiling( sqrt(d*K^(2/d)) )
  indx <- expand.grid(rep(list(0:(max_indx - 1)), d))
  eigenvalues <- apply(indx,1,FUN = function(k){sum(k^2)})
  k <- indx[order(eigenvalues)[K],] %>% as.numeric()
  eigenfunction <- function(x){1/sqrt(2)^(sum(k == 0)) * prod(sqrt(2) * cos(k*pi*(x + 1)/2))}
  return(eigenfunction)
}

get_eigenvalue <- function(d,K){
  max_indx <- ceiling( sqrt(d*K^(2/d)) )
  indx <- expand.grid(rep(list(0:(max_indx - 1)), d))
  eigenvalues <- apply(indx,1,FUN = function(k){sum(k^2)})
  k <- indx[order(eigenvalues)[K],,drop = FALSE]
  eigenvalue <- sum(k^2)
  return(eigenvalue)
}

find_best_parameters <- function(mse,thetas){
  best_params <- vector(mode = "list",length = length(mse))
  for(ii in 1:length(mse))
  {
    best_params[[ii]] <- vector(mode = "list",length = length(thetas[[1]]))
    for(jj in 1:length(mse[[ii]]))
    {
      best_params[[ii]][[jj]] <- thetas[[ii]][[jj]][which.min(rowSums(mse[[ii]][[jj]])),]
    }
  }
  best_params
}

find_best_K <- function(mse,thetas){
  best_K <- matrix(nrow = length(ns),ncol = length(thetas[[1]]))
  for(ii in 1:length(ns))
  {
    for(jj in 1:length(mse[[ii]]))
    {
      best_K[ii,jj] <- thetas[[ii]][[jj]][which.min(rowSums(mse[[ii]][[jj]])),"K"]
    }
  }
  best_K
}

find_best_radius <- function(mse,thetas){
  best_r <- matrix(nrow = length(ns),ncol = length(thetas[[1]]))
  for(ii in 1:length(ns))
  {
    for(jj in 1:length(mse[[ii]]))
    {
      if(!("r" %in% names(thetas[[ii]][[jj]]))){next}
      best_r[ii,jj] <- thetas[[ii]][[jj]][which.min(rowSums(mse[[ii]][[jj]])),"r"]
    }
  }
  best_r
}

find_best_mse <- function(mse,methods){
  best_mse <- vector(mode = "list", length = length(methods))
  names(best_mse) <- names(methods)

  for(jj in 1:length(methods))
  {
    best_mse_jj <- numeric()
    sd_best_mse_jj <- numeric()
    for(ii in 1:length(ns))
    {
      mse_ii_jj <- mse[[ii]][[jj]]
      best_mse_jj[ii] <- min(rowMeans(mse_ii_jj))
      sd_best_mse_jj[ii] <- apply(mse_ii_jj,1,sd)[which.min(rowMeans(mse_ii_jj))]/sqrt(ncol(mse_ii_jj))
    }
    best_mse[[jj]] <- data.frame(x = ns, y = best_mse_jj,sd = sd_best_mse_jj)
  }
  return(best_mse)
}

find_critical_radius <- function(err,methods, make_f0s, type_II_error = .5, sd = T, verbose = F){
  # Compute power
  power <- list()
  for(ii in 1:length(err))
  {
    power[[ii]] <- list()
    for(jj in 1:length(methods))
    {
      power[[ii]][[jj]] <- rowMeans(err[[ii]][[jj]], dims = 2)
    }
  }

  # Compute critical radius
  critical_radius <- list()
  if(sd) iters <- dim(err[[1]][[1]])[3] # Always the same number of iterations...
  for(jj in 1:length(methods))
  {
    if(sd){
      critical_radius[[jj]] <- matrix(nrow = length(power),ncol = 3)
      colnames(critical_radius[[jj]]) <- c("crit_radius","crit_radius_lb","crit_radius_ub")
    } else{
      critical_radius[[jj]] <- matrix(nrow = length(power),ncol = 1)
      colnames(critical_radius[[jj]]) <- c("crit_radius")
    }
    for(ii in 1:length(power))
    {
      f0s <- make_f0s(d, ns[[ii]]) # First function is the zero function, so uninteresting.
      l2_norm <- sapply(f0s,FUN = function(f){attributes(f)$l2_norm})

      # Find the easiest function for which we have at least 1 - type_II_error power
      critical_radius_idx <- numeric(length(type_II_error))
      critical_radius_lb_idx <- numeric(length(type_II_error))
      critical_radius_ub_idx <- numeric(length(type_II_error))
      for(tt in 1:length(type_II_error))
      {
        critical_radius_idx[tt] <- apply(power[[ii]][[jj]],2,
                                     FUN = function(p){any(p >= 1 - type_II_error[tt])}) %>%
          which() %>% max()

        if(sd){
          critical_radius_lb_idx[tt] <- apply(power[[ii]][[jj]],2,
                                          FUN = function(p){any(p + sqrt(p*(1 - p)/iters)  >= 1 - type_II_error[tt])}) %>%
            which() %>% max()
          critical_radius_ub_idx[tt] <- apply(power[[ii]][[jj]],2,
                                          FUN = function(p){any(p - sqrt(p*(1 - p)/iters)  >= 1 - type_II_error[tt])}) %>%
            which() %>% max()
        }
      }

      if(mean(critical_radius_idx) == ncol(power[[ii]][[jj]])){
        warning("Didn't achieve critical radius for n = ",ns[[ii]],"!")
      }
      if(verbose) print(critical_radius_idx)

      if(sd){
        critical_radius[[jj]][ii,] <- c(mean(l2_norm[critical_radius_idx]),
                                        mean(l2_norm[critical_radius_lb_idx]),
                                        mean(l2_norm[critical_radius_ub_idx]))
      } else{
        critical_radius[[jj]][ii,] <- l2_norm[critical_radius_idx]
      }
    }
  }
  critical_radius
}

# Compute the L2 norm of eigenfunction.
get_l2_norm <- function(g){
  stopifnot(attributes(g)$type == "eigenfunction")
  g_env <- environment(g)
  M <- g_env$M
  s <- g_env$s
  lambda <- g_env$lambda_K
  (M/lambda^(s/2))^2
}

get_bias_spectral_sobolev <- function(Ks){
  beta <- environment(make_f0(d,n))$beta
  sapply(Ks,function(K){sum(beta[-(1:K)]^2)})
}

theory_K_spectral_sobolev <- function(ns){
  max_K <- length(environment(make_f0(d,n))$beta)
  temp <- matrix(nrow = max_K,ncol = length(ns))
  for(ii in 1:length(ns))
  {
    temp[,ii] <- (1:max_K)/ns[ii] + get_bias_spectral_sobolev(1:max_K)
  }
  apply(temp,2,which.min)
}

mse_slope <- function(mse,ns){
  log_mse <- log(mse)
  log_ns <- log(ns)
  lm(log_mse ~ log_ns)$coefficients[2]
}

# Build unlabeled-labeled splits for graph Laplacian smoothing.
#
# Two twists to the usual notion of building splits.
# (1): We ensure that no connected component belongs entirely to the same fold.
# (2): If the graph has any isolets, those vertices are always treated as
# labeled; i.e. they belong to *all* folds.
# These twists ensure that the Laplacian smoothing problem will be well-posed downstream.
build_splits <- function(G,nfolds)
{
  n <- nrow(G)

  # Compute connected components and isolets
  components <- igraph::graph_from_adjacency_matrix(G, mode = "undirected") %>%
    igraph::components()
  isolets <- (1:n)[components$membership %in% which(components$csize == 1)]

  # Find folds, respecting connectivity structure as discussed above.
  folds <- numeric(n)
  folds[order(components$membership)] <- rep(1:nfolds, n)[1:n]
  unlabeled = labeled = vector(mode="list", length=nfolds)
  for(jj in 1:nfolds) {
    unlabeled[[jj]] = union(which(folds == jj),isolets)
    labeled[[jj]] = union(which(folds != jj),isolets)
  }
  splits <- list(unlabeled = unlabeled,
                 labeled = labeled)
  return(splits)
}