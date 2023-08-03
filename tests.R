make_laplacian_eigenmaps_test <- function(theta){
  r <- theta$r
  K <- theta$K
  laplacian_eigenmaps_test <- function(Y,X){
    # 1. Compute inner products.
    rs <- names(precomputed_data$inner_products)
    a_ks <- precomputed_data$inner_products[[which(names(precomputed_data$inner_products) == r)]]
    a_K <- a_ks[hh,1:K]
    
    # 2. Compute L2 norm.
    norm <- sum(a_K^2)
    
    # 3. See whether norm is greater than threshold
    norm > threshold
  }
}
attr(make_laplacian_eigenmaps_test,"precompute_data") <- function(Y,X,theta_df){
  # Precompute eigenvectors
  precomputed_data <- precompute_laplacian_eigenmaps_data(Y,X,theta_df)
  
  precomputed_inner_products <- list()
  for(kk in 1:length(precomputed_data$train))
  {
    V <- precomputed_data$train[[kk]]
    Y_mat <- do.call(rbind,Y)
    precomputed_inner_products[[kk]] <- Y_mat %*% V
  }
  names(precomputed_inner_products) <- names(precomputed_data$train)
  precomputed_data$train <- list(eigenvectors = precomputed_data$train,
                                 inner_products = precomputed_inner_products)
  return(precomputed_data)
}
attr(make_laplacian_eigenmaps_test, "compute_threshold") <- function(X,theta_df,precomputed_data){
  threshold <- numeric(nrow(theta_df))
  norm_null <- matrix(nrow = nrow(theta_df),ncol = B)
  Y_null <- matrix(rnorm(nrow(X)*B,0,1),nrow = nrow(X),ncol = B)
  rs <- unique(theta_df$r)
  for(kk in 1:length(rs))
  {
    r <- rs[kk]
    K_indx <- which(theta_df$r == r)
    V <- precomputed_data$train$eigenvectors[[which(names(precomputed_data$train$eigenvectors) == r)]]
    a_ks <- t(Y_null) %*% V
    norm_null <- apply(a_ks,1,FUN = function(a_k){cumsum(a_k^2)})
    threshold[K_indx] <- apply(norm_null,1,FUN = function(norms){quantile(norms,1 - alpha)})
  }
  return(threshold)
}

make_spectral_projection_test <- function(theta){
  K <- theta$K
  least_squares_test <- function(Y,X){
    # 1. Compute inner products.
    V_K <- precomputed_data[,1:K,drop = F]
    a_K <- 1/n * t(Y) %*% V_K
    
    # 2. Compute L2 norm.
    norm <- sum(a_K^2)
    
    # 3. See whether norm is greater than threshold
    norm > threshold
  }
}
attr(make_spectral_projection_test,"precompute_data") <- function(Y,X,theta_df){
  precompute_least_squares_data(Y,X,theta_df)
}
attr(make_spectral_projection_test, "compute_threshold") <- function(X,theta_df,precomputed_data){
  norm_null <- matrix(nrow = nrow(theta_df),ncol = B)
  Y_null <- matrix(rnorm(nrow(X)*B,0,1),nrow = nrow(X),ncol = B)

  V <- precomputed_data$train
  a_ks <- 1/n * t(Y_null) %*% V
  norm_null <- apply(a_ks,1,FUN = function(a_k){cumsum(a_k^2)})
  threshold <- apply(norm_null,1,FUN = function(norms){quantile(norms,1 - alpha)})
  return(threshold)
}