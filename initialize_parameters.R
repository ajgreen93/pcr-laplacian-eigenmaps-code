#------------------------------------------------------#
# Functions to choose initial grid of tuning parameters, for estimation.
#------------------------------------------------------#
initialize_laplacian_eigenmaps_thetas <- function(sample_X,n,
                                                  n_rs = 10,
                                                  n_Ks = max_K,
                                                  max_degree = 60,
                                                  max_K = max(initialize_spectral_projection_thetas(sample_X,n)$K),
                                                  dist_to_r = max){
  
  degrees <- seq(2,max_degree + 1,length.out = n_rs) %>% round() %>% unique()
  
  # Choose parameters in a data-dependent way using Monte Carlo
  iters <- 10 
  
  # Choose a range of radii based on different desired degree.
  degree_dependent_rs <- matrix(ncol = n_rs, nrow = iters)
  for(iter in 1:iters)
  {
    X <- sample_X(n)
    for(jj in 1:n_rs)
    {
      rneighborhood <- nn2(data = X, query = X, k = degrees[jj]) 
      degree_dependent_rs[iter,jj] <- dist_to_r(rneighborhood$nn.dists[,degrees[jj]])
    }
  }
  # Choose a range of radii intended to achieve the desired min degree.
  rs <- colMeans(degree_dependent_rs) %>% sort() %>% unique()
  
  # Choose same range of eigenvectors as used for spectral projection
  if(is.null(max_K)){
    Ks <- initialize_spectral_projection_thetas(sample_X,n) %>% pull(K)
  } else{
    Ks <- initialize_spectral_projection_thetas(sample_X,n,max_K) %>% pull(K)
  }
  thetas <- expand.grid(r = rs, K = Ks)
}

initialize_spectral_projection_thetas <- function(sample_X,n,max_K = n/50){
  thetas <- data.frame(K = 1:max_K) 
}

#------------------------------------------------------------------#
# Functions to choose initial grid of tuning parameters, for testing.
#------------------------------------------------------------------#
initialize_spectral_projection_test_thetas <- function(sample_X,n){
  max_K <- round(2*(M^2*n)^(2*d/(4*s + d)))
  initialize_spectral_projection_thetas(sample_X,n,max_K)
}

initialize_laplacian_eigenmaps_test_thetas <- function(sample_X,n){
  max_K <- round(2*(M^2*n)^(2*d/(4*s + d)))
  initialize_laplacian_eigenmaps_thetas(sample_X,n,max_K)
}