#----------------------------------------------------#
# Functions to precompute data which will be used repeatedly.
# E.g. graph eigenvectors and trigonometric basis.
#----------------------------------------------------#

precompute_least_squares_data <- function(Y,X,theta_df){
  precomputed_data <- vector(mode = "list",length = 1)
  names(precomputed_data) <- c("train")
  
  K <- max(theta_df$K)
  trig_basis <- get_trigonometric_basis(d,K) # Fourier basis
  precomputed_data[["train"]] <- 
    if(K == 1) apply(X,1,trig_basis) else apply(X,1,trig_basis) %>% t()
  return(precomputed_data)
}

precompute_laplacian_eigenmaps_data <- function(Y,X,theta_df){
  if(!is.list(Y)) Y  <- list(Y)
  rs <- unique(theta_df$r)
  precomputed_data <- vector(mode = "list",length = length(rs))
  names(precomputed_data) <- rs
  for(kk in 1:length(rs)){
    r <- rs[kk]
    K <- thetas[[ii]][[jj]] %>% filter(r == !!r) %>% pull(K)
    
    # Build G_{n,r} over X.
    G <- neighborhood_graph(X,r)
    
    # Get L.
    L <- Laplacian(G)
    
    # Compute as many eigenvectors as we will need.
    spectra <- get_spectra(L,max(K))
    L_eigenvectors <- spectra$vectors[,ncol(spectra$vectors):1]
    
    # Store eigenvectors
    precomputed_data[[kk]] <- L_eigenvectors 
  }
  precomputed_data <- list(train = precomputed_data)
  return(precomputed_data)
}