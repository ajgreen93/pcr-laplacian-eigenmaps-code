make_laplacian_eigenmaps <- function(theta)
{
  r <- theta[["r"]]
  K <- theta[["K"]]
  laplacian_eigenmaps <- function(Y,X){
    # If we've already computed the eigenvectors don't do it again
    if(exists("precomputed_data")){
      stopifnot(r %in% names(precomputed_data))
      L_eigenvectors <- precomputed_data[[which(names(precomputed_data) == r)]]
    } else{
      # Build G_{n,r} over X.
      G <- neighborhood_graph(X,r)
      
      # Get L.
      L <- Laplacian(G)
      
      # Get spectra.
      spectra <- get_spectra(L,K)
      L_eigenvectors <- spectra$vectors[,ncol(spectra$vectors):1]
    }
    V_K <- L_eigenvectors[,1:K] # Eigenvectors we need
    a_k <- t(V_K) %*% Y         # Empirical Fourier coefficients
    f_hat <- V_K %*% a_k        # Spectral projection
    
    # "Prediction" function
    predict <- function(X){f_hat}
    return(predict)
  }
}
attr(make_laplacian_eigenmaps,"precompute_data") <- precompute_laplacian_eigenmaps_data

make_spectral_projection <- function(theta)
{
  K <- theta[["K"]]
  spectral_projection <- function(Y,X)
  {
    d <- ncol(X)
    n <- length(Y)
    
    if(exists("precomputed_data")){
      psi_K <- precomputed_data[,1:K,drop = FALSE]
    } else{
      trig_basis <- get_trigonometric_basis(d,K) # Fourier basis
      psi_K <- if(K == 1) apply(X,1,trig_basis) else apply(X,1,trig_basis) %>% t()
    }
    a_k <- 1/n * (t(psi_K) %*% Y)   # Empirical Fourier coefficients
    
    # Prediction function.
    predict <- function(X){
      
      if(exists("precomputed_data")){
        psi_K <- precomputed_data[,1:K,drop = FALSE]
      } else{
        trig_basis <- get_trigonometric_basis(d,K) # Fourier basis
        psi_K <- if(K == 1) apply(X,1,trig_basis) else apply(X,1,trig_basis) %>% t()
      }
      
      preds <- psi_K %*% a_k # Spectral projection
      return(preds)
    }
    
    return(predict)
  }
} 
attr(make_spectral_projection,"precompute_data") <- precompute_least_squares_data

make_least_squares <- function(theta)
{
  K <- theta[["K"]]
  least_squares <- function(Y,X)
  {
    d <- ncol(X)
    n <- length(Y)
    
    if(exists("precomputed_data")){
      psi_K <- precomputed_data[,1:K,drop = FALSE]
    } else{
      trig_basis <- get_trigonometric_basis(d,K) # Fourier basis
      psi_K <- if(K == 1) apply(X,1,trig_basis) else apply(X,1,trig_basis) %>% t()
    }
    
    a_k <- lm.fit(x = psi_K,y = Y)$coefficients %>% as.matrix()
    
    # Prediction function.
    predict <- function(X){
      
      if(exists("precomputed_data")){
        psi_K <- precomputed_data[,1:K,drop = FALSE]
      } else{
        trig_basis <- get_trigonometric_basis(d,K) # Fourier basis
        psi_K <- if(K == 1) apply(X,1,trig_basis) else apply(X,1,trig_basis) %>% t()
      }
      
      preds <- psi_K %*% a_k # Spectral projection
      return(preds)
    }
  }
}
attr(make_least_squares,"precompute_data") <- precompute_least_squares_data

