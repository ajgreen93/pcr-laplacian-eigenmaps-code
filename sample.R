#------------------------------------#
# Design distributions.
#------------------------------------#

make_sample_uniform <- function(d)
{
  sample_X <- function(n){
    X <- matrix(runif(n*d,-1,1),ncol = d)
    return(X)
  }
}

#-----------------------------------------------#
# Regression functions
#-----------------------------------------------#

make_eigenfunction <- function(d,K,M,s){
  #----------------------------------#
  # Builds the regression function
  #
  # f0(x) = M/lambda_K^(s/2) * psi_K
  #
  # where psi_K(x) is the kth element in the trigonometric basis.
  #----------------------------------#
  lambda_K <- get_eigenvalue(d,K)
  eigenfunction_K <- get_eigenfunction(d,K)
  f0 <- function(x){
    M/lambda_K^(s/2) * eigenfunction_K(x)
  }
  attr(f0,"type") <- "eigenfunction"
  return(f0)
}

make_testing_eigenfunction_f0s <- function(d,n){
  # Configuration type parameter
  max_K <- round(2*(M^2*n)^(2*d/(4*s + d)))
  
  # Get eigenfunctions, and associated l2 norms
  g0s <- lapply(1:max_K,FUN = function(k){make_eigenfunction(d,k,M,s)})
  l2_norms <- sapply(g0s,FUN = get_l2_norm)
  
  # Remove anything with infinite l2 norm (constant functions), or
  # duplicates in l2 norm.
  g0s <- g0s[l2_norms < Inf]
  l2_norms <- l2_norms[l2_norms < Inf]
  g0s <- g0s[match(unique(l2_norms),l2_norms)]
  l2_norms <- l2_norms[match(unique(l2_norms),l2_norms)]
  
  # Construct regression functions.
  make_f0 <- function(g,l){
    f0 <- function(x){2 * g(x)}
    attr(f0,"l2_norm") <- 4 * l
    return(f0)
  }
  f0s <- mapply(FUN = make_f0, g0s,l2_norms)
  f0_null <- function(x){0}; 
  attr(f0_null,"l2_norm") <- 0
  f0s <- c(f0_null,f0s)
  return(f0s)
}

make_spectral_sobolev <- function(d,n,s = 1,M = 2^s,beta = NULL){
  #------------------------------#
  # Fixes a vector beta that belongs to a Sobolev ellipsoid,
  # then builds the regression function
  # 
  # f0(x) = sum_{k} beta_k psi_k
  #
  # where psi_k(x) is the kth element in the trigonometric basis.
  #------------------------------#
  if(is.null(beta)){
    n_coefs <- n/20
    beta <- M/sqrt(n_coefs) * c(1/sqrt( (1:n_coefs)^(2*s/d) ) ) # Sobolev norm = 2^s.
  }
  trig_basis <- get_trigonometric_basis(d,length(beta))
  f0 <- function(x){
    t(beta) %*% trig_basis(x)
  }
}