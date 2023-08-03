#------------------------------------------------#
# Configs for Laplacian Eigenmaps Figure 1.
#------------------------------------------------#

### General configs.
d <- 2
s <- 1
M <- 2^(s - 2*d + 1)
ns <- round(seq(1000,4000,length.out = 10))
iters <- 200

### Configs for sampling. ###
sample_X <- make_sample_uniform(d)
make_f0 <- function(d,n){
  g0 <- make_eigenfunction(d,round((M^2*n)^(d/(2*s + d))),M,s)
  # g0 <- make_eigenfunction(d,16,M,s)
  f0 <- function(x){2 * g0(x)}
}


### Methods. ###
methods <- list(
  laplacian_eigenmaps = make_laplacian_eigenmaps,
  spectral_projection = make_spectral_projection
)
initialize_thetas <- list(
  initialize_laplacian_eigenmaps_thetas,
  initialize_spectral_projection_thetas
)