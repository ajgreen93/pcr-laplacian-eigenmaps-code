#------------------------------------------------#
# Configs for Laplacian Eigenmaps Figure 3.
#------------------------------------------------#

### General configs.
d <- 2
s <- 1
M <- 4 * 2^(s - 2*d + 1)
ns <- 1000
iters <- 100

### Configs for sampling. ###
sample_X <- make_sample_uniform(d)
make_f0 <- function(d,n){
  g0 <- make_spectral_sobolev(d,n,s,M)
  f0 <- function(x){g0(x)}
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