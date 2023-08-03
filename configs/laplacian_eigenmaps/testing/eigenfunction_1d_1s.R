#------------------------------------------------#
# Configs for Laplacian Eigenmaps Figure 1.
#------------------------------------------------#

### General configs.
d <- 1
s <- 1
M <- 2^(s - 2*d + 1) 
ns <- round(seq(1000,4000,length.out = 10))
iters <- 200
B <- 100      # Number of permutations for permutation test.
alpha <- .05   # Accepted type I error.


### Configs for sampling. ###
sample_X <- make_sample_uniform(d)
make_f0s <- make_testing_eigenfunction_f0s


### Methods. ###
methods <- list(
  laplacian_eigenmaps = make_laplacian_eigenmaps_test,
  spectral_projection = make_spectral_projection_test
)
initialize_thetas <- list(
  initialize_laplacian_eigenmaps_test_thetas,
  initialize_spectral_projection_test_thetas
)
