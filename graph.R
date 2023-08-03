neighborhood_graph <- function(X, r, X_new = X,loop = FALSE, kernel = function(s){s <= 1})
{
  # unpack necessary parameters
  N <- nrow(X)
  
  k_max <- round(min(4*N^(.5),N - 1))
  A <- knn_to_neighborhood_graph(X,r,k_max,X_new,loop, kernel)
  return(A)
}

knn_to_neighborhood_graph <- function(X, r, k_max, X_new = X,loop, kernel)
{
  #--------------------#
  # Input: X (n x d matrix)
  #        r connectivity parameter
  #        k_max maximum degree in graph
  # Output: A (m x n sparse adjacency matrix)
  # Function: compute adjacency matrix A
  # with A = (A_ij)_{i,j in 1:n}, A_ij = 1 iff ||X[i,] - X[j,]||_2 < r
  #--------------------#
  
  # unpack necessary parameters
  n <- nrow(X)
  m <- nrow(X_new)
  d <- ncol(X)
  repeat{
    A <- Matrix(0, nrow = n, ncol = m)
    
    # compute distance matrix
    rneighborhood <- nn2(data = X, query = X_new, searchtype = 'radius', k = k_max, radius = r)
    rneighbors <- rneighborhood$nn.idx
    dist_to_neighbors <- rneighborhood$nn.dists
    if( (!loop) & all(X_new == X)){
      rneighbors <- rneighbors[,-1]
      dist_to_neighbors <- dist_to_neighbors[,-1]
    }
    
    # convert to easy format for adding to sparse matrix.
    dist_to_neighbors[rneighbors == 0] <- NA
    rneighbors[rneighbors == 0] <- NA
    
    r_neighbors_list <- as.matrix(melt(rneighbors)[,-2])
    r_neighbors_list <- r_neighbors_list[rowSums(is.na(r_neighbors_list)) == 0,]
    r_neighbors_list <- r_neighbors_list[,c(2,1)]
    
    dist_to_neighbors_vec <- c(dist_to_neighbors)
    dist_to_neighbors_vec <- dist_to_neighbors_vec[!is.na(dist_to_neighbors_vec)]
    W <- kernel(dist_to_neighbors_vec/r)
    
    # Add to matrix
    A_check <- A  # for checking convergence
    A_check[r_neighbors_list] <- 1
    A[r_neighbors_list] <- W
    
    # if k_max was not high enough, at least one vertex in A_check will have exactly k_max - 1 neighbors.
    # Increase k_max and recalculate.
    if( !(max(rowSums(A_check)) == k_max - 1) | !(k_max < n - 1)) break
    k_max <- min(2*k_max,n - 1)
  }
  
  return(A)
}

knn_graph <- function(X, k)
{
  # unpack necessary parameters
  n <- nrow(X)
  d <- ncol(X)
  A <- Matrix(0, nrow = n, ncol = n)
  
  # find neighbors
  neighbors <- nn2(data = X, query = X, k = k + 1)$nn.idx[,-1]
  neighbors_list <- as.matrix(melt(neighbors)[,-2])
  A[neighbors_list] <- 1
  
  # symmetrize
  A <- (A + t(A)) / 2
  
  return(A)
}

Laplacian <- function(A)
{
  # unpack necessary parameters
  n <- nrow(A)
  
  # Laplacian matrix
  D <- Matrix(0, ncol = n, nrow = n, sparse = T)
  diag(D) <- rowSums(A)
  if(class(D - A) == "dsyMatrix")
  {
    L <- as(D - A,"matrix")
  } else{
    L <- as(D - A,"dgCMatrix")
  }
  return(L)
}

incidence <- function(A)
{
  # unpack necessary parameters
  n <- nrow(A)
  m <- sum(A > 0)
  
  # Get adjacency representation in list form.
  adj_list <- which(G > 0, arr.ind = T)
  adj_list <- cbind(adj_list, G[adj_list], 1:m)
  
  # Convert to matrix form.
  D <- Matrix(0,nrow = m, ncol = n)
  D[adj_list[,c(4,1)]] <- adj_list[,3]
  D[adj_list[,c(4,2)]] <- -adj_list[,3]
  
  return(D)
}

# A wrapper around eigs_sym, which takes care of some failure cases when the
# matrix is poorly conditioned. 
# 
# Input: A, matrix.
#        K, number of eigenvectors.
#        retvec, do we want eigenvectors or just eigenvalues.
#        tol, tolerance for error in iterative computation.
#        sigma, for shift-and-invert
get_spectra <- function(A,K,retvec = TRUE,tol = 1e-10,sigma = tol){
  # Compute as many eigenvectors as we will need.
  maxitr <- 10000
  spectra <- tryCatch(eigs_sym(A,K,sigma = sigma,
                               which = "LM", 
                               opts = list(maxitr = maxitr,
                                           retvec = retvec,
                                           tol = tol)),
                      error = function(e){
                        test <- eigs_sym(A,K, which = "SM", opts = list(maxitr = maxitr,
                                                                        retvec = retvec,
                                                                        tol = tol))
                        return(test)
                      })
  
  # Sometimes, this will spit out garbage. 
  # You can diagnose this by the presence of ridiculously large eigenvalues, 
  # when we are looking for small eigenvalues.
  #
  # Check if there are any large eigenvalues, and if so, try again, decreasing
  # the tolerance for error; this seems to work.
  while(max(abs(spectra$values)) > 1e10)
  {
    tol <- tol * .001
    sigma <- tol
    # No need for try-catch, since we know the "LM" direction works now.
    spectra <- eigs_sym(A,K,sigma = sigma, 
                        which = "LM", 
                        opts = list(maxitr = maxitr,retvec = retvec,tol = tol))
    
                      
  }
  
  # Sometimes, this won't converge.
  # If it hasn't converged, up the number of iterations and try again.
  while(length(spectra$values) < K)
  {
    maxitr <- maxitr + 10000
    spectra <- eigs_sym(A,K,which = "SM",opts = list(maxitr = maxitr,
                                                     retvec = retvec,
                                                     tol = tol))
  }
  spectra
}
