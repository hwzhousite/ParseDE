library(PARSE)
library(Seurat)

constructNull_Guassian <- function(mat, cov_mat = NULL) {


  {
    mat <- as.matrix(mat)
    n_gene <- dim(mat)[1]
    n_cell <- dim(mat)[2]
    gene_names <- rownames(mat)
    para_feature <- rownames(mat)

    ## Marginal fitting

    para <- cbind(apply(mat,1,mean), apply(mat,1,sd))
    para <- simplify2array(para)
    names(para) <- para_feature
    if(sum(is.na(para)) > 0) {
      warning("NA produces in mean estimate; using 0 instead.")
      para[is.na(para)] <- 0
    }



    ## Covariance fitting
    if(is.null(cov_mat)){

      cov_mat <- cov(t(mat))

      for (i in 1:n_gene) {

        cov_mat[i,i] = cov_mat[i,i] + 1e-5

      }

    }


    # Null Data Generation
    if(n_gene > 1000){

      new_mvn <- mvnfast::rmvn(n_cell, mu = para[,1], sigma = cov_mat)

    }else{

      new_mvn <- mvtnorm::rmvnorm(n_cell, mean = para[,1], sigma = cov_mat)

    }

  }


  return(new_mvn)
}



# Inputs:
#   mat    : G x C numeric matrix (genes x cells)
#   mu_est : G x K numeric matrix (class means per gene)
#   label  : length-C integer vector with values in 1:K
#
# Output:
#   cov_mat: G x G pooled covariance matrix
constructNull_Guassian_pool <- function(mat, label, mu_est, cov_mat = NULL) {


  # Ensure basic types/shapes (cheap guards)
  mat    <- as.matrix(mat)
  mu_est <- as.matrix(mu_est)
  storage.mode(mat)    <- "double"
  storage.mode(mu_est) <- "double"

  G <- nrow(mat)
  C <- ncol(mat)
  K <-nrow(mu_est)
  # Cov Estimated with Heter Mean
  # Pooled Covariance Matrix
  # Build per-cell mean matrix by column indexing (G x C), no one-hot needed
  # Each column i takes the class-mean vector mu_est[, label[i]]
  mu_mat <- t(mu_est)[, label, drop = FALSE]

  # Centered data
  Xc <- mat - mu_mat

  # Pooled covariance: (Xc %*% t(Xc)) / (C - 1), using tcrossprod for speed
  cov_mat <- tcrossprod(Xc) / (C - 1)
  cov_mat <- as.matrix(cov_mat)
  # Diagonal jitter for numerical stability
  for (i in 1:G) {

    cov_mat[i,i] <- cov_mat[i,i] + 1e-5

  }


  # Use matrixStats if available (faster), else fallback
  mu <- apply(mat,1,mean)


  if(G > 1000){

    new_mvn <- mvnfast::rmvn(C, mu = mu, sigma = cov_mat)

  }else{

    new_mvn <- mvtnorm::rmvnorm(C, mean = mu, sigma = cov_mat)

  }


  return(new_mvn)
}


score_func <- function(mu, mu_mean = NULL,sigma, gene){
  K = nrow(mu)
  comb_mat = matrix(0, nrow = 1, ncol =ncol(mu))

  for (i in 1:K) {

    for (j in 1:K) {

      comb_mat[ ,gene] = comb_mat[ ,gene] + (mu[i,gene] - mu[j,gene])^2/sigma

    }

  }

  return(comb_mat/2)

}
