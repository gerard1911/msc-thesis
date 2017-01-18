# ESTIMATED PARAMETERS EM-ALGORITHM WITH MEAN-SPECIFIC LINEAR COMBINATIONS

# INPUT
#
#   y in R(nx1),                response variable
#   X in R(nxk),                matrix with expoloratory variables, default = NULL
#   J in N(1),                  number of components of mixture distribution, default = 1
#   pi in [0,1]^J,              Normalized probability vector with estimates
#   beta in R(kxJ),             Matrix with Estimated regression coefficients
#   sigma2 in R(JxJ),           Diagonal matrix with estimated variances
#   w in R(nxJ),                Matrix with estimated weights per component
#   equal.var in {FALSE,TRUE},  logical variable declaring equal variances for every component, default = FALSE

# MAXIMIZATION STEP EM-ALGORITHM
cem.maximization <- function(y.list, X.list, K, Partition.matrix, equal.variances){
  
  n = nrow(Partition.matrix)
  k = K
  J = ncol(Partition.matrix)
  pi0 = rep(NA,J)
  beta0 = matrix(NA,nr=k, nc = J)
  sigma0 = matrix(0, nrow = J, ncol = J)
  # Estimated pi
  for(j in 1:J){
    pi0[j] = sum(Partition.matrix[,j]) / n
  }
  print(pi0)
  # Estimated beta
  for(j in 1:J){
    if(is.null(X.list[[j]]) == FALSE && is.na(X.list[[j]] )==FALSE){
    beta0[,j] = ginv(t(as.matrix(X.list[[j]], nr = nrow(X.list[[j]]))) %*% as.matrix(X.list[[j]],
                                                                                     nr = nrow(X.list[[j]]))) %*% 
  t(as.matrix(X.list[[j]], nr = nrow(X.list[[j]]))) %*% as.matrix(y.list[[j]], nc = 1)
    }else{
      beta[0,j] = 0
    }
  }
  
  # Estimated sigma
  for(j in 1:J){
    sigma0[j,j] = (t(as.matrix(y.list[[j]],nc=1) - X.list[[j]] %*% beta0[,j])%*%
      (as.matrix(y.list[[j]],nc=1) - X.list[[j]] %*% beta0[,j])) / sum(Partition.matrix[,j]) 
  }
  
  # Create output list
  maxim.output = list(pi = pi0, beta = beta0, sigma = sigma0)
  return(maxim.output)
}