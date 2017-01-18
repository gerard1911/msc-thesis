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
em.maximization <- function(y,X,J,pi0,beta0,sigma0,w,equal.variances){
  
  n = length(y)
  k = ncol(X)
  # Weighted squared differences for equal variances
  totsum = 0
  
  # Fill in matrix
  for(j in 1:J){
    sum = 0     # sum of weights per component
    sum1 = 0    # sum of squared difference per component
    
    for(i in 1:n){
      sum = sum + w[i,j]
      if(pi0[j]!=0){
        sum1 = sum1 + w[i,j]*(y[i]-X[i,]%*%beta0[,j])^2
      }
    }
    
    # Update parameters
    pi0[j] = sum/n
    if(pi0[j]!=0){
      if(qr(t(X)%*% diag(w[,j]) %*% X)$rank == k ){
      beta0[,j] = solve(t(X)%*% diag(w[,j]) %*% X ) %*% t(X)%*% diag(w[,j])%*%as.matrix(y,nr=n)
    }else{
      lambda = 10^-7
      beta0[,j] = solve(t(X)%*% diag(w[,j]) %*% diag(w[,j]) %*% X + lambda * diag(k))%*% t(X)%*% diag(w[,j]) %*% diag(w[,j]) %*%as.matrix(y,nr=n) # NEED TO WORK ON THIS!!!   
    }
      if(equal.variances==FALSE){
        sigma0[j,j] = sum1/sum
      }
    }else{
      beta0[,j] = rep(NA,dim(X)[2])
      sigma0[j,j] = NA
    }
    totsum = totsum + sum1
  }
  
  if(equal.variances==TRUE){
    sigma0 = diag(as.numeric(totsum/n),nrow=J)
  }
  
  # Create output list
  maxim.output = list(pi = pi0, beta = beta0, sigma = sigma0)
  return(maxim.output)
}