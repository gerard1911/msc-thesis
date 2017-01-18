# ESTIMATED WEIGHTS EM-ALGORITHM WITH MEAN-SPECIFIC LINEAR COMBINATIONS

# INPUT
#
#   y in R(nx1),                response variable
#   X in R(nxk),                matrix with expoloratory variables, default = NULL
#   J in N(1),                  number of components of mixture distribution, default = 1
#   pi in [0,1]^J,              Normalized probability vector with estimates
#   beta in R(kxJ),             Matrix with Estimated regression coefficients
#   sigma2 in R(JxJ)            Diagonal matrix with estimated variances

# WEIGHT MATRIX FUNCTION

em.weights <- function(y,X,J,pi0,beta0,sigma0){
  n = length(y)
  w = matrix(,nrow=n,ncol = J)
  for(i in 1:n){   
    for(j in 1:J){
      if(pi0[j] != 0){
        w[i,j] = pi0[j]*dnorm(y[i],mean=(X[i,]%*%beta0[,j]),sd=sqrt(sigma0[j,j]))
      }else{
        w[i,j] = 0
      }
      cons = 0  # total sum of weights per component
      for(k in 1:J){
        if(pi0[k]!=0 ){
          cons = cons + pi0[k]*dnorm(y[i],mean=(X[i,]%*%beta0[,k]),sd=sqrt(sigma0[k,k])) 
        }  
      }
      # If there is only one observation for the group, sigma ->0. This implies the other component's probabilities
      # are zeroes. Hence, if w[i,j] == Inf, the others have to be zero and so w[i,j] = 1.
      if(w[i,j] != Inf){
      w[i,j] = w[i,j]/ cons
      }else{
        w[i,j] = 1
      }
    }
  }
  return(w)
}