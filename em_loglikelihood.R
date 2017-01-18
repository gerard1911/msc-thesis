# LOGLIKELIHOOD EM-ALGORITHM WITH MEAN-SPECIFIC LINEAR COMBINATIONS

# INPUT
#
#   y in R(nx1),                response variable
#   X in R(nxk),                matrix with expoloratory variables, default = NULL
#   J in N(1),                  number of components of mixture distribution, default = 1
#   pi in [0,1]^J,              Normalized probability vector with estimates
#   beta in R(kxJ),             Matrix with Estimated regression coefficients
#   sigma2 in R(JxJ)            Diagonal matrix with estimated variances

# Create function
em.loglikelihood <- function(y,X,J,pi,beta,sigma2){
	n = length(y)
	q = 0
	for(i in 1:n){
		t=0
		for(j in 1:J){
      if(pi[j] != 0){
			t = t + pi[j]*dnorm(y[i],mean = X[i,]%*%beta[,j], sd = sqrt(sigma2[j,j]))
      }
	}
			q = q + log(t)
	}
	return(q)
}
