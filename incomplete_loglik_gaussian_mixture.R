# The incomplete log-likelihood function
ell <- function(y,X,J,pi,beta,sigma2){
	n = length(y)
	q = 0
	for(i in 1:n){
		t=0
		for(j in 1:J){
      if(pi[j]!=0){
			t = t + pi[j]*dnorm(y[i],mean = X[i,]%*%beta[,j], sd = sqrt(sigma2[j,j]))
      }else{
        t = t
      }
	}
			q = q + log(t)
	}
	return(q)
}
