stochastic.em.step <- function(weights){
  n = nrow(weights)
  J = ncol(weights)
  P = matrix(0,nr=n,nc =J)
  
  for(i in 1:n){
    P[i,] = rmultinom(n=1,size=1,prob=weights[i,]) 
  }
  return(P)
}