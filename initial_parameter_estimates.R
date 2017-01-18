init.param.estimates <- function(y,X,J){
  
  # Determine data dimensions
  n = length(y)
  K = ncol(X)
  # Create first estimates of theta
  pi.hat = rep(1/J,J)
  
  # Determine set sizes for J different groups.
  n.hat = round(pi.hat*n)
  if(sum(n.hat)>n){
    n.hat[J]  = n.hat[J] - 1
  }else if(sum(n.hat)<n){
    n.hat[1] = n.hat[1] + 1
  }
  
  # Sample data sets
  integs = seq(1,n)
  list.int = list()
  Xlist = list()
  ylist = list()
  beta.hat = matrix(,nr=K,nc=J)
  sigma.hat = diag(J)
  for(j in 1:J){
    list.int[[j]] = sort(sample(integs, size = n.hat[j],replace=FALSE))
    Xlist[[j]] = X[list.int[[j]],]
    ylist[[j]] = y[list.int[[j]]]
    beta.hat[,j] = ginv(t(Xlist[[j]])%*%Xlist[[j]])%*%t(Xlist[[j]]) %*%ylist[[j]]
    if(K>1){
    sigma.hat[j,j] = (t(ylist[[j]]-Xlist[[j]]%*%beta.hat[,j])%*%(ylist[[j]]-Xlist[[j]]%*%beta.hat[,j]))/n.hat[j]
    }else{
      sigma.hat[j,j] = sum((ylist[[j]]-Xlist[[j]]*beta.hat[,j])^2)/n.hat[j]
    }
    integs = NULL
    for(k in 1 :n){
      if(k %in% list.int[[j]]){
      }else{
        integs = c(integs,k)
      }
  }
  }
  return(list(pi = pi.hat, beta = beta.hat, sigma = sigma.hat))
}
