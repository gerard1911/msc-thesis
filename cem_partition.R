cem.partition <- function(w){
  n = nrow(w)
  J = ncol(w)
  P = matrix(rep(0,n*J),nc = J, nr = n)
  
  # Determine partition matrix.
  for(i in 1:n){
    max.count = 1
    for(j in 1:J){
      if(w[i,j]>w[i,max.count]){
        max.count = j
      }
    }
    P[i,max.count] = 1
  }
  return(P)
}