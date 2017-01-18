# BIC score for mixtures of graphical models

# Input
#
# S.list, list with J Sample covariances matrices, 
# P.list, list with J sparse precision matrices from glasso
# Ncom, vector of size J with sample size per component


BIC.MGM <- function(S.list, P.list,Ncom){
	Jc = length(Ncom)
	p = dim(S.list)[1]
	
	# Calculate loglikelihood
	LL = 0
	for(j in 1:Jc){
		LL = LL + Ncom[j]/2*(log(det(P.list[,,j])) - sum(diag(S.list[,,j]%*%P.list[,,j]) ))
	}
	
	# Calculate no. of estimated parameters
	nmu = Jc*p
	npi = Jc-1
	nsig = NULL
	for(k in 1:Jc){
		for(i in 1:p){
			for(j in 1:p){
				if(j<i){
				S.list[i,j,k] = 0
				}
			}
		
		}
		nsig[k] = sum(S.list[,,k] != 0)
		
	}
	print(nmu + npi + sum(nsig))
	BIC = -2*LL + log(sum(Ncom)) * (nmu + npi + sum(nsig))
	return(BIC)
}