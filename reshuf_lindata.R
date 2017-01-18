# Reshuffle data

res.data <-function(datalist){
	n =nrow(datalist[[1]])+1
	m = ncol(datalist[[1]])
	dat = matrix(,nrow=n,ncol=m)
	for(j in 1:5){
		dat[n,j] = datalist[[j]][n-1,1]
	}
	dat[1:(n-1),1] = datalist[[2]][,2]
	dat[1:(n-1),2] = datalist[[1]][,2]
	dat[1:(n-1),3] = datalist[[1]][,3]
	dat[1:(n-1),4] = datalist[[1]][,4]
	dat[1:(n-1),5] = datalist[[1]][,5]
	return(dat)
}