# Reshuffle mixture data

res.mix.data <-function(datalist){
	n =nrow(datalist[[1]])+2
	m = ncol(datalist[[1]])
	n1 = round(n/2)
	dat = matrix(,nrow=n,ncol=m)
	for(j in 1:5){
		dat[25,j] = datalist[[j]][24,1]
		dat[50,j] = datalist[[j]][48,1]
	}
	dat[1:24,1] = datalist[[2]][1:24,2]
	dat[1:24,2] = datalist[[1]][1:24,2]
	dat[1:24,3] = datalist[[1]][1:24,3]
	dat[1:24,4] = datalist[[1]][1:24,4]
	dat[1:24,5] = datalist[[1]][1:24,5]
	dat[26:49,1] = datalist[[2]][25:48,2]
	dat[26:49,2] = datalist[[1]][25:48,2]
	dat[26:49,3] = datalist[[1]][25:48,3]
	dat[26:49,4] = datalist[[1]][25:48,4]
	dat[26:49,5] = datalist[[1]][25:48,5]
	return(dat)
}