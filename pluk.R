# YEAST DATA ANALYSIS
library(igraph)
library(MASS)

yeast.average.data.set = list()
yeast.full.data.set = list()
gen.names = c('CBF1','GAL4','SWI5','GAL80','ASH1')


for(i in 1:5){
  yeast.average.data.set[[i]] = yeast.time.dependent[[9+i]]  
  yeast.full.data.set[[i]] = yeast.time.dependent[[i]]
}

# Create lists
conf.score = list()
gen.score = list()
auroc.score = list()

# Apply for EM-algorithm, different number of components, non-fixed weights, gene(t+1) ~gene (t).
for(j in 1:2){
conf.score[[j]] = config.scores(dataset = yeast.average.data.set, FUN = em.algorithm, nit =1 ,n.monte=1000,J=j, print.step = TRUE,) 
gen.score[[j]] = gene.scores(conf.score[[j]])
colnames(gen.score[[j]]) = gen.names
rownames(gen.score[[j]]) = gen.names
auroc.score[[j]] = auroc(gen.score[[j]], (yeast.true.relation))
print(auroc.score[[j]])
par(mfrow=c(1,2))
draw.gene.graph(relation.matrix=yeast.true.relation,gen.matrix=gen.score[[j]],k=j)
plot(auroc.score[[j]]$ROC[,1],auroc.score[[j]]$ROC[,2],main = 'ROC-curve', xlab = 'false positive',
     ylab = 'true positive')
lines(auroc.score[[j]]$ROC[,1],auroc.score[[j]]$ROC[,2])
}

# Apply for fixed weights
j=3
conf.score[[j]] = config.scores(dataset = yeast.average.data.set, FUN = em.algorithm,
                                true.weights=matrix(c(rep(1,14),rep(0,19),rep(0,14),rep(1,19)),nc = 2),
                                nit =1 ,n.monte=1000, J=2, print.step = TRUE) 
gen.score[[j]] = gene.scores(conf.score[[j]])
colnames(gen.score[[j]]) = gen.names
rownames(gen.score[[j]]) = gen.names
auroc.score[[j]] = auroc(gen.score[[j]], (yeast.true.relation))
print(auroc.score[[j]])
par(mfrow=c(1,2))
draw.gene.graph(relation.matrix=yeast.true.relation,gen.matrix=gen.score[[j]],k=j)
plot(auroc.score[[j]]$ROC[,1],auroc.score[[j]]$ROC[,2],main = 'ROC-curve', xlab = 'false positive',
     ylab = 'true positive')
lines(auroc.score[[j]]$ROC[,1],auroc.score[[j]]$ROC[,2])

# gene(t) ~gene(t)

for(j in 4:5){
  conf.score[[j]] = config.scores(dataset = yeast.time.independent, FUN = em.algorithm, nit =1 ,
                                  n.monte=1000, J=j, print.step = TRUE) 
  gen.score[[j]] = gene.scores(conf.score[[j]])
  colnames(gen.score[[j]]) = gen.names
  rownames(gen.score[[j]]) = gen.names
  auroc.score[[j]] = auroc(gen.score[[j]], (yeast.true.relation))
  print(auroc.score[[j]])
  par(mfrow=c(1,2))
  draw.gene.graph(relation.matrix=yeast.true.relation,gen.matrix=gen.score[[j]],k=j)
  plot(auroc.score[[j]]$ROC[,1],auroc.score[[j]]$ROC[,2],main = 'ROC-curve', xlab = 'false positive',
       ylab = 'true positive')
  lines(auroc.score[[j]]$ROC[,1],auroc.score[[j]]$ROC[,2])
}

# Apply for fixed weights
j=6
conf.score[[j]] = config.scores(dataset = yeast.time.independent, FUN = em.algorithm,
                                true.weights= matrix(c(rep(1,15),rep(0,20),rep(0,15),rep(1,20)),nc=2),
                                                    nit =1 ,n.monte=1000, J=2, print.step = TRUE) 
gen.score[[j]] = gene.scores(conf.score[[j]])
colnames(gen.score[[j]]) = gen.names
rownames(gen.score[[j]]) = gen.names
auroc.score[[j]] = auroc(gen.score[[j]], (yeast.true.relation))
print(auroc.score[[j]])
par(mfrow=c(1,2))
draw.gene.graph(relation.matrix=yeast.true.relation,gen.matrix=gen.score[[j]],k=j)
plot(auroc.score[[j]]$ROC[,1],auroc.score[[j]]$ROC[,2],main = 'ROC-curve', xlab = 'false positive',
     ylab = 'true positive')
lines(auroc.score[[j]]$ROC[,1],auroc.score[[j]]$ROC[,2])


#_________________________________________________________________________________________________________________


# Apply for CEM-algorithm, different number of components, non-fixed weights, gene(t+1) ~gene (t).
for(j in 7:8){
  conf.score[[j]] = config.scores(dataset = yeast.average.data.set, FUN = cem.algorithm, nit =1 ,n.monte=1000,
                                  J=j, print.step = TRUE) 
  gen.score[[j]] = gene.scores(conf.score[[j]])
  colnames(gen.score[[j]]) = gen.names
  rownames(gen.score[[j]]) = gen.names
  auroc.score[[j]] = auroc(gen.score[[j]], (yeast.true.relation))
  print(auroc.score[[j]])
  par(mfrow=c(1,2))
  draw.gene.graph(relation.matrix=yeast.true.relation,gen.matrix=gen.score[[j]],k=j)
  plot(auroc.score[[j]]$ROC[,1],auroc.score[[j]]$ROC[,2],main = 'ROC-curve', xlab = 'false positive',
       ylab = 'true positive')
  lines(auroc.score[[j]]$ROC[,1],auroc.score[[j]]$ROC[,2])
}

# Apply for fixed weights
j=9
conf.score[[j]] = config.scores(dataset = yeast.average.data.set, FUN = cem.algorithm,
                                true.weights=matrix(c(rep(1,14),rep(0,19),rep(0,14),rep(1,19)),nc = 2),
                                nit =1 ,n.monte=1000, J=2, print.step = TRUE) 
gen.score[[j]] = gene.scores(conf.score[[j]])
colnames(gen.score[[j]]) = gen.names
rownames(gen.score[[j]]) = gen.names
auroc.score[[j]] = auroc(gen.score[[j]], (yeast.true.relation))
print(auroc.score[[j]])
par(mfrow=c(1,2))
draw.gene.graph(relation.matrix=yeast.true.relation,gen.matrix=gen.score[[j]],k=j)
plot(auroc.score[[j]]$ROC[,1],auroc.score[[j]]$ROC[,2],main = 'ROC-curve', xlab = 'false positive',
     ylab = 'true positive')
lines(auroc.score[[j]]$ROC[,1],auroc.score[[j]]$ROC[,2])

# gene(t) ~gene(t)

for(j in 10:11){
  conf.score[[j]] = config.scores(dataset = yeast.time.independent, FUN = cem.algorithm, nit =1 ,
                                  n.monte=1000, J=j, print.step = TRUE) 
  gen.score[[j]] = gene.scores(conf.score[[j]])
  colnames(gen.score[[j]]) = gen.names
  rownames(gen.score[[j]]) = gen.names
  auroc.score[[j]] = auroc(gen.score[[j]], (yeast.true.relation))
  print(auroc.score[[j]])
  par(mfrow=c(1,2))
  draw.gene.graph(relation.matrix=yeast.true.relation,gen.matrix=gen.score[[j]],k=j)
  plot(auroc.score[[j]]$ROC[,1],auroc.score[[j]]$ROC[,2],main = 'ROC-curve', xlab = 'false positive',
       ylab = 'true positive')
  lines(auroc.score[[j]]$ROC[,1],auroc.score[[j]]$ROC[,2])
}

# Apply for fixed weights
j=12
conf.score[[j]] = config.scores(dataset = yeast.time.independent, FUN = cem.algorithm,
                                true.weights= matrix(c(rep(1,15),rep(0,20),rep(0,15),rep(1,20)),nc=2),
                                nit =1 ,n.monte=1000, J=2, print.step = TRUE) 
gen.score[[j]] = gene.scores(conf.score[[j]])
colnames(gen.score[[j]]) = gen.names
rownames(gen.score[[j]]) = gen.names
auroc.score[[j]] = auroc(gen.score[[j]], (yeast.true.relation))
print(auroc.score[[j]])
par(mfrow=c(1,2))
draw.gene.graph(relation.matrix=yeast.true.relation,gen.matrix=gen.score[[j]],k=j)
plot(auroc.score[[j]]$ROC[,1],auroc.score[[j]]$ROC[,2],main = 'ROC-curve', xlab = 'false positive',
     ylab = 'true positive')
lines(auroc.score[[j]]$ROC[,1],auroc.score[[j]]$ROC[,2])




#___________________________
# Apply for graphical model
j=13
S = var(yeast.data.list$average.data[,1:5])*34/35
N = 35

# Select max rho
x = glassopath(S)
i=1
while(sum(x$wi[,,i] != 0) != 5){
	i = i + 1
}
rhomax = x$rholist[i]
# Create list of sparse matrix estimates
x = glassopath(S,rholist = seq(from =0,by = rhomax/1000,to = rhomax))$wi
# Create BIC values for every estimate
BIClist = NULL
Slist = array(,dim = c(5,5,1))
Plist = Slist
Slist[,,1] = S
BICmin = Inf
idmin = 0
for(i in 1:1001){
	Plist[,,1] = x[,,i] 
	BIClist[i] = BIC.MGM(S.list = Slist, P.list = Plist,Ncom = N)
	if(BICmin > BIClist[i]){
		BICmin = BIClist[i]
		idmin = i
	}
}
diag(x[,,idmin]) = NA
gen.score[[j]] = (x[,,idmin]!=0)
colnames(gen.score[[j]]) = gen.names
rownames(gen.score[[j]]) = gen.names
auroc.score[[j]] = auroc(as.matrix(gen.score[[j]]), (yeast.true.relation.und))
print(auroc.score[[j]])
par(mfrow=c(1,2))
draw.gene.graph(relation.matrix=yeast.true.relation.und,gen.matrix=gen.score[[j]],m='undirected')
plot(auroc.score[[j]]$ROC[,1],auroc.score[[j]]$ROC[,2],main = 'ROC-curve', xlab = 'false positive',
     ylab = 'true positive')
lines(auroc.score[[j]]$ROC[,1],auroc.score[[j]]$ROC[,2])




#______________
# MGM
j=14
S1 = var(yeast.data.list$average.data[1:15,1:5])*14/15
S2 = var(yeast.data.list$average.data[1:15,1:5])*19/20
N1 = c(15,20)

# Select max rho
x = glassopath(S1)
i=1
while(sum(x$wi[,,i] != 0) != 5){
	i = i + 1
}
rhomax1 = x$rholist[i]
# Select max rho
x = glassopath(S2)
i=1
while(sum(x$wi[,,i] != 0) != 5){
	i = i + 1
}
rhomax2 = x$rholist[i]

# Create list of sparse matrix estimates
x1 = glassopath(S1,rholist = seq(from =0,by = rhomax1/1000,to = rhomax1))$wi
x2 = glassopath(S2,rholist = seq(from =0,by = rhomax2/1000,to = rhomax2))$wi

# Create BIC values for every estimate
BIClist = matrix(,nrow=1001,ncol=1001)
Slist = array(,dim = c(5,5,2))
Slist[,,1] = S1
Slist[,,2] = S2
Plist = array(,dim = c(5,5,2))
Slist[,,1] = S
BICmin = Inf
idmin = c(0,0)
for(i in 1:1001){
	for(j in 1:1001){
	Plist[,,1] = x1[,,i] 
	Plist[,,2] = x2[,,j]
	BIClist[i,j] = BIC.MGM(S.list = Slist, P.list = Plist,Ncom = N)
	if(BICmin > BIClist[i,j]){
		BICmin = BIClist[i,j]
		idmin=c(i,j)
	}
}
}
j=14
diag(x1[,,idmin[1]]) = 0
diag(x2[,,idmin[2]]) = 0
x = (x1[,,idmin[1]] != 0) + (x2[,,idmin[2]] != 0)
diag(x) = NA
gen.score[[j]] = (x!=0)
colnames(gen.score[[j]]) = gen.names
rownames(gen.score[[j]]) = gen.names
auroc.score[[j]] = auroc(gen.score[[j]], (yeast.true.relation.und))
print(auroc.score[[j]])
par(mfrow=c(1,2))
draw.gene.graph(relation.matrix=yeast.true.relation.und,gen.matrix=gen.score[[j]],m='undirected')
plot(auroc.score[[j]]$ROC[,1],auroc.score[[j]]$ROC[,2],main = 'ROC-curve', xlab = 'false positive',
     ylab = 'true positive')
lines(auroc.score[[j]]$ROC[,1],auroc.score[[j]]$ROC[,2])

