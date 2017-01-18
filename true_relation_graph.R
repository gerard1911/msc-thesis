# Graph of relations
library(igraph)
gen.names = c('CBF1','GAL4','SWI5','GAL80','ASH1')
yeast.true.relation = matrix(c(0,1,0,0,0,0,0,1,1,0,1,0,0,1,1,0,1,0,0,0,1,0,0,0,0), nr =5)
yeast.true.relation.und = (yeast.true.relation + t(yeast.true.relation))!=0
colnames(yeast.true.relation) = gen.names
rownames(yeast.true.relation) = gen.names

g = graph.adjacency(t(yeast.true.relation))
plot(g, vertex.size = 45, main = 'True Network',layout = matrix(c(1,2,3,4,5,1,0,2,0,1),nc = 2), edge.lty =c(1,1,2,1,1,1,2,1,1), edge.curved = c(0,0,1,0,0,0,1,0,0), edge.width = 2, edge.color = 'black'
)