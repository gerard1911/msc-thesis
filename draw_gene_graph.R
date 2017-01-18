draw.gene.graph <- function(relation.matrix,gen.matrix,k){
  
  # Run through matrix
  edge.type = NULL
  edge.curve = NULL
  edge.color = NULL
  for(j in 1:5){
    for(i in 1:5){
      if(i!=j){
      
        if(gen.matrix[i,j] == 1){ 
          if(relation.matrix[i,j] ==1){  # TRUE POSITIVES, DRAWN IN BLACK
          
            if(gen.matrix[i,j] == gen.matrix[j,i]){  # If i<->j, use curved edges.
              edge.curve = c(edge.curve,1)
            }else{
              edge.curve= c(edge.curve,0)
            }
          
            if((i ==4 && j ==2) || (i==2 && j ==4)){ # For indirect relation, use striped edges.
              edge.type = c(edge.type,2)   
            }else{
              edge.type = c(edge.type,1)
            }
            edge.color = c(edge.color,'black') # Use black as color.
          
          }else{ #FALSE POSITIVES, DRAWN IN GREY
            if(gen.matrix[i,j] == gen.matrix[j,i]){ # If i<->j, use curved edges.
              edge.curve = c(edge.curve,1) 
            }else{
              edge.curve= c(edge.curve,0)
            }
          
            if((i ==4 && j ==2) || (i==2 && j ==4)){ # For indirect relation, use striped edges.
              edge.type = c(edge.type,2)   
            }else{
              edge.type = c(edge.type,1)
            }
            edge.color = c(edge.color,'grey')
          } # end false positives
        }      # end genmatrix == 1
      } 
    }
  }
  # Create graph data.
  g = graph.adjacency(t(gen.matrix))
  
  # Use the plot function.
  plot(g, vertex.size = 45, main = paste('Estimated network, EM with J=',k),layout = matrix(c(1,2,3,4,5,1,0,2,0,1),nc = 2),
     edge.lty =edge.type, edge.curved = edge.curve, edge.width = 2, edge.color = edge.color)
}