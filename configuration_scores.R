# CONFIGURATION SCORES OF ALGORITHM WITH MEAN-SPECIFIC LINEAR COMBINATIONS

# INPUT
#
#   dataset                     List with dataframe of interest
#   FUN                         Function to be used (EM, CEM, SEM)
#   nit in N(1),                number of requested iterations, default = 1
#   J in N(1),                  number of components of mixture distribution, default = 1
#   n.monte,                    number of required Monte Carlo iterations, default = 1
#   equal.var in {FALSE,TRUE},  logical variable declaring equal variances for every component, default = FALSE
#   print.step,                 logical variable declaring to print the currently running configuration.
#   true.weights in R(nxJ),     matrix with true weights, default = NULL.               


# Create function
config.scores <- function(dataset, true.weights = NULL, FUN, nit=1,J=1,n.monte=1,equal.var=FALSE, print.step = FALSE){
  
  # Create Configuration Score Matrix
  score_em = matrix(rep(0,75),nrow=5,ncol=15)
  bicscores = matrix(rep(0,75),nrow=5,ncol=15)
  loglscores = matrix(rep(0,75),nrow=5,ncol=15)
  parscores = matrix(rep(0,75),nrow=5,ncol=15)
  
  
  
  for(h in 1:nit){
    print(h)
    
    # EM-algorithm on all possible configurations, create empty datalist per configuration.
    yeast.analysis  = list()
    k = 1
    t=1
    
    # Iterate for all 5 response genes
    while(k <=5){
      data = dataset[[k]]
      
      # Single dependencies
      for(j in 2:5){
        names = c(colnames(data)[1],colnames(data)[j])
        emp = FUN(y=data[,1],X=data[,j],true.weights=true.weights,J=J, n.monte = n.monte,equal.var=equal.var )
        emp$covariates = paste(names[1],'<-', names[2])
        yeast.analysis[[t]] = emp
        if(print.step == TRUE){
          print(paste('configuration at',round(100*t/75)))
          print(emp$covariates)
        }
        t = t+1
      }
      
      # Double dependencies
      
      for(j in 2:5){
        for(l in 2:5){
          if( l>j){
            names = c(colnames(data)[1],colnames(data)[j],colnames(data)[l])      
            emp = FUN(y=data[,1],X=cbind(data[,j],data[,l]),true.weights=true.weights,J=J, n.monte = n.monte,equal.var=equal.var) 
            emp$covariates = paste(names[1],'<-', names[2],',',names[3])          
            yeast.analysis[[t]] = emp
            if(print.step == TRUE){
              print(paste('configuration at',round(100*t/75)))
              print(emp$covariates)
            }
            t = t+1          
          }
          
        }
        
      }
      
      # Triple dependencies
      
      for(j in 2:5){
        for(l in 2:5){
          for(m in 2:5){          
            if(l>j && m>l){  
              names = c(colnames(data)[1],colnames(data)[j],colnames(data)[l],colnames(data)[m])
              emp = FUN(y=data[,1],X=cbind(data[,j],data[,l],data[,m]),true.weights=true.weights,J=J, n.monte = n.monte,equal.var=equal.var)            
              emp$covariates = paste(names[1],'<-', names[2],',',names[3],',',names[4])            
              yeast.analysis[[t]] = emp
              if(print.step == TRUE){
                print(paste('configuration at',round(100*t/75)))
                print(emp$covariates)
              }
              t = t+1
            }
          }
        }
      }
      
      #Quadruple dependencies
      
      for(j in 2:5){  
        for(l in 2:5){
          for(m in 2:5){
            for(n in 2:5){            
              if(l>j && m>l && n>m){              
                names = c(colnames(data)[1],colnames(data)[j],colnames(data)[l],colnames(data)[m],colnames(data)[n])  
                emp = FUN(y=data[,1],X=cbind(data[,j],data[,l],data[,m],data[,n]),true.weights=true.weights,J=J, n.monte = n.monte,equal.var=equal.var)
                emp$covariates = paste(names[1],'<-', names[2],',',names[3],',',names[4],',',names[5])          
                yeast.analysis[[t]] = emp
                if(print.step == TRUE){
                  print(paste('configuration at',round(100*t/75)))
                  print(emp$covariates)
                }
                t = t+1              
              }            
            }          
          }        
        }      
      }
      k = k+1 
    }    
    
    # Determine best BIC score
    p = 0
    for(j in 1:5){
      k=1
      i=1
      BIC = Inf
      while(i<=15){
        if(yeast.analysis[[i+p]]$Bic < BIC){
          BIC = yeast.analysis[[i+p]]$Bic
          k = i
        }
        i = i+1
      }
      score_em[j,k]=score_em[j,k]+1 
      p = p + 15
    }
  }
  return(conf.score = score_em)
}