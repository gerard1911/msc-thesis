# msc-thesis
R-code of MSc-thesis: 'Comparing graphical models with pseudo-likelihood methods and their Bayesian counterparts on gene regulatory networks'. at University of Groningen.

Be sure to install the MASS, MVTNORM, igraph and glasso package.

For the simulation study, open artificial_data_analysis.r. 
  - For creating mixture linear data, use create.mixture.linear.data(N,J) (by default one component)
  - For creating mixture multivariate data, use create.mvn.data(N,J) (by default one component)
  
For the complete yeast data analysis, open pluk.r
