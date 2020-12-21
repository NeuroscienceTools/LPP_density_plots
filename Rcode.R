# /******************************************************************************
#   ****************************************************************************** 
#   ** 
#   ** Copyright (c) 2020 VRVis Zentrum für Virtual Reality und Visualisierung
# ** Forschungs-GmbH All rights reserved.
# ** 
#   ************************************************************ 
#   ** 
#   ** THIS IS PUBLISHED PROPRIETARY SOURCE CODE OF VRVis GmbH The copyright
# ** notice above does not evidence any actual or intended publication of such
# ** code.
# ** 
#   ******************************************************************************
#   ******************************************************************************/
#   
#   // author: Florian Ganglberger <ganglberger@vrvis.at>


library(geigen)
library(MASS)

########################## GENERATE RANDOM DATA ###########################
############# YOU CAN IGNORE THIS IF YOU USE YOUR OWN DATA ################

#Generate random connectivity (correlation matrices) that change over time
#Note that this randomly genenerated data only simulates intra cluster variability (no inter cluster variablity!)

amount_of_nodes<-600   #Amount of nodes the connectivity/correlation matrix should have
amount_of_timepoints<-8   #Amount of timepoints 
amount_of_clusters_to_observe<-5  #The amount of clusters you want to visualize via the density map (colors in the plot)
ratio_of_connectivity_change<-0.2 #ratio of nodes that will change in the connectivity matrix

#generate an artificial correlation matrix as basis for the connectivity matrix
artificial_connectivity_matrix<-cor(matrix(cumsum(rnorm(amount_of_nodes*amount_of_nodes)),amount_of_nodes,amount_of_nodes))

#compute clusters that will have intra cluster variability over time
clusters<-kmeans(artificial_connectivity_matrix,floor(amount_of_clusters_to_observe/ratio_of_connectivity_change))$cluster

#connectivity_timeseries is a list that contains the connectiviy matrix over time
connectivity_timeseries<-list()
for(timepoint in 1:amount_of_timepoints){
  connectivity_timeseries[[timepoint]]<-artificial_connectivity_matrix
}

for(cluster in 1:amount_of_clusters_to_observe){
  random_connectivity_change<-cumsum(rnorm(amount_of_timepoints))/5 #compute a random change over all timepoints
  
  for(timepoint in 1:amount_of_timepoints){
    #adapt the connectivity for clusters
    connectivity_timeseries[[timepoint]][clusters==cluster,clusters==cluster]<-artificial_connectivity_matrix[clusters==cluster,clusters==cluster]+random_connectivity_change[timepoint]
  }
}

#compute mean connectiviy matrix
mean_connectivity_matrix <- matrix(0,amount_of_nodes,amount_of_nodes)

for(timepoint in 1:amount_of_timepoints){
  mean_connectivity_matrix<-mean_connectivity_matrix+connectivity_timeseries[[timepoint]]
}
mean_connectivity_matrix<-mean_connectivity_matrix/amount_of_timepoints







########## Visualising the Transition of Large Networks via Dimensionality Reduction ############


#Generate an affinity matrix out of the matrix X and fold it with a kernel function with sigma=t
getAffinityMatrix<-function(X,t=3){
  #Leave only the top 10% connections per row
  isTopRow<-t(apply(X,1,function(x){
    x>quantile(x,0.9)
  }))
  X[!(isTopRow)]<-0
  X[X<0]<-0
  
  #COSIM FUNCTION to make top10% matrix symetric again
  mat=tcrossprod(X, X)
  t1=sqrt(apply(X, 1, crossprod))
  X <- mat / outer(t1,t1)
  
  # distances need to be transformed according to 
  # original LPP paper by a kernel function. 
  # He, Xiaofei & Niyogi, Partha. (2002). Locality Preserving Projections (LPP).
  # IEEE Transactions on Reliability - TR. 16. 
  # The following therm does
  # this similarly: It creates a distance out of a similarity (1-X) and transforms it
  # like in the kernel function where t=3
  X<-exp(-((1-X)^2)/t)
  
  return(X)
}

#Perform LPP
X<-getAffinityMatrix(mean_connectivity_matrix)
D <- diag(rowSums(X))
L <- (D - X)
DP <- t(X) %*% D %*% X
LP <- t(X) %*% L %*% X
projectionGeigen <-  geigen(LP, DP,symmetric=TRUE)
projection_to_mean_space<-projectionGeigen$vectors[,order(projectionGeigen$values)][,-1] 
# smallest eigenvalue is numerical error, so we remove it as in 
# Vos de Wael R, Benkarim O, Paquola C, et al. BrainSpace: a toolbox for the analysis
# of macroscale gradients in neuroimaging and connectomics datasets. Commun Biol. 2020

#Compute the range for the plot so that every timepoint uses the same range
range_for_plotting<-c(range((X %*% projection_to_mean_space)[,1]),range((X %*% projection_to_mean_space)[,2]))

#colors that will be used for density maps
clusterColors<-rainbow(amount_of_clusters_to_observe)


#Project all time points to mean space
for(timepoint in 1:amount_of_timepoints){
  
  #compute affinity matrix for timepoint
  projection_at_timepoint <- getAffinityMatrix(connectivity_timeseries[[timepoint]]) %*% projection_to_mean_space
  
  for(cluster in 1:amount_of_clusters_to_observe){
   #generate the density map on a 200x200 grid
   densityMap <- kde2d(projection_at_timepoint[cluster==clusters,1], projection_at_timepoint[cluster==clusters,2], n=200, lims = range_for_plotting)

   #plot the density map
   image(densityMap, col=colorRampPalette(c('#FFFFFF00',clusterColors[cluster]),alpha=TRUE)(100),add=(cluster >1))
         
  }
  #add projected data points
  points(projection_at_timepoint[,c(1,2)],pch=".")
  title(main=paste("Timepoint",timepoint))
 
}
                     
