#upper sphere
library(scatterplot3d)
library(rgl)
library(rerf)
library(umap)
library(vegan)
library(pracma)
library(MASS)

#this generate ellipsoid data
upper_sphere<-function(r=3, N=200, plot=TRUE){
  
  set.seed(1)
  data=c()
  
  parameters=c()
  
  index=1
  while(index<=N){
    u_v=c(runif(1, 0, pi), runif(1, 0, pi))
    parameters=rbind(parameters, u_v)
    index=index+1
  }
  rownames(parameters) <- c()
  colnames(parameters) <- c()
  
  data=c()
  row=1
  while (row<=dim(parameters)[1]){
    x = r*cos(parameters[row, 1])*sin(parameters[row,2])
    y = r*sin(parameters[row, 1])*sin(parameters[row,2])
    z = r*cos(parameters[row,2])
    data1=cbind(cbind(x,y), z)
    data=rbind(data, data1)
    row=row+1
  }
  
  #data=(data[order(data[,1]),])
  if (plot){
    plot3d(data[1:(num_of_points/5),1], data[1:(num_of_points/5),2], data[1:(num_of_points/5),3], col='red', size=3, xlim = c(-r-1,r+1), ylim=c(-r-1,r+1), zlim = c(-r-1, r+1))
    plot3d(data[(num_of_points/5+1):(2*num_of_points/5),1], data[(num_of_points/5+1):(2*num_of_points/5),2], data[(num_of_points/5+1):(2*num_of_points/5),3], col='green', size=3, add=TRUE)
    plot3d(data[(2*num_of_points/5+1):(3*num_of_points/5),1], data[(2*num_of_points/5+1):(3*num_of_points/5),2], data[(2*num_of_points/5+1):(3*num_of_points/5),3], col='orange', size=3, add=TRUE)
    plot3d(data[(3*num_of_points/5+1):(4*num_of_points/5),1], data[(3*num_of_points/5+1):(4*num_of_points/5),2], data[(3*num_of_points/5+1):(4*num_of_points/5),3], col='blue', size=3, add=TRUE)
    plot3d(data[(4*num_of_points/5+1):(num_of_points),1], data[(4*num_of_points/5+1):(num_of_points),2], data[(4*num_of_points/5+1):(num_of_points),3], size=3, add=TRUE)
  }
  return(list(data, parameters))
}

#we will generate the geodesic distance
f = function(data1, data2, r){ 
  value=(t(data1)%*%data2/(r^2))
  phi= acos(pmin(pmax(value,-1.0),1.0))
  d=r*phi
  return(d)
}

sphere_geodesic<-function(data, r, num_of_points){
  D_geo = matrix(rep(0, num_of_points*num_of_points), nrow = num_of_points, ncol = num_of_points)
  j=1
  while (j <= num_of_points){
    i=1
    while( i <= num_of_points){
      D_geo[i,j] = f(data[i,], data[j,], r)
      i = i +1
    }
    j=j+1
  }
  diag(D_geo)=0
  return(D_geo)
}

#this generate noise
generate_high_dim_gaussian_noise<-function(num_of_points, noise_dim=6, const=70){
  matrix_of_0 = matrix(rep(0, num_of_points*(noise_dim)), nrow = num_of_points, ncol = noise_dim)
  cov_matrix = matrix(rep(0, noise_dim*noise_dim), nrow = noise_dim, ncol = noise_dim)
  diag(cov_matrix) = c(rep(const, noise_dim))
  Sig1 = cov_matrix
  noise = mvrnorm(n = num_of_points, (rep(0, noise_dim)), Sig1, tol = 1e-7, empirical = FALSE, EISPACK = FALSE)
  rownames(noise) <- c()
  colnames(noise) <- c()
  return(noise)
}

normalizeData <- function(X) {
  X <- sweep(X, 2, apply(X, 2, min), "-")
  sweep(X, 2, apply(X, 2, max), "/")
}
