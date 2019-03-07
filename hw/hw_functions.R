# plot 5 curves in 1 plot
library(stats)
library(MASS)
library(Matrix)
library(scatterplot3d)
library(rerf)
library(vegan)
source("precision_recall.R")
library(umap)

normalizeData <- function(X) {
  X <- sweep(X, 2, apply(X, 2, min), "-")
  sweep(X, 2, apply(X, 2, max), "/")
}


hw_data<-function(num_of_points=500){
  set.seed(1)
  t=runif(num_of_points, min = 0, max = 1)
  t=sort(t)
  x1=t^2
  x2=2*t*(1-t)
  x3=(1-t)^2
  data = cbind(x1, x2, x3)
  return(list(data, t))
}



# generate D_geo, this is the geodesic distance
hw_geodesic<-function(t, num_of_points){
  f = function(x){sqrt(8*(3*x^2-3*x+1))}
  D_geo = matrix(rep(0, num_of_points*num_of_points), nrow = num_of_points, ncol = num_of_points)
  j=1
  while (j <= num_of_points){
    i=1
    while( i <= num_of_points){
      D_geo[i,j] = abs(integrate(f,t[i],t[j])$value)
      i = i +1
    }
    j=j+1
  }
  return(D_geo)
}

# #### We generate the first noise 

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

generate_high_dim_uniform_noise<-function(num_of_points, noise_dim=6, const=70){
  #matrix_of_0 = matrix(rep(0, num_of_points*(noise_dim)), nrow = num_of_points, ncol = noise_dim)
  #cov_matrix = matrix(rep(0, noise_dim*noise_dim), nrow = noise_dim, ncol = noise_dim)
  #diag(cov_matrix) = c(rep(const, noise_dim))
  #Sig1 = cov_matrix
  #noise = mvrnorm(n = num_of_points, (rep(0, noise_dim)), Sig1, tol = 1e-7, empirical = FALSE, EISPACK = FALSE)
  noise=replicate(noise_dim, runif(num_of_points, (-1)*const, const)) 
  rownames(noise) <- c()
  colnames(noise) <- c()
  return(noise)
}

