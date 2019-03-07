# plot 5 curves in 1 plot
library(stats)
library(MASS)
library(Matrix)
library(scatterplot3d)
library(rerf)
library(vegan)
library(umap)

normalizeData <- function(X) {
  X <- sweep(X, 2, apply(X, 2, min), "-")
  sweep(X, 2, apply(X, 2, max), "/")
}


helix_data<-function(num_of_points=500){
  set.seed(1)
  p = (3 * pi / 2) * (1 + 2*sort(runif(num_of_points, 0, 1)));  
  samples = as.matrix(cbind(cbind(2*p*cos(2*p), 2*p*sin(2*p)), 2*p))
  return(list(samples, p))
}



# generate D_geo, this is the geodesic distance for helix
helix_geodesic<-function(t, num_of_points){
  f = function(p){sqrt( (2*cos(2*p)-4*p*sin(2*p))^2 + (2*sin(2*p) + 4*p*cos(2*p))^2  + 4)}
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


