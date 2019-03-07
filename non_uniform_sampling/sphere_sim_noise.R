# sphere noise simulation
source("sphere_function.R")
source("precision_recall.R")
library(stats)
library(MASS)
library(Matrix)
library(scatterplot3d)
library(rerf)
library(vegan)
library(umap)

num_of_points=2000
sphere_data_object=upper_sphere(N=num_of_points)
sphere_data=sphere_data_object[[1]]
parameters=sphere_data_object[[2]]

D_truth=sphere_geodesic(sphere_data, r=3, num_of_points)


#generate noise
noise_dim=6
noise=generate_high_dim_gaussian_noise(num_of_points, noise_dim, const = 70)


sphere_noise_data=cbind(sphere_data, noise)
#non_normalize vs normalize (without noise)
#first we do non_normalized
#URerF
g_noise=Urerf(sphere_noise_data, trees = 300, Progress = TRUE, LinearCombo=FALSE, splitCrit = "bicfast", normalizeData = FALSE)
W_noise=g_noise$similarityMatrix
D_rf_noise=1-W_noise
# isomap
D_eucd_noise = as.matrix(dist(sphere_noise_data))
iso_dist_noise = as.matrix(isomapdist(D_eucd_noise, k=9))
#UMAP
custom.settings = umap.defaults
custom.settings$n_neighbors=num_of_points
a_noise = umap(sphere_noise_data, config = custom.settings)
D_umap_noise=as.matrix(dist(a_noise$layout))

############################################################################
#now we do normalized
norm_sphere_noise_data=normalizeData(sphere_noise_data)
#URerF norm
g_noise_norm=Urerf(norm_sphere_noise_data, trees = 300, Progress = TRUE, splitCrit = "bicfast", normalizeData = FALSE)
W_noise_norm=g_noise_norm$similarityMatrix
D_rf_noise_norm=1-W_noise_norm
# isomap norm
D_eucd_noise_norm = as.matrix(dist(norm_sphere_noise_data))
iso_dist_noise_norm = as.matrix(isomapdist(D_eucd_noise_norm, k=9))
#UMAP norm
a_noise_norm = umap(norm_sphere_noise_data, config = custom.settings)
D_umap_noise_norm=as.matrix(dist(a_noise_norm$layout))

at_K=seq(10, 300, 30)

#########################actually generate the p-r list#################################
D_rf_noise_p_r_list = p_r_list(D_rf_noise, D_truth, at_K, num_of_points)
D_rf_noise_precision_list= D_rf_noise_p_r_list$precisionList

#####################################################################################
D_iso_noise_p_r_list = p_r_list(iso_dist_noise, D_truth, at_K, num_of_points)
D_iso_noise_precision_list= D_iso_noise_p_r_list$precisionList

####################################################################
D_umap_noise_p_r_list = p_r_list(D_umap_noise, D_truth, at_K, num_of_points)
D_umap_noise_precision_list= D_umap_noise_p_r_list$precisionList


######normalized p################
#########################actually generate the p-r list#################################
D_rf_noise_norm_p_r_list = p_r_list(D_rf_noise_norm, D_truth, at_K, num_of_points)
D_rf_noise_norm_precision_list= D_rf_noise_norm_p_r_list$precisionList

#####################################################################################
D_iso_noise_norm_p_r_list = p_r_list(iso_dist_noise_norm, D_truth, at_K, num_of_points)
D_iso_noise_norm_precision_list= D_iso_noise_norm_p_r_list$precisionList

####################################################################
D_umap_noise_norm_p_r_list = p_r_list(D_umap_noise_norm, D_truth, at_K, num_of_points)
D_umap_noise_norm_precision_list= D_umap_noise_p_r_list$precisionList

save(D_rf_noise_precision_list, D_rf_noise_norm_precision_list, 
     D_iso_noise_precision_list, D_iso_noise_norm_precision_list,
     D_umap_noise_precision_list, D_umap_noise_norm_precision_list, 
     at_K, file="sphere_noise_precision_list.Rdata")

