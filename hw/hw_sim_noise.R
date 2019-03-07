#hw_sim_noise
source("hw_functions.R")
source("precision_recall.R")
library(stats)
library(MASS)
library(Matrix)
library(scatterplot3d)
library(rerf)
library(vegan)
library(umap)

num_of_points=1000

hw_data_object=hw_data(num_of_points)
hw_data=hw_data_object[[1]]
t=hw_data_object[[2]]

D_truth=hw_geodesic(t, num_of_points)

#generate noise
noise_dim=2
#noise=generate_high_dim_gaussian_noise(num_of_points, high_dim, const = 70)
noise_1=generate_high_dim_uniform_noise(num_of_points, noise_dim, const = 70)
noise_2=generate_high_dim_gaussian_noise(num_of_points, noise_dim, const=70)


hw_noise_data=cbind(hw_data, noise_2)
#non_normalize vs normalize (without noise)
#first we do non_normalized
#URerF
g_noise=Urerf(hw_noise_data, trees = 300, Progress = TRUE, LinearCombo=FALSE, splitCrit = "bicfast", normalizeData = FALSE)
W_noise=g_noise$similarityMatrix
D_rf_noise=1-W_noise
# isomap
D_eucd_noise = as.matrix(dist(hw_noise_data))
iso_dist_noise = as.matrix(isomapdist(D_eucd_noise, k=9))
#UMAP
custom.settings = umap.defaults
custom.settings$n_neighbors=num_of_points
a_noise = umap(hw_noise_data, config = custom.settings)
D_umap_noise=as.matrix(dist(a_noise$layout))

############################################################################
#now we do normalized
norm_hw_noise_data=normalizeData(hw_noise_data)
#URerF norm
g_noise_norm=Urerf(norm_hw_noise_data, trees = 300, Progress = TRUE, splitCrit = "bicfast", normalizeData = FALSE)
W_noise_norm=g_noise_norm$similarityMatrix
D_rf_noise_norm=1-W_noise_norm
# isomap norm
D_eucd_noise_norm = as.matrix(dist(norm_hw_noise_data))
iso_dist_noise_norm = as.matrix(isomapdist(D_eucd_noise_norm, k=9))
#UMAP norm
a_noise_norm = umap(norm_hw_noise_data, config = custom.settings)
D_umap_noise_norm=as.matrix(dist(a_noise_norm$layout))

at_K=seq(50, 300, 50)

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
     at_K, file="hw_noise_precision_list.Rdata")

