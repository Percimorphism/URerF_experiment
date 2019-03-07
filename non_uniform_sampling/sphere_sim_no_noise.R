# sphere no noise simulation
# this is no longer uniform sampling from a sphere
# sampling 4 patches on the surface and see how that works
source("sphere_function.R")
source("old_precision_recall.R")
library(stats)
library(MASS)
library(Matrix)
library(scatterplot3d)
library(rerf)
library(vegan)
library(umap)

N=200
num_of_points=4*N
sphere_data_object=non_uniform_sampleing_upper_sphere(N=N)
sphere_data=sphere_data_object[[1]]
parameters=sphere_data_object[[2]]

#D_truth=sphere_geodesic(sphere_data, r=3, num_of_points)
data_label= c(rep('1', num_of_points/4), rep('2', num_of_points/4), rep('3', num_of_points/4), rep('4', num_of_points/4))


#non_normalize vs normalize (without noise)
#first we do non_normalized
#URerF
g=Urerf(sphere_data, trees = 300, Progress = TRUE, splitCrit = "bicfast", normalizeData = FALSE)#, min.parent = 20)
W=g$similarityMatrix
D_rf=1-W
# isomap
D_eucd = as.matrix(dist(sphere_data))
iso_dist = as.matrix(isomapdist(D_eucd, k=100))
#UMAP
custom.settings = umap.defaults
custom.settings$n_neighbors=num_of_points
a = umap(sphere_data, config = custom.settings)
D_umap=as.matrix(dist(a$layout))

# ############################################################################
# #now we do normalized
# norm_sphere_data=normalizeData(sphere_data)
# #URerF norm
# g_norm=Urerf(norm_sphere_data, trees = 300, Progress = TRUE, splitCrit = "bicfast", normalizeData = FALSE)
# W_norm=g_norm$similarityMatrix
# D_rf_norm=1-W_norm
# # isomap norm
# D_eucd_norm = as.matrix(dist(norm_sphere_data))
# iso_dist_norm = as.matrix(isomapdist(D_eucd_norm, k=20))
# #UMAP norm
# a_norm = umap(norm_sphere_data, config = custom.settings)
# D_umap_norm=as.matrix(dist(a_norm$layout))

at_K=seq(1, 100, 5)

#########################actually generate the p-r list#################################
D_rf_p_r_list = p_r_list(D_rf, data_label, at_K, num_of_points)
D_rf_precision_list= D_rf_p_r_list$precisionList

#####################################################################################
D_iso_p_r_list = p_r_list(iso_dist, data_label, at_K, num_of_points)
D_iso_precision_list= D_iso_p_r_list$precisionList

####################################################################
D_umap_p_r_list = p_r_list(D_umap, data_label, at_K, num_of_points)
D_umap_precision_list= D_umap_p_r_list$precisionList


######normalized p################
#########################actually generate the p-r list#################################
# D_rf_norm_p_r_list = p_r_list(D_rf_norm, D_truth, at_K, num_of_points)
# D_rf_norm_precision_list= D_rf_norm_p_r_list$precisionList
# 
# #####################################################################################
# D_iso_norm_p_r_list = p_r_list(iso_dist_norm, D_truth, at_K, num_of_points)
# D_iso_norm_precision_list= D_iso_norm_p_r_list$precisionList
# 
# ####################################################################
# D_umap_norm_p_r_list = p_r_list(D_umap_norm, D_truth, at_K, num_of_points)
# D_umap_norm_precision_list= D_umap_p_r_list$precisionList

save(D_rf_precision_list, 
     D_iso_precision_list, 
     D_umap_precision_list,
     at_K, file="sphere_no_noise_precision_list.Rdata")
