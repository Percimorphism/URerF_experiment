# helix simulations
source("helix_functions.R")
source("precision_recall.R")
library(stats)
library(MASS)
library(Matrix)
library(scatterplot3d)
library(rerf)
library(vegan)
library(umap)

num_of_points=1000

helix_data_object=helix_data(num_of_points)
#plot3d(helix_data[,1], helix_data[,2], helix_data[,3])
helix_data=helix_data_object[[1]]
t=helix_data_object[[2]]

D_truth=helix_geodesic(t, num_of_points)

#non_normalize vs normalize (without noise)
#first we do non_normalized
#URerF
g=Urerf(helix_data, trees = 300, Progress = TRUE, splitCrit = "bicfast", normalizeData = FALSE)
W=g$similarityMatrix
D_rf=1-W
# isomap
D_eucd = as.matrix(dist(helix_data))
iso_dist = as.matrix(isomapdist(D_eucd, k=9))
#UMAP
custom.settings = umap.defaults
custom.settings$n_neighbors=num_of_points
a = umap(helix_data, config = custom.settings)
D_umap=as.matrix(dist(a$layout))

############################################################################
#now we do normalized
norm_helix_data=normalizeData(helix_data)
#URerF norm
g_norm=Urerf(norm_helix_data, trees = 300, Progress = TRUE, splitCrit = "bicfast", normalizeData = FALSE)
W_norm=g_norm$similarityMatrix
D_rf_norm=1-W_norm
# isomap norm
D_eucd_norm = as.matrix(dist(norm_helix_data))
iso_dist_norm = as.matrix(isomapdist(D_eucd_norm, k=9))
#UMAP norm
a_norm = umap(norm_helix_data, config = custom.settings)
D_umap_norm=as.matrix(dist(a_norm$layout))

at_K=seq(50, 300, 50)

#########################actually generate the p-r list#################################
D_rf_p_r_list = p_r_list(D_rf, D_truth, at_K, num_of_points)
D_rf_precision_list= D_rf_p_r_list$precisionList

#####################################################################################
D_iso_p_r_list = p_r_list(iso_dist, D_truth, at_K, num_of_points)
D_iso_precision_list= D_iso_p_r_list$precisionList

####################################################################
D_umap_p_r_list = p_r_list(D_umap, D_truth, at_K, num_of_points)
D_umap_precision_list= D_umap_p_r_list$precisionList


######normalized p################
#########################actually generate the p-r list#################################
D_rf_norm_p_r_list = p_r_list(D_rf_norm, D_truth, at_K, num_of_points)
D_rf_norm_precision_list= D_rf_norm_p_r_list$precisionList

#####################################################################################
D_iso_norm_p_r_list = p_r_list(iso_dist_norm, D_truth, at_K, num_of_points)
D_iso_norm_precision_list= D_iso_norm_p_r_list$precisionList

####################################################################
D_umap_norm_p_r_list = p_r_list(D_umap_norm, D_truth, at_K, num_of_points)
D_umap_norm_precision_list= D_umap_p_r_list$precisionList

save(D_rf_precision_list, D_rf_norm_precision_list, 
     D_iso_precision_list, D_iso_norm_precision_list,
     D_umap_precision_list, D_umap_norm_precision_list, 
     at_K, file="helix_no_noise_precision_list.Rdata")
