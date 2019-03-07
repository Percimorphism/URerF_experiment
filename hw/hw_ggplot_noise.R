library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
load("hw_noise_precision_list.Rdata")


precision_df1=data.frame(at_K, prec=D_rf_noise_precision_list, Algo = as.factor("URerF"), Norm =as.factor("no normalization"))
precision_df2=data.frame(at_K, prec=D_iso_noise_precision_list, Algo = as.factor("Isomap"), Norm =as.factor("no normalization"))
precision_df3=data.frame(at_K, prec=D_umap_noise_precision_list, Algo = as.factor("umap"), Norm =as.factor("no normalization"))
precision_df4=data.frame(at_K, prec=D_rf_noise_norm_precision_list, Algo = as.factor("URerF"), Norm =as.factor("normalization"))
precision_df5=data.frame(at_K, prec=D_iso_noise_norm_precision_list, Algo = as.factor("Isomap"), Norm =as.factor("normalization"))
precision_df6=data.frame(at_K, prec=D_umap_noise_norm_precision_list, Algo = as.factor("umap"), Norm =as.factor("no normalization"))

precision_df=rbind(precision_df1, precision_df2, precision_df3, precision_df4, precision_df5, precision_df6)
p <- ggplot(precision_df, aes(at_K, prec, colour = Algo, shape = Norm)) + geom_line(alpha=0.9, show.legend = TRUE) + geom_point(alpha=0.9, show.legend = TRUE)  + xlab('@K') + ylab('Precision/Recall') + scale_color_brewer(palette="Dark2")
plot(p)                