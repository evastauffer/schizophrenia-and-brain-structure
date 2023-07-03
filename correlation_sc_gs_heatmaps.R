#### The genetic relaionships between brain structure and schizophrenia
#### Stauffer et al., 2023

# This script correlates SC and GS, creates scatterplots and heatmaps, correlates GS and SC with geodesic distance

rm(list=ls())
options(stringsAsFactors = FALSE)

library(data.table); library(reshape); library(ggplot2);library(gtools);library(gtools);
library(ComplexHeatmap);library(egg);library(ggpubr)

### Import relevant data 
region_number = fread("~/data/raw/region_number.csv") #names and region numbers
region_number$S_no = paste0("ROI_",region_number$S_no)

### Distance, genetic similarity and structural covariance ###
distance = fread("~/data/raw/geodesic_distance_HCP.csv") # file from richard
distance = distance[1:180,1:180]
distance = as.matrix(distance)
distance_upper = as.data.frame(distance[upper.tri(distance)])

### Create colors
cmap_sns = c("#3f7f93","#7ba7b5","#b8d1d9","#f2f2f2","#e8c0b5","#d58a76","#c25539")
col_fun = colorRamp2(c(-1, 0, 1), c(cmap_sns[1],cmap_sns[4],cmap_sns[7]))


### Set up loop

# List of phenotypes
modality_list = c("SA","CT","ICVF")

# Create empty list for plots
correlation_sc_gs = NULL
correlation_gs_sc_list = NULL

distance_table_gs=NULL
distance_table_sc=NULL

distance_genetic_similarity = NULL
distance_structural_covariance = NULL

heatmap_genetics = NULL
heatmap_imaging = NULL

for (modality in modality_list) {
  #### Correlation between SC and GS  #### 
  #Import structural covariance matrix 
  cormat_imaging = fread(paste0("~/data/processed/imaging_data/cormat_imaging/structural_covariance_",modality,".txt"))
  cormat_imaging = as.matrix(cormat_imaging)

  #Import genetic similaritymatrix 
  cormat_genetics = fread(paste0("~/data/processed/cormat_genes/FB_",modality,"_cormat_spear.txt"))
  cormat_genetics = as.matrix(cormat_genetics)

  #Combine genetic and structural matrices
  structural_upper = as.data.frame(cormat_imaging[upper.tri(cormat_imaging)])
  colnames(structural_upper) = "structural_covariance"
  
  genetics_upper = as.data.frame(cormat_genetics[upper.tri(cormat_genetics )])
  colnames(genetics_upper) = "genetic_similarity"
  
  df = cbind(structural_upper, genetics_upper)

  #Save correlation results as table
  correlation = cor.test(df$structural_covariance,df$genetic_similarity, method = "spearman")
  cor_gs_sc = tidy(correlation)
  cor_gs_sc$modality = modality
  correlation_gs_sc_list = rbind(correlation_gs_sc_list,cor_gs_sc)

  #Save correlation results as plot
  correlation_sc_gs[[modality]] = ggplot(df, aes(x=structural_covariance, y=genetic_similarity))+ 
    geom_point(alpha = 0.5,size = 1.5, color = "grey") +
    geom_density2d(size=0.4)+
    stat_cor(method = "spearman",label.x.npc = "left",label.y.npc = "top",cor.coef.name = c("rho"), size =6)+
    geom_smooth(method=loess,colour = "blue")+ 
    theme_bw()+
    theme(axis.text=element_text(size=16),axis.title=element_blank())

  #### Correlation between SC, GS and distance  #### 
  df_distance = cbind(df, distance_upper)
  colnames(df_distance)= c("structural_covariance", "genetic_similarity","distance")
  
  # Correlation distance and genetic similarity
  correlation = cor.test(df_distance$distance,df_distance$genetic_similarity, method = "spearman")
  cor_dist_gs = tidy(correlation)
  cor_dist_gs$modality = modality
  distance_table_gs = rbind(distance_table_gs,cor_dist_gs)
  
  #Plot
  distance_genetic_similarity[[modality]] = ggplot(df_distance, aes(x=distance, y=genetic_similarity))+ 
    geom_point(alpha = 0.5,size = 1.5, color = "grey") +
    geom_density2d(size=0.2, color = "red")+
    geom_smooth(method=loess,colour = "red")+ 
    stat_cor(method = "spearman",label.x= 75,label.y.npc = "top",cor.coef.name = c("rho"), size =6,p.accuracy = 0.001)+
    theme_bw()+
    theme(axis.text=element_text(size=16),axis.title=element_blank())
  
  # Correlation distance and structural covariance
  correlation = cor.test(df_distance$distance,df_distance$structural_covariance, method = "spearman", exact = F)
  cor_dist_sc = tidy(correlation)
  cor_dist_sc$modality = modality
  distance_table_sc = rbind(distance_table_sc,cor_dist_sc)
  
  #Plot
  distance_structural_covariance[[modality]]=ggplot(df_distance, aes(x=distance, y=structural_covariance))+ 
    geom_point(alpha = 0.5,size = 1.5, color = "grey") +
    geom_density2d(size=0.2, color = "red")+
    geom_smooth(method=loess,colour = "red")+ 
    stat_cor(method = "spearman",label.x= 75,label.y.npc = "top",cor.coef.name = c("rho"), size =6,p.accuracy = 0.001)+
    theme_bw()+
    theme(axis.text=element_text(size=16),axis.title=element_blank())

  #### Heatmaps of genetic similarity and structural covariance #### 
  ## For genetic data
  # Arrange accoring to modules; import modules
  modules_g = fread(paste0("~/output/louvain/louvain_module_affiliation_genetics_",modality,".txt"))
  
  # Annotation for heatmap
  annot_g= HeatmapAnnotation(Module = modules_g$module,show_legend = FALSE, simple_anno_size = unit(0.25, "cm"),simple_anno_size_adjust = TRUE,
                            col=list(Module = c("1"="cyan", "2"= "red", 
                                               "3" = "gold1", "4"="green", 
                                               "5" ="blue")
                            )
  )
  cormat_genetics2 = cormat_genetics
  
  # Plot and impose structure
  axisfont = 3
  cormat_genetics2[order(modules_g$module),order(modules_g$module)] # reorder
  heatmap_genetics[[modality]] = grid.grabExpr(draw(Heatmap(cormat_genetics2, name = "correlation", col = col_fun, top_annotation= annot_g,
                                               row_title = paste0(modality),
                                               row_order =order(modules_g$module), column_order = order(modules_g$module),
                                               show_row_names = FALSE, show_column_names = FALSE, show_heatmap_legend = FALSE)))
  
  
  ## For imaging data
  # Arrange accoring to modules; import modules
  modules_i = fread(paste0("~/output/louvain/louvain_module_affiliation_imaging_",modality,".txt"))
  
  # Annotation for heatmap
  annot_i= HeatmapAnnotation(Module = modules_i$module,show_legend = FALSE, simple_anno_size = unit(0.25, "cm"),simple_anno_size_adjust = TRUE,
                             col=list(Module = c("1"="cyan", "2"= "red", 
                                                 "3" = "gold1", "4"="green", 
                                                 "5" ="blue")
                             )
  )
  cormat_imaging2 = cormat_imaging
  
  # Plot and impose structure
  axisfont = 3
  cormat_imaging2[order(modules_i$module),order(modules_i$module)] # reorder
  heatmap_imaging[[modality]] = grid.grabExpr(draw(Heatmap(cormat_imaging2, name = "correlation", col = col_fun, top_annotation= annot_i,
                                                           row_title = paste0(modality),
                                                           row_order =order(modules_i$module), column_order = order(modules_i$module),
                                                           show_row_names = FALSE, show_column_names = FALSE, show_heatmap_legend = FALSE)))
  
  
  
}  


# Export output
## Scatterplots SC and GS
correlation_all = ggpubr::ggarrange(plotlist=correlation_sc_gs, ncol =1 ,nrow =3)
ggsave(correlation_all,filename = '~/output/plots/subset/FB_correlation_gs_sc.pdf',height = 11, width = 4)

## Correlation stats SC and GS
correlation_gs_sc_list = correlation_gs_sc_list[,c(6,4,1,3)]
write.xlsx(correlation_gs_sc_list, '~/output/sc_gs_correlation_subset.xlsx')

## Distance and genetic similarity
p_distance_genetic_similarity_all = ggpubr::ggarrange(plotlist =distance_genetic_similarity, ncol =1 ,nrow =3,common.legend = F)
ggsave(p_distance_genetic_similarity_all,filename = '~/output/plots/distance_genetic_similarity_all.pdf',height = 11, width = 4)

## Distance and structural covariance
p_distance_structural_covariance_all = ggpubr::ggarrange(plotlist =distance_structural_covariance, ncol =1 ,nrow =3,common.legend = F)
ggsave(p_distance_structural_covariance_all,filename = '~/output/plots/distance_structural_covariance_all.pdf',height = 11, width = 4)

## Table distance gs and sc
distance_table_gs$pfdr = p.adjust(distance_table_gs$p.value, method = "fdr")
distance_table_sc$pfdr = p.adjust(distance_table_sc$p.value, method = "fdr")
table_distance = rbind(distance_table_gs,distance_table_sc)
table_distance  = table_distance [,c(6,1,3,7)]
write.xlsx(table_distance ,'~/output/distance_table_gs_sc.xlsx')

## Heatmaps genetic similarity
heatmap_genetics_all = ggpubr::ggarrange(
  ggpubr::ggarrange(plotlist =heatmap_genetics, ncol =1 ,nrow =4,common.legend = F))
ggsave(heatmap_genetics_all,filename = '~/output/plots/subset/FB_heatmap.png',height = 11, width = 4)

## Heatmaps structural covariance
heatmap_imaging_all = ggpubr::ggarrange(
  ggpubr::ggarrange(plotlist =heatmap_imaging, ncol =1 ,nrow =4,common.legend = F))
ggsave(heatmap_imaging_all,filename = '~/output/plots/subset/heatmap_imaging.png',height = 11, width = 4)



