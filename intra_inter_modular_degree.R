#### The genetic relaionships between brain structure and schizophrenia
#### Stauffer et al., 2023

# This script calculates participation coefficients. Once within module degree centrality and once between module degree centrality
# scatterplot between pls and participation

library(data.table);library(dplyr);library(ggplot2);
library(ggpubr);library(igraph);library(tibble);library(brainGraph);
library(ggseg); library(ggsegGlasser);library(gtools);library(network);
library(ggnet); library(broom);library(cmocean);library(xlsx)


modality_list = c("SA", "CT", "ICVF")

# Prepare some lists to store output
within_all = NULL
between_all = NULL

# Set up loop
for (modality in modality_list) {
  #Import structural covariance matrix 
  cormat_imaging = fread(paste0("~/data/processed/imaging_data/cormat_imaging/structural_covariance_",modality,".txt"))
  cormat_p = cormat_imaging                # Duplicate data frame
  cormat_p[cormat_p < 0] = 0   
  cormat_imaging2 = as.matrix(cormat_p)
  
  modules_i = fread(paste0("~/output/louvain/louvain_module_affiliation_imaging_",modality,".txt"))
  module_vec = c(modules_i$module)
  cormat_imaging2=as.data.frame(cormat_imaging2)
  
  c = unique(modules_i$module)
  
  # Calculate within and between degree for each module
  within = between = vector(length = 180)
  
  for (i in c) {
    within[which(module_vec==i)] = rowMeans(cormat_imaging2[which(module_vec==i),which(module_vec==i)])
    between[which(module_vec==i)] = rowMeans(cormat_imaging2[which(module_vec==i),which(module_vec!=i)])
  }
  
  within = as.data.frame(within)
  within$module = module_vec
  within$ROI = names(cormat_imaging2)

  
  ### add PLS weights for component 1 ###
  pls_hub = fread(paste0("~/output/SCZ_",modality,"_pls_weight_hubs.txt"))
  within = merge(within, pls_hub, by.x = "ROI", by.y = "name")
  within_all = rbind(within_all,within)
  
  between = as.data.frame(between)
  between$module = module_vec
  between$ROI = names(cormat_imaging2)
  between = merge(between, pls_hub, by.x = "ROI", by.y = "name")
  between_all = rbind(between_all,between)
  

}


#### Visualise output
within_all$modality  = factor(within_all$modality,levels = c("SA","CT", "ICVF"))

## Plot pls weights and within connectivity
p_pls_within_i= ggplot(within_all , aes(x=boot.weight, y = within))+ 
  geom_point(size = 0.8) +
  geom_smooth(method=lm,colour = "red")+ 
  facet_wrap(~ modality, ncol = 1,nrow = 3) +
  theme_bw()+
  stat_cor(method = "spearman",label.x.npc = "left",label.y.npc = "top",cor.coef.name = c("rho"), size =5,aes(label = ..r.label..))+
  theme(axis.text=element_text(size=14),axis.title=element_blank(),plot.title = element_blank(),strip.background = element_blank(),strip.text.x = element_blank())

ggsave(p_pls_within_i,filename = "~/output/plots/subset/SCZ_scatterplot_plsweight_within_i.pdf",height = 4.5, width = 3)

## Plot pls weights and between connectivity
p_pls_within_i= ggplot(between_all , aes(x=boot.weight, y=between))+ 
  geom_point(size = 0.8) +
  geom_smooth(method=lm,colour = "red")+ 
  facet_wrap(~ modality, ncol = 1,nrow = 3) +
  theme_bw()+
  stat_cor(method = "spearman",label.x.npc = "left",label.y.npc = "top",cor.coef.name = c("rho"), size =5,aes(label = ..r.label..))+
  theme(axis.text=element_text(size=14),axis.title=element_blank(),plot.title = element_blank(),strip.background = element_blank(),strip.text.x = element_blank())

ggsave(p_pls_within_i,filename = "~/output/plots/subset/SCZ_scatterplot_plsweight_between_i.pdf",height = 4.5, width = 3)




