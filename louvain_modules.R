#### The genetic relaionships between brain structure and schizophrenia
#### Stauffer et al., 2023

# This script identified modules (Louvain) and plots modules onto brain 

## Load package
library(igraph); library(data.table); library(gtools);
library(ggseg);library(ggsegGlasser);library(ggsegExtra);
library(RColorBrewer); library(ggplot2); library(dplyr);
library(grid);library(ggnet);library(xlsx);library(reshape2);
library(dplyr);library(tibble)

# Import relevant data
mapping = fread("/mnt/b2/home4/arc/ems206/cortex_project_04_03_2021/data/raw/mapping_glasser.txt")
mapping = mapping[mixedorder(mapping$S_no),]


#### Louvain clustering ####

#Set up loop
modality_list = c("SA","CT","ICVF")

# Make some storage
louvain_communities = NULL
louvain_imaging_list = NULL
louvain_genetics = NULL
set.seed(34)

for (modality in modality_list) {
  
  ##### Networks for genetic data
  #Read in genetic matrix
  cormat_genes = fread(paste0("~/data/processed/cormat_genes/FB_",modality,"_cormat_spear.txt"))
  cormat_genes = as.matrix(cormat_genes)
  cormat_genes_t =cormat_genes
  diag(cormat_genes)<-0
  
  ## Modules for unthresholded matrix 
  # Create network
  net<-graph.adjacency(abs(cormat_genes),weighted=TRUE,mode="undirected")
  
  # Louvain clusters
  louvain = cluster_louvain(net, weights = NULL)
  V(net)$community <- louvain$membership
  louvain_communities[[modality]] = louvain
  
  # save louvain membership
  community = data.frame("S_no" = louvain$names,"module" = louvain$membership)
  write.table(community,paste0("~/output/louvain/louvain_module_affiliation_genetics_",modality,".txt"), row.names = F, col.names = T, quote = F)
  community$modality = modality
  community = merge(community,mapping, by = "S_no")
  louvain_genetics = rbind(louvain_genetics,community)
  
  ##### Networks for imaging data
  #Read in genetic matrix
  cormat_imaging = fread(paste0("~/data/processed/imaging_data/cormat_imaging/structural_covariance_",modality,".txt"))
  cormat_imaging = as.matrix(cormat_imaging)
  diag(cormat_imaging)<-0
  
  ## Modules for unthresholded matrix 
  # Create network
  net_imaging<-graph.adjacency(abs(cormat_imaging),weighted=TRUE,mode="undirected")
  
  # Louvain clusters
  louvain_imaging = cluster_louvain(net_imaging, weights = NULL)
  V(net_imaging)$community <- louvain_imaging$membership
  
  # save louvain membership
  community_imaging = data.frame("S_no" = louvain_imaging$names,"module" = louvain_imaging$membership)
  write.table(community_imaging ,paste0("~/output/louvain/louvain_module_affiliation_imaging_",modality,".txt"), row.names = F, col.names = T, quote = F)
  community_imaging$modality = modality
  community_imaging = merge(community_imaging,mapping, by = "S_no")
  louvain_imaging_list = rbind(louvain_imaging_list,community_imaging)
  
}

### Plot modules on brain ###
# Plot genetic modules
louvain_genetics$label <- paste0('rh_R_',louvain_genetics$region)
louvain_genetics$module = as.factor(louvain_genetics$module)
louvain_genetics$modality = as.factor(louvain_genetics$modality)
louvain_genetics$modality  = factor(louvain_genetics$modality,levels = c("SA","CT", "ICVF"))

p_genes = louvain_genetics%>% 
  group_by(modality)%>% 
  ggplot() +
  geom_brain(atlas = glasser, hemi="right",aes(fill = module), colour="black")+ 
  labs(fill = "Module")+
  scale_fill_manual(values = c("1"="cyan", "2"= "red", 
                               "3" = "gold1", "4"="green", 
                               "5" ="blue"),na.value = "transparent")+ 
  theme_void()+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),legend.position = "bottom", legend.title =element_blank())+
  facet_wrap("modality", ncol = 1)
ggsave(p_genes,filename = "~/output/plots//louvain_clusters_genes.pdf",height = 6, width =4)

# Plot imaging modules
louvain_imaging_list$label <- paste0('rh_R_',louvain_imaging_list$region)
louvain_imaging_list$module = as.factor(louvain_imaging_list$module)
louvain_imaging_list$modality = as.factor(louvain_imaging_list$modality)
louvain_imaging_list$modality  = factor(louvain_imaging_list$modality,levels = c("SA","CT", "ICVF"))

p_imaging = louvain_imaging_list%>% 
  group_by(modality)%>% 
  ggplot() +
  geom_brain(atlas = glasser, hemi="right",aes(fill = module), colour="black")+ 
  labs(fill = "Module")+
  scale_fill_manual(values = c("1"="cyan", "2"= "red", 
                               "3" = "gold1", "4"="green", 
                               "5" ="blue"),na.value = "transparent")+ 
  theme_void()+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),legend.position = "bottom", legend.title =element_blank())+
  facet_wrap("modality", ncol = 1)
ggsave(p_imaging,filename = "~/output/plots//louvain_clusters_imaging.pdf",height = 6, width =4)
