#### The genetic relaionships between brain structure and schizophrenia
#### Stauffer et al., 2023

# This script performs enrichments of PLS weights in Mesulam classes


rm(list=ls())
options(stringsAsFactors = FALSE)

library(ggplot2);library(data.table);library(tidyr);library(tidyverse);library(xlsx);
library(ggridges); library(httr); library(dplyr);library(matrixStats);require(reshape2)

### Preparation
# Set path
inpath = '~/output/pls/' # Path to output directory

# mapping of glasser regions
mapping = fread("~/data/raw/mapping_glasser.txt")

# custom colormap for mesulam
cmap=c("#f8a796","#6182ac", "#7ca840","#fbd779")

# Load HCP to mesulam
hcp2mesulam <- read.csv("https://raw.githubusercontent.com/ucam-department-of-psychiatry/maps_and_parcs/master/Maps/map2map_hcp_to_mesulam.csv", header=T)
hcp2mesulam = hcp2mesulam[2:181,]
hcp2mesulam$label1 = gsub("^L_","",hcp2mesulam$label1)
hcp2mesulam$label1 = gsub("_ROI","",hcp2mesulam$label1)

#adapt order
mapping = mapping[match(hcp2mesulam$label1,mapping$region ),]

#load all things for running spin- permutation
coord <- read.table("https://raw.githubusercontent.com/rb643/rotate_parcellation/master/sphere_HCP.txt", header=F)
source('https://raw.githubusercontent.com/rb643/rotate_parcellation/master/R/rotate.parcellation.R')
perms <- rotate.parcellation(coord.l = as.matrix(coord[1:180,]), coord.r = as.matrix(coord[181:360,]), nrot = 1000)
perms = read.table("~/data/processed/permutation.txt")

### Set up Loop 
modality_list = c("SA","CT","ICVF")
  enrichmes_all = NULL
  stats_all = NULL
  for (modality in modality_list) {
    ### Import PLS weights for component 1 and 2 ###
    pls_weights_component_1 = fread(paste(inpath,"SCZ_",modality,"_PLS_Weights_Component_1_ascending.csv", sep = ""))
    pls_weights_component_1$component = "Component 1"
    
    # Combine pls weights with Mesulam mappingg
    weights_mapping = left_join(mapping,pls_weights_component_1, by= c("S_no"= "name"))
    weights_mapping1 = weights_mapping2 = as.data.frame(weights_mapping)
    weights_mapping1$name = paste0("lh_L_", weights_mapping1$region)
    weights_mapping2$name = paste0("rh_R_", weights_mapping2$region)
    weights_mapping3= rbind(weights_mapping1,weights_mapping2)
    
    # Calculate mean PLS weight per class
    real = weights_mapping3 %>% group_by(Class) %>% summarise(mean_weight = mean(boot.weight))
    
    # populate the null model by looping through the permuted indices and recomputing the mean
    null <- real
    colnames(null) =  c("Class","Real")
    for (i in 1:1000){
      tempweights = weights_mapping3
      tempweights$Class =  tempweights$Class[perms[,i]]
      tempnull =  tempweights %>% group_by(Class) %>% summarise(mean_tstat = mean(boot.weight))
      null =  merge(null,tempnull,by='Class')
    }
    
    # need some reshaping for plotting
    null$Real = NULL
    null =  t(as.matrix(null))
    colnames(null) =  null[1,]
    null =  null[-1,]
    null =  melt(null)
    colnames(real) =  c("Var2","realvalue")
    null =  merge(null,real,by='Var2',no.dups = F)
    null$value =  as.numeric(as.character(null$value))
    null$modality = modality

    mu = null %>% group_by(Var2) %>% summarise(meanV = mean(value))
    std =  null %>% group_by(Var2) %>% summarise(sdV = sd(value))
    x =  null %>% group_by(Var2) %>% summarise(x = mean(realvalue))
    
    z =  left_join(mu,std,by = "Var2")
    z =  left_join(z,x,by = "Var2")
    z$z =  (z$x - z$meanV)/z$sdV
    z$p =  (1-pnorm(abs(z$z)))*2
    z = z[z$Var2 %in% c("heteromodal","idiotypic","paralimbic","unimodal"),]
    z$pfdr = p.adjust(z$p, method = "fdr")
    z$modality =  modality
    
    stats_all =  rbind(stats_all,z)  
    enrichmes_all =  rbind(enrichmes_all,null)
  }
  ### Export results
  write.xlsx2(x = as.data.frame(stats_all), file= paste(inpath, "enrichment_mesulam_subset_SCZ.xlsx",sep = ""), row.names = F)
  
  ### Get significant enrichments
  sig = as.data.frame(stats_all[stats_all$pfdr<0.05,])
  write.xlsx2(x = sig, file= paste(inpath, "enrichment_mesulam_SCZ.xlsx",sep = ""), row.names = F)

  ### Plot distributions
  enrichmes_all$modality <- as.factor(enrichmes_all$modality)
  enrichmes_all <- enrichmes_all[enrichmes_all$Var2 %in% c("heteromodal","idiotypic","paralimbic","unimodal"),]
  enrichmes_all<- droplevels(enrichmes_all)
  
  # null distribution as ridgeplot and real values on top.
  enrichments_all_plot = enrichmes_all
  levels(enrichments_all_plot$modality)[match("ICVF",levels(enrichments_all_plot$modality))] <-"NDI"
  
  sig$modality = as.factor(sig$modality)
  levels(sig$modality)[match("ICVF",levels(sig$modality))] <-"NDI"
  
  enrichments_all_plot = left_join(enrichments_all_plot, sig, by = c("modality"="modality", "Var2"="Var2"))
  enrichments_all_plot$ast = add.significance.stars(enrichments_all_plot$pfdr, cutoffs = c(0.05, 0.01, 0.001))
  
  p_mesulam <- ggplot(data=enrichments_all_plot, aes(y=Var2,x=value,fill=Var2)) +
    geom_density_ridges(alpha = 0.2) +
    geom_point(aes(y=Var2,x=realvalue, fill=Var2),size=3,shape=21) +
    scale_fill_manual(values = cmap) + 
    scale_colour_manual(values = cmap) + 
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      strip.background = element_blank(),
      panel.spacing = unit(2, "lines")
    ) +  
    facet_wrap(~modality,nrow = 3, ncol = 3, scales ="free" )+
    geom_text(aes(label = ast, x = realvalue + 0.5), size =4.5, angle = 90)+
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          title = element_blank(),axis.text=element_text(size=16),strip.text = element_text(size=14)) +
    coord_flip() 
  ggsave( mesulam ,filename = "~/output/plots/subset/SCZ_PLS_enrichment_mesulam_sig.pdf", sep = "",height = 4, width = 6)
  





