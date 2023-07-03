#### The genetic relaionships between brain structure and schizophrenia
#### Stauffer et al., 2023

# This script plots PLS weights on brain and extracts t and u scores for genes 

# 1) Loop through component 1 for each modality and plot weights on brain
# 2) Get number of significant regions per modality
# 3) Get mean PLS weight per modality
# 4) Export t and u components for later

rm(list=ls())
options(stringsAsFactors = FALSE)

library(plsdepot);library(ggplot2);library(svMisc);library(data.table);
library(tidyr);library(tidyverse);library(xlsx); library(ggplot2);
library(ggplot2);library(ggseg);library(ggsegGlasser);library(scico);
library(ggrepel); library(Hmisc)

######### Prepare what's needed for loop ######### 
# Set path
inpath = '~/output/pls/' # Path to directory

#Import relevant data
mapping = fread("~/data/raw/mapping_glasser.txt")
glasser_areas = fread("~/data/raw/glasser_areas.txt") # read in region number for plotting

# Set colours
fake_scico <- scico(3, palette = "vik") #make colour palette

######### 1) Brain plot with PLS weights ######### 

#### Define modalitiy
modality_list = c("SA","CT","ICVF")

####Set up lists to store output

  pls_weights_all = NULL 
  
  ### Ímport disorder data ###
  disorder_df = fread("~/data/processed/h_magma/disorders_traits/FB_SCZ.genes.out")
  disorder_df = disorder_df[,c("GENE","ZSTAT")]
  names(disorder_df)[2] = disorder
  
  for (modality in modality_list) {
    ### Import PLS weights for component 1 and 2 ###
    pls_weights_component_1 = fread(paste(inpath,"SCZ_",modality,"_PLS_Weights_Component_1_ascending.csv", sep = ""))
    pls_weights_component_1$component = "Component 1"
    pls_weights=pls_weights_component_1
    pls_weights$phenotype = modality
    
    ### Collect output for brainplot across all modalities ###
    pls_weights_all = rbind(pls_weights_all,pls_weights)
    
  } 
  
  # Combine with info needed for plotting
  pls_weights_all = merge(pls_weights_all,mapping, by.x="name", by.y = "S_no")
  
  # Prepare data for brain plot
  pls_weights_all$label <- paste0('rh_R_',gsub("_ROI","",pls_weights_all$region)) #adapt label names
  pls_weights_all$boot.weight = ifelse(pls_weights_all$fdr.pvalue<0.05,pls_weights_all$boot.weight,NA)
  pls_weights_all$phenotype = as.factor(pls_weights_all$phenotype)
  pls_weights_all$phenotype = factor(pls_weights_all$phenotype, levels =c("SA","CT","ICVF"))
  
  # and plot
  brainblot = pls_weights_all %>% 
    group_by(phenotype)%>% 
    ggplot() + 
    geom_brain(atlas = glasser, hemi="right",
                          aes(fill = boot.weight), colour="black")+
    scale_fill_viridis_c(na.value = "transparent", name="PLS weights")+
    facet_wrap(c("phenotype"),ncol = 1)+
    theme_void()+ 
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),legend.position = "bottom", legend.title =element_blank())

  ggsave(brainblot,filename = "~/output/plots/subset/SCZ_PLS_weights_brainplot.pdf",height = 6, width = 4)

######### 2)  Get number of significant regions per modality ######### 
  n_sig =  pls_weights_all %>%
    na.omit() %>%
    group_by(phenotype) %>%
    dplyr::summarise(number_sig = n()) %>%
    arrange(desc(number_sig))
  write.xlsx2(n_sig,file = "~/output/pls/SCZ_nsig_pls_weights_subset.xlsx")

######### 3)  Get mean PLS weight regions per modality ######### 
  mean_sig =  pls_weights_all %>%
    na.omit() %>%
    group_by(phenotype) %>%
    dplyr::summarise(mean = mean(boot.weight)) %>%
    arrange(desc(mean))
  write.xlsx2(  mean_sig,file ="~/output/pls/SCZ_mean_pls_weights_subset.xlsx")
  

######### 4) Correlation between PLS weights 
  cor_list<- pls_weights_all %>%
    select(c(phenotype, boot.weight, name)) %>% 
    #select(where((is.factor())~recode(., ICVF = "NDI")))%>% 
    pivot_wider(names_from = phenotype, values_from = boot.weight) %>% 
    rename("NDI"="ICVF")%>%
    select(-name) %>% 
    as.matrix() %>%
    rcorr(type = "spearman")


