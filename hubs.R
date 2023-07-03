#### The genetic relaionships between brain structure and schizophrenia
#### Stauffer et al., 2023


# 1) calculated hubness for imaging and genetic data (hubness is the rowsum for each node)
# 2) Plots hubness onto brain
# 3) Correlat hubness based on genetic data with hubness based on imaging data
# 4) Compare Hubness to PLS weights

library(data.table);library(dplyr);library(ggplot2);
library(ggpubr);library(igraph);library(tibble);library(brainGraph);
library(ggseg); library(ggsegGlasser);library(gtools);library(network);
library(ggnet); library(broom);library(cmocean);library(xlsx)

inpath = '/mnt/b2/home4/arc/ems206/cortex_project_04_03_2021/output/pls/'

# Import info on regions for plotting
mapping = fread("/mnt/b2/home4/arc/ems206/cortex_project_04_03_2021/data/raw/mapping_glasser.txt")
mapping = mapping[mixedorder(mapping$S_no),]


# Prepare some lists to store output
data_pls_all=NULL

table_cor_hub_gs = NULL
table_cor_hub_sc = NULL

# Set up loop
modality_list = c("SA","CT","ICVF")
  for (modality in modality_list) {
    
    ##### Hubness for genetic data #####
    #Read in genetic matrix
    cormat_genes = fread(paste0("~1/data/processed/cormat_genes/FB_",modality,"_cormat_spear.txt"))
    cormat_genes_p = cormat_genes                
    cormat_genes_p[cormat_genes_p < 0] = 0    
    cormat_genes = as.matrix(cormat_genes_p)
    
    # Hubs
    df_hs_g = as.data.frame(rowMeans(cormat_genes))
    names(df_hs_g) = c("hs_g")
    df_hs_g$ROI = paste("ROI",1:180, sep = "_")

    ##### Hubness  for imaging #####
    #Read in imaging matrix
    cormat = fread(paste0("~/data/processed/imaging_data/cormat_imaging/structural_covariance_",modality,".txt"))
    cormat_p = cormat                 # Duplicate data frame
    cormat_p[cormat_p < 0] = 0     # Set negative values to 0
    cormat_genes = as.matrix(cormat_p)
    
    # Hubs
    df_hs_i = as.data.frame(rowMeans(cormat))
    names(df_hs_i) = c("hs_i")
    df_hs_i$ROI = paste("ROI",1:180, sep = "_")

    ##### Combine genetic and imaging #####
    df_hs = merge(df_hs_g,df_hs_i, by = "ROI")
    df_hs$modality = modality
    
    ### Import PLS weights for component 1 ###
    pls_weights_component_1 = fread(paste(inpath,"SCZ_",modality,"_PLS_Weights_Component_1_ascending.csv", sep = ""))
    pls_weights_component_1$component = "Component 1"
    pls_weights=pls_weights_component_1
    
    # Combine with hubs and strength
    data_pls = merge(pls_weights, df_hs, by.x=c("name"), by.y=c("ROI"))
    data_pls_all = rbind(data_pls_all, data_pls)
    data_pls2 = data_pls[,c(1:5,7:9)]
    write.table(data_pls2,paste("/~/output/SCZ_",modality,"_pls_weight_hubs.txt",sep = ""),row.names = F,col.names = T,quote = F)
 
    ### Correlation table hubness, sc, gs
    correlation_g = cor.test(data_pls$boot.weight,data_pls$hs_g, method = "spearman", exact = F)
    cor_hub_gs = tidy(correlation_g)
    cor_hub_gs$modality = modality
    table_cor_hub_gs = rbind(table_cor_hub_gs,cor_hub_gs)
    
    correlation_i = cor.test(data_pls$boot.weight,data_pls$hs_i, method = "spearman", exact = F)
    cor_hub_sc = tidy(correlation_i)
    cor_hub_sc$modality = modality
    table_cor_hub_sc = rbind(table_cor_hub_sc,cor_hub_sc)
    
  }
  
  ## Export table
  table_cor_hub_sc$pfdr = p.adjust(table_cor_hub_sc$p.value, method = "fdr")
  table_cor_hub_gs$pfdr = p.adjust(table_cor_hub_gs$p.value, method = "fdr")
  write.xlsx(table_cor_hub_sc, "~/output/SCZ_cor_hub_sc.xlsx")
  write.xlsx(table_cor_hub_gs, "~/output/SCZ_cor_hub_gs.xlsx")

  
  ## Visualise results
  
  # Set levels for modality
  data_pls_all$modality  = factor(data_pls_all$modality,levels = c("SA","CT", "ICVF"))
  
  # Plot hub score genetic X PLS weight
  plot_list_hub_g = ggplot(data_pls_all, aes(x=boot.weight, y=hs_g))+ 
    geom_point(size = 0.8) +
    geom_smooth(method=lm,colour = "red")+ 
    facet_wrap(~ modality, ncol = 1,nrow = 3) +
    theme_bw()+
    stat_cor(method = "spearman",label.x.npc = "left",label.y.npc = "top",cor.coef.name = c("rho"), size =5,aes(label = ..r.label..))+
    theme(axis.text=element_text(size=14),axis.title=element_blank(),plot.title = element_blank(),strip.background = element_blank(),strip.text.x = element_blank())
  
   ggsave(plot_list_hub_g ,filename = "~/output/plots/SCZ_scatterplot_hubs_g_pls.pdf",height = 4.5, width = 3)
  
  
  # Plot hub score imaging X PLS weight
  plot_list_hub_i= ggplot(data_pls_all, aes(x=boot.weight, y=hs_i))+ 
    geom_point(size = 1) +
    geom_smooth(method=lm,colour = "red")+ 
    facet_wrap(~ modality, ncol = 1,nrow = 3) +
    theme_bw()+
    stat_cor(method = "spearman",label.x.npc = "left",label.y.npc = "top",cor.coef.name = c("rho"), size =5,aes(label = ..r.label..))+
    theme(axis.text=element_text(size=14),axis.title=element_blank(),plot.title = element_blank(),strip.background = element_blank(),strip.text.x = element_blank())
  ggsave(plot_list_hub_i ,filename = "~/output/plots/SCZ_scatterplot_hubs_i_pls.pdf",height = 4.5, width = 3)
  

  # Plot hub score genetic versus imaging
  plot_list_hub_scatter=  ggplot(data_pls_all, aes(x=hs_g, y=hs_i))+ 
    geom_point(size = 1) +
    geom_smooth(method=lm,colour = "red")+ 
    facet_wrap(~ modality, ncol = 1,nrow = 3) +
    theme_bw()+
    stat_cor(method = "spearman",label.x.npc = "left",label.y.npc = "top",cor.coef.name = c("rho"), size =5,aes(label = ..r.label..))+
    theme(axis.text=element_text(size=14),axis.title=element_blank(),plot.title = element_blank(),strip.background = element_blank(),strip.text.x = element_blank())
  
  ggsave(plot_list_hub_scatter,filename = "~/output/plots/SCZ_scatterplot_hubs_g_i.pdf",height = 4.5, width = 3)
  
  
  # Plot hubness on brain
  hub_brain = merge(data_pls_all, mapping, by.x = "name", by.y = "S_no")
  hub_brain$label = paste0('rh_R_',hub_brain$region)
  
  hub_brain_sc = hub_brain
  names(hub_brain_sc)[8] = "hub"
  hub_brain_sc$hs_g = NULL
  hub_brain_sc$type_cor = "SC"
  
  hub_brain_gs = hub_brain
  names(hub_brain_gs)[7] = "hub"
  hub_brain_gs$hs_i = NULL
  hub_brain_gs$type_cor = "GS"
  
  hub_brain_plot = rbind(hub_brain_sc, hub_brain_gs)
  hub_brain_plot$type_cor  = factor(hub_brain_plot$type_cor ,levels = c("GS","SC"))
  hub_brain_plot$modality  = factor(hub_brain_plot$modality,levels = c("SA","CT", "ICVF"))
  
  p_brain_hub = hub_brain_plot%>% 
    group_by(type_cor,modality)%>% 
    ggplot() +
    geom_brain(atlas = glasser, hemi="right",aes(fill = hub), colour="black")+ 
    labs(fill = "Hubness \n imaging")+ scale_fill_cmocean(alpha = 1, start = 0, end = 1, direction = 1, discrete = FALSE, name = "haline",na.value="transparent")+
    theme_void()+
    facet_wrap(~ modality + type_cor, ncol = 2,nrow = 3) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),legend.position = "bottom", legend.title =element_blank())
 ggsave(p_brain_hub,filename ="~/output/plots/SCZ_brainplot_hubbness.pdf",height = 6, width = 8)
 
 

