#### The genetic relaionships between brain structure and schizophrenia
#### Stauffer et al., 2023

# This script calculates cophonetic correlations within each modality between GS and SC

rm(list=ls())
options(stringsAsFactors = FALSE)

library(data.table);library(ComplexHeatmap);library(scico);library(circlize);library(data.table);
library(gtools);library(ggpubr);library(dendextend);library(xlsx);library(dplyr);

### Import relevant data 
region_number = fread("~/data/raw/region_number.csv") #names and region numbers
region_number$S_no = paste0("ROI_",region_number$S_no)

GWAS_ID = fread("/mnt/b2/home4/arc/ems206/cortex_project_04_03_2021/data/raw/GWAS_ID.txt") # SUbjects that went into GWAS


### Set up loop
out = NULL
plots = NULL

modality_list = c("SA","CT","ICVF")

for (modality in modality_list) {
  ### Import and scale datafiles with genetic effects for each region
  genetic_data = fread(paste0("~/data/processed/h_magma/FB_combined/FB_",modality,"_combined_regions.txt"))
  genetic_data$GENE = NULL
  genetic_data = scale(genetic_data, center = T, scale = T) # scale data
  
  # Create dendrograms
  genetic_data_t = t(genetic_data) 
  dend_gene <- genetic_data_t %>% 
    dist(method = "euclidian") %>% # calculate a distance matrix based on euclidian distance, 
    hclust(method = "ward.D2") %>% # on it compute hierarchical clustering using the "ward" method, 
    as.dendrogram # and lastly, turn that object into a dendrogram.

  ### Repeat for imaging data
  imaging_data = fread(paste0("~/data/processed/imaging_data/",modality,"_regional.txt")) # full UKB imaging set
  imaging_data = imaging_data[imaging_data$FID %in% GWAS_ID$IID,] # restrict to subjects that went into GWAS
  imaging_data = imaging_data[,3:182]
  setnames(imaging_data, old = c(names(imaging_data)), new = c(region_number$S_no))
  imaging_data = scale(imaging_data, center = T, scale = T) # scale data
  
  # Create dendrograms
  imaging_data_t = t(imaging_data) dend_imaging <- imaging_data_t %>% # take the a vector from 1 to 5
    dist(method = "euclidian") %>% # calculate a distance matrix based on euclidian distance, 
    hclust(method = "ward.D") %>% # on it compute hierarchical clustering using the "ward" method, 
    as.dendrogram # and lastly, turn that object into a dendrogram.
  
  ### Calculate cophenteic correlation, permutation test
  
  ## "The value can range between -1 to 1. With near 0 values meaning that the two trees are not statistically similar. For exact p-value one should result to a permutation test. One such option will be to permute over the labels of one tree many times, and calculating the distribution under the null hypothesis (keeping the trees topologies constant)." 
  ## Null hypothesis is correlation = 0 
  ## two tailed:  H1  cor != 0
  ## one tailed: H1 cor > 0 OR H1 cor < 0
  
  cor_coph <- cor_cophenetic(dend_imaging,dend_gene, method = "spearman")
  # Permutation
  set.seed(2248)
  R <- 1000
  cor_coph_results <- numeric(R)
  dend_mixed <- dend_gene
  for (i in 1:R){
    print(i)
    dend_mixed <- sample.dendrogram(dend_mixed, replace = F) #just swapping labels for a permutation p-value
    cor_coph_results[i] <- cor_cophenetic(dend_imaging, dend_mixed)
  }
  cor_coph_results <- cor_coph_results %>% as.data.frame() 
  mean = mean(cor_coph_results$.)
  sd = sd(cor_coph_results$.)
  z = (cor_coph-mean)/sd
  p = 2*pnorm(q=z, lower.tail=FALSE)
  
  p = as.data.frame(p)
  p$phenotype = paste(modality)
  p$corr = cor_coph
  p$sum = sum(abs(cor_coph_results)>=cor_coph)
  out = rbind(out,p)
  
  plots[[modality]] = ggplot(cor_coph_results, aes(x=.)) +
    geom_density(fill = "#59c9a5", colour = "#59c9a5", alpha = 0.8)  +
    geom_vline(xintercept=1, size = 1.2, linetype = "dashed", colour="#2e1f27") +
    geom_vline(xintercept=0, size = 1.2, linetype = "dashed", colour="#2e1f27") +
    geom_vline(xintercept=cor_coph, size=1.2, colour="blue") +
    ylab("Density") +
    xlab("Cophenetic correlation coefficent") + 
    ggtitle(paste0(modality))+
    theme_bw()
  
}

  ### Save plots
  p_cophonetic = ggpubr::ggarrange(
    ggpubr::ggarrange(plotlist =plots, ncol =3 ,nrow =1,common.legend = F))
  ggsave(p_cophonetic,filename = '~/output/plots/FB_cophonetic_permutation.png',height = 3, width = 12)

  ### Store cophonetic correlations for all phenos
  out2 = out[,c(2,3,1)]
  out2$pfdr = p.adjust(out2$p, method = "fdr")
  out2$corr = round(out2$corr,2)
  out2 = out2 %>% mutate_at(c(3,4),.funs = format, scientific=TRUE, digits = 2)
  write.table(out2,  file="~output/FB_cophonetic_correlations_permuted_subset.txt",row.names = F, col.names = T, quote = F)
  write.xlsx(out2, "~/output/FB_cophonetic_correlations_permuted_subset.xlsx")

