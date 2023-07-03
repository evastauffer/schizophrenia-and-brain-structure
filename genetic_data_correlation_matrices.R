#### The genetic relaionships between brain structure and schizophrenia
#### Stauffer et al., 2023

## This script calculates correlation matrices for genetic data (adult and fetal brain)

# Prepare dataframes
  # restrict to protein coding genes, remove mhc region
  # for each region extract genes and z scores
  # combine these results within each modality
# Construct correlation matrix

rm(list=ls())
options(stringsAsFactors = FALSE)

library(data.table); library(biomaRt); library(reshape); library(ggplot2);
library(gridExtra); library(GenomicRanges); library(dplyr); library(tidyr);library(gtools);library(gplots);
library(corrplot); library(gtools); library(purrr); library(scico);library(pander); library(sna); library(psych);
library(heatmaply);library(easyGgplot2);library(ComplexHeatmap);library(egg);library(ggpubr);library(circlize);library(reshape2);
library(tidyverse);library(grid);library(plyr)


### Import relevant data 
region_number = fread("/region_number.csv") #names and region numbers
region_number$S_no = paste0("ROI_",region_number$S_no)
Protein_coding = fread("protein_coding_IDs.txt") # protein coding genes
GWAS_ID = fread("GWAS_ID.txt") # SUbjects that went into GWAS


###-----------------First for Adult brain genetic data, not controlled for global values----------------###

# List of phenotypes
modality_list = c("ICVF","CT","SA")

for (modality in modality_list) {

  datalist = list()
  ### First extract association values for each gene and each region ###
   for (i in 1:180){
    adult = fread(paste0("Out_AB_",modality,"_plinkmeta",i,".genes.out"))
    adult = adult[adult$GENE %in% Protein_coding$Gene_stable_ID,] # exclude non protein coding genes
    mhcrange = GRanges(seqnames=6, IRanges(25000000,35000000)) # MHC region
    ranges = GRanges(seqnames=adult$CHR, IRanges(adult$START, adult$STOP), gene=adult$GENE) #selecting out genes that are located in MHC region
    olap = findOverlaps(ranges, mhcrange)
    mhcgene = ranges[queryHits(olap)]$gene
  
    adult = adult[!(adult$GENE %in% mhcgene),] # remove MHC genes
  
    gene_zstats = adult[,c("GENE","ZSTAT")] # Extract gene list 
    names(gene_zstats)[names(gene_zstats) == 'ZSTAT'] <- paste0("ROI_", i)
    datalist[[i]] = gene_zstats
  
  }
  # Combine results from all 180 regions
  data_adult = join_all(datalist, by='GENE', type='left')
  data_adult = data.frame(column_to_rownames(data_adult, var = "GENE"))
  
  write.table(data_adult, file = paste0("AB_",modality,"_combined_regions.txt"), row.names = F, col.names = T, quote = F)

  ### Calculate correlation matrix ###
  data_adult = scale(data_adult, center = T, scale = T) # scale data
  cormat_spear_genes <- cor(as.matrix(data_adult),method = 'spearman', use = "pairwise.complete.obs") # spearman correlation
  write.table(cormat_spear_genes, file = paste0("AB_",modality,"_cormat_spear.txt"), row.names = F, col.names = T, quote = F)
}



###-----------------Now for fetal brain genetic data, not controlled for global values----------------###

modality_list = c("ICVF","CT","SA")
for (modality in modality_list) {
  print(modality)
   datalist = list()
  ### First extract association values for each gene and each region ###
    for (i in 1:180){
     fetal = fread(paste0("/Out_FB_",modality,"_plinkmeta",i,".genes.out"))
     fetal = fetal[fetal$GENE %in% Protein_coding$Gene_stable_ID,] # exclude non protein coding genes
     mhcrange = GRanges(seqnames=6, IRanges(25000000,35000000)) # MHC region
     #selecting out genes that are located in MHC region
     ranges = GRanges(seqnames=fetal$CHR, IRanges(fetal$START, fetal$STOP), gene=fetal$GENE)
     olap = findOverlaps(ranges, mhcrange)
     mhcgene = ranges[queryHits(olap)]$gene
  
     fetal = fetal[!(fetal$GENE %in% mhcgene),] # remove MHC genes
  
     gene_zstats = fetal[,c("GENE","ZSTAT")] # Extract gene list 
     names(gene_zstats)[names(gene_zstats) == 'ZSTAT'] <- paste0("ROI_", i)
     datalist[[i]] = gene_zstats
  
   }
  # Combine results from all 180 regions
  data_fetal = join_all(datalist, by='GENE', type='left')
  data_fetal = data.frame(column_to_rownames(data_fetal, var = "GENE"))
  write.table(data_fetal, file = paste0("FB_",modality,"_combined_regions.txt"), row.names = T, col.names = T, quote = F)

  
  ### Calculate correlation matrix ###
  data_fetal = scale(data_fetal, center = T, scale = T) # scale data
  cormat_spear_genes <- cor(as.matrix(data_fetal),method = 'spearman', use = "pairwise.complete.obs") # spearman correlation
  write.table(cormat_spear_genes, file = paste0("FB_",modality,"_cormat_spear.txt"), row.names = F, col.names = T, quote = F)
}




