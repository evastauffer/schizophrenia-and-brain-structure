#### The genetic relaionships between brain structure and schizophrenia
#### Stauffer et al., 2023

# This script extract significant genes for schizophrenia, runs GO analysis and devvelopmental expresison profiles; plot with MRI and SCZ

rm(list=ls())
options(stringsAsFactors = FALSE)

library(data.table);library(GenomicRanges);library(dplyr);
library(ggplot2);library(xlsx);library(gprofiler2)

# Read in some data that we need for selecting genes
  Protein_coding = fread("~/data/raw/protein_coding_IDs.txt") # protein coding genes

# Read in data for developmental trajectories
  # Load relevant data 
  ensg_genename = fread("~/data/raw/ensg_genename.txt")
  
  ID = fread("~/Psychencode_brainIDfiles.txt")
  data1 = fread("~/psychencode_scaledlogtransformed_genexpr.txt")
  
  ID$period = ifelse(ID$Window < 5, 0, NA) 
  ID$period = ifelse(ID$Window > 5, 1, ID$period)
  
  # names for developmental windows
  windownames<- c("8-9","12-13","16-17","19-22","35pcw \n 4mos","0.5-2.5","3-11","13-19","21-40")


  #Ímport data fetal
  disorder_df = fread("~/data/processed/h_magma/disorders_traits/FB_SCZ.genes.out")
  disorder_df = disorder_df[disorder_df$GENE %in% Protein_coding$Gene_stable_ID,] # exclude non protein coding genes
  mhcrange = GRanges(seqnames=6, IRanges(25000000,35000000)) # MHC region
  #selecting out genes that are located in MHC region
  ranges = GRanges(seqnames=disorder_df$CHR, IRanges(disorder_df$START, disorder_df$STOP), gene=disorder_df$GENE)
  olap = findOverlaps(ranges, mhcrange)
  mhcgene = ranges[queryHits(olap)]$gene
  
  disorder_df = disorder_df[!(disorder_df$GENE %in% mhcgene),] # remove MHC genes
  disorder_df$pbon = 0.05/nrow(disorder_df)
  disorder_df_sig_fetal = disorder_df[disorder_df$P<disorder_df$pbon,]
  write.table(disorder_df_sig_fetal, "~/output/genes_sig_scz_fetal.txt", quote = F, row.names = F)
  
  ## Run GO enrichment on significant fetal h_magma
  background = as.data.frame(unique(disorder_df$GENE))
  gostres = gprofiler2::gost(query = disorder_df_sig_fetal$GENE, organism = "hsapiens",
                             ordered_query = F,
                             correction_method = c("g_SCS"),
                             custom_bg = disorder_df$GENE,
                             sources = c("GO:BP"))
  
  gostres_df= gostres$result

  write.xlsx2(gostres_df, "/mnt/b2/home4/arc/ems206/cortex_project_04_03_2021/output/gostres_df_SCZ.xlsx")
  
  ## Add in adult h-magma for developmental trajectories
  disorder_df = fread("~/data/processed/h_magma/AB_SCZ.genes.out")
  disorder_df = disorder_df[disorder_df$GENE %in% Protein_coding$Gene_stable_ID,] # exclude non protein coding genes
  mhcrange = GRanges(seqnames=6, IRanges(25000000,35000000)) # MHC region
  #selecting out genes that are located in MHC region
  ranges = GRanges(seqnames=disorder_df$CHR, IRanges(disorder_df$START, disorder_df$STOP), gene=disorder_df$GENE)
  olap = findOverlaps(ranges, mhcrange)
  mhcgene = ranges[queryHits(olap)]$gene
  
  disorder_df = disorder_df[!(disorder_df$GENE %in% mhcgene),] # remove MHC genes
  disorder_df$pbon = 0.05/nrow(disorder_df)
  disorder_df_sig_adult = disorder_df[disorder_df$P<disorder_df$pbon,]
  
  union_sig = rbind(disorder_df_sig_fetal, disorder_df_sig_adult)
  union_sig = unique((union_sig)[,c("GENE","START", "STOP")])
  
  ## developmental trajectories in union gene list
  sig_in_data1 = data1[data1$GENE %in% union_sig$GENE,] 
  
  Means = colMeans(sig_in_data1[,2:608])
  Names = names(sig_in_data1[,2:608])
  
  Means = as.data.frame(cbind(Names, Means))
  Means_with_ID = merge(ID, Means, by.x = "ID", by.y = "Names")

  Cortex_only = subset(Means_with_ID, region_superbroad == "Neocortex")
  
  # Run linear regression analyses to investigate diff between adult and fetal brains
  s = summary(lm(Means ~ period , data = Cortex_only,na.action = na.omit))
  coef  = coefficients(s)[2,]
  coef = as.data.frame(coef)
  write.xlsx(coef,"~/output/development_SCZ.xlsx")
  
  Cortex_only$name = "SCZ"
  
  # plot expression trajectories with trajectories of MRI metrics
  Cortex_only_MRI = fread("~/output/Cortex_only_allregions_sig_modality.txt")
  Cortex_only_MRI$regiontype =NULL
  Cortex_only_MRI$V12=NULL
  Cortex_only_MRI$V13 =NULL
  Cortex_only_MRI$name = ifelse(Cortex_only_MRI$name == "ICVF", "NDI", Cortex_only_MRI$name)
  
  cortex_only_mri_scz = rbind(Cortex_only_MRI,Cortex_only)
  p_dev_trajectory_cortex = ggplot(cortex_only_mri_scz,aes(x=Window, y= as.numeric(as.character(Means)), fill=name, color=name)) + 
    xlim(1,10) + ylab("Normalized expression") + geom_smooth(method = "loess", formula = y ~ x, se = TRUE) + 
    xlab("Developmental Stages") + scale_x_continuous(breaks=seq(1,14,1)) + theme_classic() + geom_vline(xintercept = 5) +
    theme_classic() + geom_vline(xintercept = 5) + 
    scale_fill_discrete(name = "Phenotype")+scale_color_discrete(name = "Phenotype")+
    labs(linetype = "Phenotype group")+
    theme(axis.text.y = element_text(size = 15),axis.text.x = element_text(size = 15),
          axis.title = element_text(size = 15))+
    scale_x_continuous(breaks=seq(1,9,1),labels= windownames)
  ggsave(p_dev_trajectory_cortex,file="~/output/plots/subset/dev_trajectory_cortex_SCZ.pdf",height = 5,width=9)
  

