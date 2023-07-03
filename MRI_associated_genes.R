#### The genetic relaionships between brain structure and schizophrenia
#### Stauffer et al., 2023

# This script identifies significant genes for each MRI metric and characterise them.
# multiple comparison correction: Matrix decomposition within modalities; additional script

rm(list=ls())
options(stringsAsFactors = FALSE)

library(data.table);library(biomaRt); library(reshape); library(WGCNA); library(ggplot2);
library(gridExtra); library(GenomicRanges); library(dplyr); library(tidyr);library(gtools);library(gplots)
library(corrplot); library(gtools); library(purrr); library(scico);library(ggpubr); library(pander); library(xlsx); 
library("ggVennDiagram");library(ggseg);library(ggsegGlasser);library(ggsegExtra); 

# Import some base data files
Protein_coding = fread("~/data/raw/protein_coding_IDs.txt") # protein coding genes
n_independent_tests_modality= fread("~/output/limatrixdecompmodality.txt") # number of independent tests for significance threshold
n_independent_tests_modality$modality = as.character(n_independent_tests_modality$modality)
mapping = fread("~/data/raw/mapping_glasser.txt") # mapping of brain regions based on HCP


###### Step 1: Restrict genetic data to protein-coding genes, exclude mhc region and combine gene values from 180 regions per MRI metric into one

  # Based on adult brain HI-c data:
  modality_list = c("SA", "CT","ICVF")
  
  p_thresholds_all = NULL
  for (m in modality_list) {
    combined = NULL
    for (i in 1:180){
      adult = fread(paste0("~/data/processed/h_magma/AB/Out_AB_",m,"_plinkmeta",i,".genes.out"))
      adult = adult[adult$GENE %in% Protein_coding$Gene_stable_ID,] # exclude non protein coding genes
      mhcrange = GRanges(seqnames=6, IRanges(25000000,35000000)) # MHC region
      #selecting out genes that are located in MHC region
      ranges = GRanges(seqnames=adult$CHR, IRanges(adult$START, adult$STOP), gene=adult$GENE)
      olap = findOverlaps(ranges, mhcrange)
      mhcgene = ranges[queryHits(olap)]$gene
      
      adult = adult[!(adult$GENE %in% mhcgene),] # remove MHC genes
      adult$ROI = paste0("ROI_",i)
      combined = rbind(combined, adult) #combine output from all ROIs for each modality
    }
   
    # Modality specific threshold
    n_tests_modality = as.vector(n_independent_tests_modality[n_independent_tests_modality$modality==m,1])
    combined$sig_modality = ifelse(combined$P < 0.05/(nrow(adult)*n_tests_modality[[1]]),1,0)
    write.table(combined, file = paste0("~/data/processed/h_magma/AB_combined/AB_",m,"_combined.txt"), row.names = F, col.names = T, quote = F)
    
    # Save thresholds
    modspecific = 0.05/(nrow(adult)*n_tests_modality[[1]])
    p_thresholds$modality = m
    p_thresholds_all = rbind(p_thresholds_all, p_thresholds)
  }
  p_thresholds_all$modspecific= signif(p_thresholds_all$modspecific, digits = 3)
  write.table(p_thresholds_all, file = paste0("~/output/AB_modalityspecific_gene_thresholds.txt"), row.names = F, col.names = T, quote = F)
  
  
  # Based on fetal brain HI-c data:
  for (m in modality_list) {
    combined = NULL
    for (i in 1:180){
      fetal = fread(paste0("~/data/processed/h_magma/FB/Out_FB_",m,"_plinkmeta",i,".genes.out"))
      fetal = fetal[fetal$GENE %in% Protein_coding$Gene_stable_ID,] # exclude non protein coding genes
      mhcrange = GRanges(seqnames=6, IRanges(25000000,35000000)) # MHC region
      #selecting out genes that are located in MHC region
      ranges = GRanges(seqnames=fetal$CHR, IRanges(fetal$START, fetal$STOP), gene=fetal$GENE)
      olap = findOverlaps(ranges, mhcrange)
      mhcgene = ranges[queryHits(olap)]$gene
      
      fetal = fetal[!(fetal$GENE %in% mhcgene),] # remove MHC genes
      fetal$ROI = paste0("ROI_",i)
      combined = rbind(combined, fetal) 
    }
    # Modality specific threshold
    n_tests_modality = as.vector(n_independent_tests_modality[n_independent_tests_modality$modality==m,1])
    combined$sig_modality = ifelse(combined$P < 0.05/(nrow(fetal)*n_tests_modality[[1]]),1,0)
    write.table(combined, file = paste0("~/data/processed/h_magma/FB_combined/FB_",m,"_combined.txt"), row.names = F, col.names = T, quote = F)
    
    # Save thresholds
    modspecific = 0.05/(nrow(fetal)*n_tests_modality[[1]])
    p_thresholds$modality = m
    p_thresholds_all = rbind(p_thresholds_all, p_thresholds)
  }
  p_thresholds_all$modspecific= signif(p_thresholds_all$modspecific, digits = 3)
  write.table(p_thresholds_all, file = paste0("~/output/FB_modalityspecific_gene_thresholds.txt"), row.names = F, col.names = T, quote = F)


###### Step 2: Identify significant genes per MRI metric and 
  # (i) investigate developmental expression trajectories 
  # (ii) GO enrichment 
  # (iii) Plot number of significant genes per brain region
  
  # Load developmental expression data
  ensg_genename = fread("~/~/data/raw/ensg_genename.txt")
  ID = fread("~/data/raw/Psychencode_brainIDfiles.txt")
  data1 = fread("~/data/raw/psychencode_scaledlogtransformed_genexpr.txt")
  
  ID$period = ifelse(ID$Window < 5, 0, NA) 
  ID$period = ifelse(ID$Window > 5, 1, ID$period)
  
  # names for developmental windows
  windownames<- c("8-9","12-13","16-17","19-22","35pcw \n 4mos","0.5-2.5","3-11","13-19","21-40")
  
  modality_list = c("SA", "CT","ICVF")
  
  # Extract number of significant genes, developmental expression profiles, run GO
    AB_sig_2 = NULL
    FB_sig_2 = NULL
    union_sig_2 = NULL
    Means_with_ID2 = NULL
    development = NULL
    union_allphenos = NULL
    gostres_df_all = NULL
    n_sig_per_region_all=NULL
    
    for (modality in modality_list){
      print(modality)
      ## Read in data
      fetal = fread(paste0("~/data/processed/h_magma/FB_combined/FB_",modality,"_combined.txt"))
      adult = fread(paste0("~/data/processed/h_magma/AB_combined/AB_",modality,"_combined.txt"))
      union = rbind(fetal,adult)
      union$Phenotype = modality
      union_allphenos = rbind(union_allphenos, union)
      
      ## Restrict to significant genes
      FB_sig = subset(fetal, fetal[[sig_modality]]  == 1)
      AB_sig  = subset(adult, adult[[sig_modality]] == 1)

      ## Extract unique significant genes
      AB_sig = unique((AB_sig)[,c("GENE","START", "STOP")]) 
      AB_sig$name = paste0("Phenotype_", modality)
      AB_sig_2 = rbind(AB_sig_2, AB_sig) 
      
      FB_sig = unique((FB_sig)[,c("GENE","START", "STOP")]) 
      FB_sig$name = paste0("Phenotype_", modality)
      FB_sig_2 = rbind(FB_sig_2, FB_sig) 
      
      union_sig = rbind(AB_sig, FB_sig)
      union_sig = unique((union_sig)[,c("GENE","START", "STOP")])
      union_sig$name = paste0("Phenotype_",modality)
      union_sig_2 = rbind(union_sig_2, union_sig)
      
      ## developmental trajectories in union gene list
      sig_in_data1 = data1[data1$GENE %in% union_sig$GENE,] 
      
      Means = colMeans(sig_in_data1[,2:608])
      Names = names(sig_in_data1[,2:608])
      
      Means = as.data.frame(cbind(Names, Means))
      Means_with_ID = merge(ID, Means, by.x = "ID", by.y = "Names")
      Means_with_ID$name = paste0(modality)
      Means_with_ID2 = rbind(Means_with_ID2, Means_with_ID)
      
      Cortex_only = subset(Means_with_ID, region_superbroad == "Neocortex")

      # Run linear regression analyses to investigate diff between adult and fetal brains
      s = summary(lm(Means ~ period , data = Cortex_only,na.action = na.omit))
      coef  = coefficients(s)[2,]
      coef$name = modality
      coef = as.data.frame(coef)
      development = rbind(development, coef)
      
      #Run GO enrichment on significant genes for each modality (fetal hi-c)
      background = as.data.frame(unique(combined$GENE))
      gostres = gprofiler2::gost(query = FB_sig$`unique(FB_sig$GENE)`, organism = "hsapiens",
                                 ordered_query = F,
                                 correction_method = c("g_SCS"),
                                 custom_bg = background$`unique(combined$GENE)`,
                                 sources = c("GO:BP"))
      
      gostres_df= gostres$result
      gostres_df$modality = modality
      gostres_df_all = rbind(gostres_df_all,gostres_df)
      
      ## Extract number of significant genes per region 
      n_sig_genes_per_region = FB_sig$`unique(FB_sig$GENE)` %>%
        group_by(ROI) %>%
        dplyr::summarise(number_sig_genes = n()) %>%
        arrange(desc(number_sig_genes))
      
      n_sig_genes_per_region = as.data.frame(merge(n_sig_genes_per_region,mapping[,1:2], by.x = "ROI",by.y="S_no",all=TRUE))
      n_sig_genes_per_region$modality = modality
      n_sig_per_region_all = rbind(n_sig_per_region_all,n_sig_genes_per_region)
    }
    
    ## Store results / Plot results
    
    #Export significant genes per modality
    write.table(AB_sig_2,file = "~/output/AB_significant_genes_sig_modality.txt",row.names = F, col.names = T, quote = F)
    write.table(FB_sig_2,file = "~/output/FB_significant_genes_sig_modality.txt",row.names = F, col.names = T, quote = F)
    write.table(union_sig_2,file = "~/output/union_significant_genes_sig_modality.txt",row.names = F, col.names = T, quote = F)
    
    # Plot developmental trajectories neocortex
    Cortex_only_allregions = subset(Means_with_ID2, region_superbroad == "Neocortex")
    p_dev_trajectory_cortex = ggplot(Cortex_only_allregions,aes(x=Window, y= as.numeric(as.character(Means)), fill=name, color=name, linetype =regiontype)) + 
      xlim(1,10) + ylab("Normalized expression") + geom_smooth(method = "loess", formula = y ~ x, se = TRUE) + 
      xlab("Developmental Stages") + scale_x_continuous(breaks=seq(1,14,1)) + theme_classic() + geom_vline(xintercept = 5) +
      theme_classic() + geom_vline(xintercept = 5) + 
      scale_fill_discrete(name = "Phenotype")+scale_color_discrete(name = "Phenotype")+
      labs(linetype = "Phenotype group")+
      theme(axis.text.y = element_text(size = 15),axis.text.x = element_text(size = 15),
            axis.title = element_text(size = 15))+
      scale_x_continuous(breaks=seq(1,9,1),labels= windownames)
    ggsave(p_dev_trajectory_cortex,file="~/output/dev_trajectory_cortex_sig_modality.pdf",height = 5,width=9)
    
    ## store developmental lme output
    development2 = development
    development2$pfdr = p.adjust(development2$Pr...t.., method = "fdr")
    development2 = development2[,c(5,1:4,6)]
    development2$sig = add.significance.stars(development2$pfdr, cutoffs = c(0.05, 0.01, 0.001))
    development2$pfdr = format(development2$pfdr, digits = 3)
    development2$Pr...t.. = format(development2$Pr...t.., digits = 3)
    write.table(development2, "~/output/development_sig_modality.txt", row.names = F, col.names = T, quote = F)

    # Store GO output
    gostres_df_all = gostres_df_all[gostres_df_all$p_value<0.05,]
    write.xlsx2(gostres_df_all, "~/output/gostres_df_all_modality_sig_genes.xlsx")
  
    # Plot number of significant genes per region onto brain 
    n_sig_per_region_all$label <- paste0('rh_R_',n_sig_per_region_all$region)
    brainplot = brainplot = n_sig_per_region_all%>% 
      mutate(across(modality, factor, levels=c("SA", "CT", "ICVF"))) %>%
      group_by(modality)%>% 
      ggplot() +
      geom_brain(atlas = glasser, hemi="right",aes(fill = number_sig_genes), colour="black")+ 
      scale_fill_continuous(trans = 'reverse',na.value = "transparent")+
      theme_void()+
      theme(legend.position = "none")+
      facet_wrap("modality")
    ggsave(brainplot,filename = "~/output/nsig_genes_brainplot.pdf",height = 10, width = 12)

    
###### Step 3: Investigate genetic overlap between MRI metrics
    FB_sig_2_ven = FB_sig_2[,c(1,4)]
    FB_sig_2_ven <- split(FB_sig_2_ven, f = FB_sig_2_ven$name )
    FB_sig_2_ven <- lapply(FB_sig_2_ven, subset, select = -name)
    FB_sig_2_ven <- lapply(FB_sig_2_ven, FUN = function(x) x[[1]])
    
    # Plot number of overlapping and unique genes
    cols = c("CT" = "#F8766D", "SA" = "#00BFC4", "NDI" = "#7CAE00")
    p_ven = ggVennDiagram(x =FB_sig_2_ven,category.names = c("CT","NDI", "SA"), label_alpha = 0)+
      scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
      scale_colour_manual(values =c("#F8766D","#7CAE00","#00BFC4"))
    ggsave(p_ven,filename = "~/output/venn_between_phenos_sig_modality.png",height = 3.5, width = 3.5)
    
    # Extract overlapping genes between all phenos
    shared_genes =  process_region_data(Venn(FB_sig_2_ven))
    shared_genes2 = shared_genes
    shared_genes2$item = toString(shared_genes2$item )
    shared_genes = as.data.frame(unlist(shared_genes[[7,3]]))
    shared_genes = merge(shared_genes, ensg_genename, by.x = "unlist(shared_genes[[7, 3]])", by.y ="Gene_stable_ID")
    names(shared_genes)[1] = "GENE"
    write.table(shared_genes,file ="~/output/shared_genes_between_phenos_sig_modality.txt",row.names = F, col.names = T, quote = F)
    