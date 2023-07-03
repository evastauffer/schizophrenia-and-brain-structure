#### The genetic relaionships between brain structure and schizophrenia
#### Stauffer et al., 2023

# Extract overlapping genes between MRI modality and SCZ, test whether overlap is significant

library(data.table);library(GenomicRanges);library(dplyr);library(ggVennDiagram);
library(ggplot2);library(xlsx)

### Import gene names
ensg_genename = fread("~/data/raw/ensg_genename.txt")

### Ímport significant genes for scz
  disorder_df_sig = fread("~/output/genes_sig_scz_fetal.txt")

### Import significant genes per modality, get overlapping genes
  modality_genes = fread("~/output/FB_significant_genes_sig_modality.txt")
  modality_genes$name = gsub("Phenotype_","",modality_genes$name)
  
  modality_list = c("CT","SA", "ICVF")
  sig_genes_modality_disorder = NULL
  
  for (m in modality_list) {
    temp_genes = modality_genes[modality_genes$name==m,]
    temp_genes = merge(temp_genes,disorder_df_sig, by.x = "GENE",by.y ="GENE")
    sig_genes_modality_disorder = rbind(sig_genes_modality_disorder,temp_genes)
  }
  
  sig_genes_modality_disorder = merge(sig_genes_modality_disorder,ensg_genename,by.x="GENE",by.y="Gene_stable_ID")
  sig_genes_modality_disorder = sig_genes_modality_disorder[order(sig_genes_modality_disorder$name,sig_genes_modality_disorder$CHR),]                
  sig_genes_modality_disorder_sup = sig_genes_modality_disorder[,c(4,5,1,14,2,3)]
  write.xlsx(sig_genes_modality_disorder_sup, "~/output/genes_modality_scz_overlap_subset.xlsx")

### Extract number of shared significant genes per phenotype
  n_genes = sig_genes_modality_disorder %>%
    group_by(name) %>%
    dplyr::summarise(number_sig_genes = n()) %>%
    arrange(desc(number_sig_genes))

### Extract overlapping genes between all phenos and scz
  shared_genes =  process_region_data(Venn(sig_genes_modality_disorder_ven))
  shared_genes2 = shared_genes
  shared_genes2$item = toString(shared_genes2$item )
  shared_genes = as.data.frame(unlist(shared_genes[[15,3]]))

### Test whether overlap is higher than expected by chance by resampling 
  Protein_coding = fread("~/data/raw/protein_coding_IDs.txt") # protein coding genes

  # Import h-magma data on scz
  disorder_df = fread("~/data/processed/h_magma/disorders_traits/FB_SCZ.genes.out")
  disorder_df = disorder_df[disorder_df$GENE %in% Protein_coding$Gene_stable_ID,] # exclude non protein coding genes
  mhcrange = GRanges(seqnames=6, IRanges(25000000,35000000)) # MHC region
  #selecting out genes that are located in MHC region
  ranges = GRanges(seqnames=disorder_df$CHR, IRanges(disorder_df$START, disorder_df$STOP), gene=disorder_df$GENE)
  olap = findOverlaps(ranges, mhcrange)
  mhcgene = ranges[queryHits(olap)]$gene
  disorder_df = disorder_df[!(disorder_df$GENE %in% mhcgene),] # remove MHC genes

  # set up loop
  modality_list = c("CT","SA", "ICVF")
  df_odds_all = NULL
  
  for (m in modality_list) {
    temp_genes = modality_genes[modality_genes$name==m,]
  
    # Hypergeometric test on real data
    go.obj <- newGeneOverlap(listA = disorder_df_sig$GENE, 
                             listB = temp_genes$GENE, 
                             genome.size = length(disorder_df$GENE))
    go.obj <- testGeneOverlap(go.obj)
    odds = go.obj@odds.ratio
    
    # Hypergeometric test permuted, randomly select genes for MRI
    odds_rand_list = NULL
    n_genes = nrow(temp_genes)
    for (i in 1:10000) {
      temp_genes_perm = sample_n(disorder_df,n_genes, replace = T)
      go.obj.rand <- newGeneOverlap(listA = disorder_df_sig$GENE, 
                                    listB = temp_genes_perm$GENE, 
                                    genome.size = length(disorder_df$GENE))
      go.obj.rand <- testGeneOverlap(go.obj.rand)
      odds_rand = go.obj.rand@odds.ratio
      odds_rand_list = rbind(odds_rand_list,odds_rand)
    }
    
    df_odds = as.data.frame(cbind(odds_rand_list, odds))
    df_odds$modality = m
    df_odds_all = rbind(df_odds_all,df_odds)
  }
  
  # Calculate statistics 
  mu = df_odds_all %>% group_by(modality) %>% summarise(meanV = mean(V1))
  std =  df_odds_all %>% group_by(modality) %>% summarise(sdV = sd(V1))
  x =  df_odds_all %>% group_by(modality) %>% summarise(x = mean(odds))
  z =  left_join(mu,std,by = "modality")
  z =  left_join(z,x,by = "modality")
  z$z =  (z$x - z$meanV)/z$sdV
  z$p =  (1-pnorm(abs(z$z)))*2
  z$pfdr = p.adjust(z$p, method = "fdr")
  # export results
  write.xlsx(x = z, file="~output/gene_overlap_scz_permuted.xlsx")

