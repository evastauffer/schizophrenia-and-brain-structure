#### The genetic relaionships between brain structure and schizophrenia
#### Stauffer et al., 2023


# this script runs lava on positionally enriched regions (FUMA)
# Here shown for a specific region measured using SA


rm(list=ls())
options(stringsAsFactors = FALSE)


library(xlsx);library(LAVA);library(data.table);library(tidyverse)

###------------- Prepare DATA files-------------###

  # Import gene location from magma website, built 37
  gene_locations =fread("NCBI37.3.gene.loc")
  gene_locations = gene_locations[,c(2,3,4,6)]
  colnames(gene_locations) = c("CHR", "START", "STOP","GENE")
  
  ### We need to combine FUMA enrichments with LAVA locations
  
  # Locations of enriched positions for all phenotypes (chromosome, start,stop)
  position = fread("~/data/processed/position_enriched_region.txt")
  position = position[order(position$chromosome),]
  
  # Enriched positions
  enriched = read.xlsx("~/data/processed/fuma_enriched.xlsx",1)
  enriched$GeneSet = gsub("chr",'',enriched$GeneSet)
  
  #  Combine files
  enriched_position = left_join(position,enriched, by = c("location"="GeneSet"))
  
  # Create temporary names for merging
  enriched_position$CHR2 = ifelse(enriched_position$location=="1q21", "1q21",enriched_position$chromosome)
  enriched_position$CHR2 = ifelse(enriched_position$location=="1p34", "1p34",enriched_position$CHR2)
  enriched_position$CHR2 = ifelse(enriched_position$location=="2q33", "2q33",enriched_position$CHR2)
  enriched_position$CHR2 = ifelse(enriched_position$location=="2p23", "2p23",enriched_position$CHR2)
  enriched_position$CHR2 = ifelse(enriched_position$location=="2p22", "2p22",enriched_position$CHR2)
  enriched_position$CHR2 = ifelse(enriched_position$location=="3p21", "3p21",enriched_position$CHR2)
  enriched_position$CHR2 = ifelse(enriched_position$location=="3q22", "3q22",enriched_position$CHR2)
  
  ## LAVA locations
  loci = read.loci("~/data/raw/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
  
  # Create temporary names for merging
  loci$CHR2 = ifelse(loci$CHR==1 & loci$START >= position[position$location=="1q21",start] & loci$STOP <= position[position$location=="1q21",stop], "1q21",loci$CHR)
  loci$CHR2 = ifelse(loci$CHR==1 & loci$START >= position[position$location=="1p34",start] & loci$STOP <= position[position$location=="1p34",stop], "1p34",loci$CHR2)
  loci$CHR2 = ifelse(loci$CHR==2 & loci$START >= position[position$location=="2q33",start] & loci$STOP <= position[position$location=="2q33",stop], "2q33",loci$CHR2)
  loci$CHR2 = ifelse(loci$CHR==2 & loci$START >= position[position$location=="2p23",start] & loci$STOP <= position[position$location=="2p23",stop], "2p23",loci$CHR2)
  loci$CHR2 = ifelse(loci$CHR==2 & loci$START >= position[position$location=="2p22",start] & loci$STOP <= position[position$location=="2p22",stop], "2p22",loci$CHR2)
  loci$CHR2 = ifelse(loci$CHR==3 & loci$START >= position[position$location=="3p21",start] & loci$STOP <= position[position$location=="3p21",stop], "3p21",loci$CHR2)
  loci$CHR2 = ifelse(loci$CHR==3 & loci$START >= position[position$location=="3q22",start] & loci$STOP <= position[position$location=="3q22",stop], "3q22",loci$CHR2)
  
  # Combine enriched positions with lava positions
  enriched_positions_lava = left_join(enriched_position, loci, by = c("CHR2"="CHR2"))


####-------------------- LAVA: Process input-------------------- ####

  ## Set up for SA region 178
  input = process.input(input.info.file="input.info.txt",
                        sample.overlap.file="sample.overlap.txt",
                        ref.prefix="~/data/raw/g1000_eur/g1000_eur",
                        phenos=c("SA_plinkmeta178Nadapted","scz")
  )
  
  
  ## Create loop to run analysis for each locus
  #Store output for each enriched position

    temp = enriched_positions_lava[enriched_positions_lava$Phenotype=="SA",]
    temp = subset(temp,temp$START>=temp$start & temp$STOP<=temp$stop)
    n_loc = nrow(temp)
    #Unique loci
    locus_list = c(temp$LOC)
  
    # Univariate threshold is 0.05/number of tested loci
    univ_threshold = 0.05/n_loc
    
    # Store resluts for all loci 
    univ_all =NULL
    bivar_all = NULL
  
    for (l in locus_list) {
      #### Process locus ####
      # -------------------------------------------------
      locus = process.locus(loci[l,], input)
      print(l)
      #### Test local h2 and run rg ####
      # with univariate threshold
      rg = run.univ.bivar(locus, adap.thresh=NULL, univ.thresh=univ_threshold)
      univ = rg$univ
      univ$locus = l
      bivar =rg$bivar
      bivar$locus = l
      
      bivar_all = dplyr::bind_rows(bivar_all,bivar)
      univ_all = dplyr::bind_rows(univ_all,univ)
    }
    univ_all$bonferroni_threshold = univ_threshold
    univ_all$nloci = n_loc
    
    # Export univariate and bivariate analysis
    write.table(univ_all, paste0("~/output/lava_univ_SA_178.txt"), row.names = F)
    write.table(bivar_all, paste0("~output/lava_bivar_SA_178.txt"), row.names = F)


  ### Restrict to significant genetic correlations, add locus info
  enriched_positions_lava_temp = enriched_positions_lava[enriched_positions_lava$Phenotype=="SA",c(7:10)]
  
  data_biv = data_biv %>% 
    na.omit() %>% 
    dplyr::mutate(p_bon = 0.05/nrow(cur_data())) %>%
    subset(p<p_bon) %>%
    left_join(enriched_positions_lava_temp,c("locus"="LOC"))
  
 write.xlsx(data_biv,"~/output/lava_bivar_SA_178.xlsx")




