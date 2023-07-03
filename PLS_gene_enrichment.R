#### The genetic relaionships between brain structure and schizophrenia
#### Stauffer et al., 2023

# This script plots T and U scores and performs gene enrichments 

rm(list=ls())
options(stringsAsFactors = FALSE)

library(ggplot2);library(data.table);library(gprofiler2);library(xlsx);library(nlme);
library(ggrepel);library(Hmisc);library(dplyr);library(ensembldb);
library(EnsDb.Hsapiens.v86);library(scico)


######### Prepare what's needed for loop ######### 
# Set path
inpath = '~/output/pls/' # Path to output directory

# Loss of function genes downoaded from https://gnomad.broadinstitute.org/downloads#v2-constraint file: pLoF Metrics by Gene TSV
LoF = read.table("/mnt/b2/home4/arc/ems206/cortex_project_04_03_2021/data/raw/gnomad.v2.1.1.lof_metrics.by_gene.txt", sep = "\t", header = T)
LoF = LoF[,c("gene_id","pLI" , "oe_lof_upper_bin", "gene_length","start_position","end_position",  "chromosome" )]
LoF = LoF[!is.na(LoF$oe_lof_upper_bin), ]
LoF$oe_lof_upper_bin2= ifelse(LoF$oe_lof_upper_bin == 0,1,0) #code first bin as 1, those are constraint, rest as zero
length = LoF[,c("gene_id","gene_length")]

# Cell types 
load("~/data/celltype_geschwind/Mean_data.RData")
celltype = merge(specific_mean,length, by.x = "GENE", by.y = "gene_id")

# Get genetic position information 
edb <- EnsDb.Hsapiens.v86
g <- genes(edb)
EnsDb = as.data.frame(g)
EnsDb = EnsDb[EnsDb$gene_biotype=="protein_coding",] 
EnsDb = EnsDb[,c(1:3,6,7)]
names(EnsDb)[5] = "Gene_name"


######### Loop for enrichment analysis ######### 
threshold_list = c(.01, .03)
modality_list = c("SA","CT","ICVF")

for (threshold in threshold_list) {
    ### Empty list for results
    genes_per_modality = NULL
    df_plot_all = NULL
    gostres1_df_all = NULL
    all_LoF = NULL
    celltypes_all_modalities = NULL
    
    for (modality in modality_list) {
      ### Import T and U values ###
      scores = fread(paste(inpath,"SCZ_",modality,"_t_u_component.txt",sep = ''))
      
      # Select T and U scores for component 1
      T1= scores$t1
      U1=scores$u1
      gene_names = scores[,1, drop=F]
      
      ### Leave-one-out: for each gene check whether exclusion decreases the correlatin between T and U
      cor.loo = vector(length=length(T1))
      for(i in 1:length(T1)){
        cor.loo[i] = cor.test(T1[-i],U1[-i])$estimate
      }
      
      diff.loo = cor.test(T1,U1)$estimate - cor.loo
      df = data.frame(GENE=gene_names$GENE,t=T1, u= U1, loo = diff.loo)

     
      ### Selecting genes based on thresholds
      df.thresho= df %>% top_frac(wt=loo,threshold)
      genes_impacting_correlation = left_join(df.thresho,EnsDb, by= c("GENE"="gene_id"))
      genes_impacting_correlation = genes_impacting_correlation[order(genes_impacting_correlation$loo,decreasing = T1),]
      
      # Save them for later
      n_genes = genes_impacting_correlation[,c("Gene_name","loo")]
      n_genes = n_genes[order(n_genes$loo,decreasing=T),"Gene_name",drop=F]
      n_genes$modality = modality
      genes_per_modality = rbind(genes_per_modality,n_genes)
      
      ### Prepare data for T and U score Plot further down
      df_plot = left_join(df,genes_impacting_correlation, by = "GENE")
      df_plot$Gene_Name = df_plot$Gene_name
      df_plot$Gene_name = ifelse(is.na(df_plot$Gene_name),0,1)
      df_plot$Gene_name = as.factor(df_plot$Gene_name)
      df_plot$t.x = scale(df_plot$t.x)
      df_plot$u.x = scale(df_plot$u.x)
      df_top = df_plot %>% dplyr::arrange(desc(loo.x)) %>% 
        dplyr::slice(1:10)
      df_plot$label = if_else(df_plot$Gene_Name %in% df_top$Gene_Name,  
                              df_plot$Gene_Name, NULL)
      df_plot$modality = modality
      df_plot_all = rbind(df_plot_all,df_plot)
      
      ### Run functional enrichment analysis GO ###
      df2  = df[order(df$loo, decreasing = T),]
      gostres1 = gprofiler2::gost(query = df.thresho$GENE, organism = "hsapiens",
                                  ordered_query =F,
                                  user_threshold = 0.05,
                                  custom_bg = df2$GENE,
                                  correction_method = c("gSCS"),
                                  sources = c("GO:BP"))
      
      
      gostres1_df= gostres1$result
      gostres1_df$modality = modality
      gostres1_df_all = rbind(gostres1_df_all,gostres1_df)
      
      #### Enrichment for loss of function intolerant genes  ###
      #Prepare data
      scores$relevant = ifelse(scores$GENE %in% genes_impacting_correlation$GENE, 1, 0) # Set pleiotropic genes to 1
      scores = scores[,c("GENE","Gene_name","relevant")]
      df_LoF = merge(scores,LoF, by.x="GENE", by.y ="gene_id")
      # run model
      model = glm(relevant ~ oe_lof_upper_bin2 + log(gene_length), data = df_LoF, family = binomial)
      s=summary(model)
      coef  =coefficients(s)[2,]
      coef$disorder = disorder
      coef = as.data.frame(coef)
      
      odds = exp(cbind("Odds ratio" = coef(model), confint.default(model, level = 0.95)))[2,]
      odds$disorder = disorder
      odds = as.data.frame(odds)
      
      out_LoF = cbind(coef,odds)
      out_LoF = out_LoF[,c(5,1:4,6:8)]
      out_LoF$modality = modality
      all_LoF = rbind(all_LoF, out_LoF)
      
      ### Enrichment in celltypes for pleiotropic genes ###
      # Prepare data
      celltype_specific = celltype[celltype$GENE %in% scores$GENE,]
      celltype_specific$significant = ifelse(celltype_specific$GENE %in% genes_impacting_correlation$GENE, 1, 0)
      
      # Run enrichment
      table_celltypes = NULL
      celltypes = c("vRG","oRG","PgS","PgG2M", "IP","ExN","ExM", "ExM-U", "ExDp1","ExDp2","InMGE","InCGE","OPC","End","Per","Mic")
      
      for (type in celltypes){
        beta_2 = summary(lm(log(celltype_specific[[type]]/celltype_specific$total_mean) ~ celltype_specific$significant + log(celltype_specific$gene_length)))
        estimate = as.data.frame(cbind(type, t(beta_2$coefficients[2,]), disorder))
        table_celltypes = rbind(table_celltypes, estimate)
      }
      
      table_celltypes$pfdr = p.adjust(table_celltypes$`Pr(>|t|)`, method = "fdr")
      table_celltypes$modality = modality
      celltypes_all_modalities = rbind(celltypes_all_modalities,table_celltypes)
    }
  

    ### Export genes per modality and threshold
    write.xlsx2(x = genes_per_modality, file = paste(inpath,  "SCZ_top_",threshold,"_genesper_modality_subset.xlsx", sep = ""))

    ### Plot T and U scores, genes colored by delta correlation and top ten genes annotated
    df_plot_all$modality = factor(df_plot_all$modality, levels = c("SA","CT","ICVF"))
    df_plot_all$loo.x = scale(df_plot_all$loo.x)

    p_annotated = ggplot(df_plot_all, aes(x=t.x, y=u.x)) +
      geom_point(aes(color=loo.x),size=2.5)+
      scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                            high = "red", space = "Lab" ,name = "Change in r")+
      geom_smooth(method=lm, se=FALSE, size = 0.5)+
      geom_text_repel(aes(label = label), show.legend = FALSE,
                      size =4 ,colour = "black",nudge_x = .15,
                      box.padding = 0.6,
                      segment.ncp = 1,
                      segment.linetype = 2,
                      min.segment.length = unit(0, 'lines'),
                      fontface = "bold",
                      segment.size = 0.3,
                      segment.angle = 20,
                      max.overlaps = 30,
                      force = 20,
                      direction = "both")+
      guides(shape="none")+
      theme_bw()+theme(legend.key.size = unit(0.5, 'cm'),strip.background = element_blank(),panel.spacing = unit(1, "lines"),
                       strip.text.x = element_blank(),axis.text=element_text(size=18),axis.title = element_blank())+
      facet_wrap(~modality,ncol = 3)
    ggsave(p_annotated,filename = paste("~/output/plots/subset/SCZ_PLS_annotated_changeinr_",threshold,".pdf", sep=""),height = 2.5, width = 14)
  

    ### Export LOF results and make plot
    write.xlsx2(x = all_LoF, file = paste(inpath, "SCZ_PLS_LoF_enrichment_subset_top_",threshold,".xlsx", sep = ""), row.names = F)
    pLI_scico <- scico(2, palette = "devon", begin = 0.3, end = 0.75)
    
    p_LoF<- ggplot(data=all_LoF, aes(x=modality, y=Odds.ratio, ymin=X2.5..,ymax= X97.5..)) +
      geom_pointrange(size = 1, colour = ifelse(all_LoF$Pr...z..<0.05,pLI_scico[1],pLI_scico[2]), fatten = 3) + 
      geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
      coord_flip() +  # flip coordinates (puts labels on y axis)
      xlab("Modality") + ylab("OR (95% CI)") +
      theme_bw()
    ggsave(p_LoF,filename = paste0("~/output/plots/subset/SCZ_top_",threshold,"_PLS_LoF_enrichment.pdf"),height = 6, width = 4)
    
    
    ### Export celltype enrichment results and make plot
    celltypes_all_modalities$modality = factor(celltypes_all_modalities$modality, levels = c("ICVF","CT","SA"))
    celltypes_all_modalities$type = factor(celltypes_all_modalities$type, levels = c("vRG","oRG", "PgS","PgG2M","IP","ExN","ExM","ExM-U","ExDp1","ExDp2","InMGE", "InCGE", "OPC","End", "Per","Mic"))
    celltypes_all_modalities$pfdrlog = -log10(celltypes_all_modalities$pfdr)
    celltypes_all_modalities$stars <- cut(celltypes_all_modalities$pfdr, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 
    write.xlsx2(x = celltypes_all_modalities, file = paste(inpath,  "SCZ_top_",threshold,"_PLS_celltype_enrichment_subset.xlsx", sep = ""), row.names = F)
    
    p_celltype = ggplot(celltypes_all_modalities, aes(x=type,y=modality)) +
      geom_tile(aes(fill = `t value`)) + 
      xlab("Cell type") + ylab("Metric")+
      geom_text(aes(label = stars)) +
      scale_fill_gradient2(low ="blue",high="red")+ 
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    ggsave(p_celltype,filename = paste0("~/output/plots/subset/SCZ_top_",threshold,"_PLS_celltype_enrichment.pdf"),height = 3, width = 8)
    
    
    ## Export GO enrichment
    write.xlsx2(x = gostres1_df_all, file = paste(inpath,  "SCZ_top_",threshold,"_PLS_GO_enrichment.xlsx", sep = ""), row.names = F)
    gostres_sig = gostres1_df_all[gostres1_df_all$p_value<0.05,]
    write.xlsx2(x = gostres_sig, file = paste(inpath,  "SCZ_top_",threshold,"_PLS_GO_enrichment_sig.xlsx", sep = ""), row.names = F)
}
