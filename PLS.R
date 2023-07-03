#### The genetic relaionships between brain structure and schizophrenia
#### Stauffer et al., 2023

# This script runs PLS on brain structure and SCZ

rm(list=ls())
options(stringsAsFactors = FALSE)


library(plsdepot);library(ggplot2);library(svMisc);library(data.table);
library(tidyr);library(tidyverse);library(xlsx)

ensg_genename = fread("/mnt/b2/home4/arc/ems206/cortex_project_04_03_2021/data/raw/ensg_genename.txt")

### Prepare SCZ data
  disorder_df = fread("~/data/processed/h_magma/disorders_traits/FB_SCZ.genes.out")
  disorder_df = disorder_df[,c("GENE","ZSTAT")]
  names(disorder_df)[2] = disorder
  
### Create storage for loop
  screeplots = NULL
  variance_all = NULL
  
### Set up loop for PLS analysis
  modality_list = c("SA","CT","ICVF")
  for (modality in modality_list) {
    print(modality)
    
    ### Prepare Data
    #Import brain region data
    roi_df = fread(paste0("~/data/processed/h_magma/FB_combined/FB_",modality,"_combined_regions.txt"))
    
    # Combine data
    datatable = merge(roi_df, disorder_df, by = "GENE")
    datatable = na.omit(datatable)
    datatable = datatable %>% remove_rownames %>% column_to_rownames(var="GENE")
    datatable = as.matrix(datatable)
    
    # Define predictors and outcome variables
    X = datatable[,1:180]
    Y = datatable[,181,drop=FALSE] # Disorder 
    
    nperm = 1000 # Define number of permutation iterations
    nboot = 1000 # Define number of bootstrap iterations
    
    ### Step 1: Basic PLS Model  
    pls.model = plsreg1(predictors = X, response = Y, comps=5, crosval = TRUE) 
    
    # Save number of components from original model
    ncomp = 5
    
    # Variance Explained
    variance= as.data.frame(pls.model$R2)
    variance$Component = 1:5
    names(variance)[1] = "R2"
    
    # Extract T and U scores for later
    x.scores = as.data.frame(pls.model$x.scores) # get T scores
    y.scores = as.data.frame(pls.model$y.scores) # get U scores
    scores = cbind(x.scores,y.scores)
    scores=rownames_to_column(scores,"GENE")
    scores = left_join(scores,ensg_genename, by= c("GENE"="Gene_stable_ID"))
    scores = scores[,c(1,12, 2:11)]
    write.table(scores, paste("~/output/pls/SCZ_",modality,"_t_u_component.txt",sep = ''),row.names = F)
    
    
    ### Step 2: Significance of PLS Components
  
    # Does the model explain more variance than expected by chance?
    perm.var.expl  = NULL
    for (perm in 1:nperm) {
      print(paste("SCZ_",modality,"_permute_",perm, sep = ""))
      # Permute Y randomly
      order = sample(1:dim(Y)[1]) 
      # Recalculate PLS with permuted Y
      perm.model = plsreg1(X, Y[order,], comps = 5, crosval = TRUE)
      # Save explained variance in Y
      var.expl = perm.model$R2
      perm.var.expl = rbind(var.expl,perm.var.expl)
    }
    # Calculate PLS component p-values
    component.pvalues = vector(length = ncomp)
    for (c in 1:ncomp) {
      # Is the real variance explained higher than in the permutation distribution?
      component.pvalues[c] = sum(pls.model$R2[c]<perm.var.expl[,c])/nperm
    }
    

    #save results
    components_var_p = cbind(variance,component.pvalues)
    components_var_p = components_var_p[,c(2,1,3)]
    write.table(components_var_p, paste0("~/output/FB_",modality,"_SCZ_components_var_p.txt"), row.names = F)
    
    variance2 = components_var_p
    variance2$phenotype = modality
    variance_all = rbind(variance_all, variance2)
    
    # Export results
    write.table(perm.var.expl, paste0("~/output/FB_",modality,"_SCZ_perm.var.expl.txt"), row.names = F)
    
    ###  Step 3: Significance Testing of PLS Weights 
    # We want the PLS components to correlate (not anti-correlate) with X
    # Therefore we realign each component if it anticorrelates
    ncomp=5
    
    R1 = cor(pls.model$x.scores,X)
    for (c in 1:ncomp) {
      if (abs(max(R1[c,])) < abs(min(R1[c,]))) {
        pls.model$x.scores[,c] = -pls.model$x.scores[,c]
        pls.model$mod.wgs[,c] = -pls.model$mod.wgs[,c]
      }
    }
    # Bootstrapping loop
    npred = dim(X)[2]
    boot.weights = array(NA,dim=c(npred,ncomp,nboot))
    pb = txtProgressBar(min = 0, max = length(nboot), initial = 0)
    for (b in 1:nboot) {
      print(paste(disorder,"_",modality,"_bootstrap_",b, sep = ""))
      # Random order
      set.seed(b)
      bootorder = sample(1:nrow(Y),nrow(Y),TRUE)  #in order to get stable weights we change the order of the rows, order of rows can impact PLS slightly
      # Recalculate PLS with bootstrapped X and Y
      Xrand = X[bootorder,]  #same order for X and Y
      Yrand = Y[bootorder,]
      boot.model = plsreg1(Xrand,Yrand, comps = ncomp, crosval = FALSE)
      for (c in 1:ncomp){
        # Realign each bootstrap component so it correlates positively with X
        if (cor(boot.model$mod.wgs[,c],pls.model$mod.wgs[,c])<0){
          boot.model$mod.wgs[,c] = -boot.model$mod.wgs[,c]
        }
      }
      # Save bootstrap weights
      boot.weights[,,b] = boot.model$mod.wgs
    }
    
    corr.weights = weight.pvals = weights.pvals.fdr = matrix(NA,nrow=npred,ncol=ncomp)
    pval = 0.05 # Desired significance level of weights
    for (c in 1:ncomp) {
      boot.sd = apply(boot.weights[,c,],1,sd)  #calculate standard deviation for each component
      corr.weights[,c] = pls.model$mod.wgs[,c]/boot.sd # translate real values into standard deviations
      weight.pvals[,c] = 2*pnorm(-abs(corr.weights[,c]))
      weights.pvals.fdr[,c] = p.adjust(weight.pvals[,c] ,method = 'fdr')
      
      df_boot = data.frame(name = colnames(X), boot.weight = corr.weights[,c], 
                           pvalue = weight.pvals[,c], fdr.pvalue = weights.pvals.fdr[,c])
      
      df.asc = df_boot[order(df_boot$boot.weight),]
      df.asc$order = order(df_boot$boot.weight)
      write.csv(df.asc, paste("~/output/pls/SCZ_",modality,'_PLS_Weights_Component_',c,'_ascending.csv', sep = ""),
                row.names = FALSE)
    }
    # boot.weight: those are PLS values, they indicate how much a region is associated wit SCZ in standard deviations
    
  } 
  
  # Save variance explained for all phenotypes per disorder
  variance_all = variance_all[,c(4,1:3)]
  write.xlsx(variance_all,"~/output/pls/SCZ_variance_permuted.xlsx", sep = "")
  write.table(variance_all,"~/output/pls/SCZ_variance_permuted.txt", row.names = F)
  


