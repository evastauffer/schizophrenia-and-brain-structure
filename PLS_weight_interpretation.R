#### The genetic relaionships between brain structure and schizophrenia
#### Stauffer et al., 2023

# This script plots Z scores from SCZ and brain regions with highest PLS weight

library(data.table); library(ggplot2); library(ggpubr)

# Set path
inpath = '/mnt/b2/home4/arc/ems206/cortex_project_04_03_2021/output/pls/'

##### Import PLS weights for component 
pls_weights_component_1 = fread(paste(inpath,"SCZ_SA_PLS_Weights_Component_1_ascending.csv", sep = ""))
pls_weights_component_1$component = "Component 1"

##### Import Z scores for schizophrenia ######### 
disorder_df = fread("/mnt/b2/home4/arc/ems206/cortex_project_04_03_2021/data/processed/h_magma/disorders_traits/FB_SCZ.genes.out")
disorder_df = disorder_df[,c("GENE","ZSTAT")]
names(disorder_df)[2] = "SCZ"

##### Import brain region data for SA
roi_df = fread("/mnt/b2/home4/arc/ems206/cortex_project_04_03_2021/data/processed/h_magma/FB_combined/FB_SA_combined_regions.txt")

##### Combine data
datatable = merge(roi_df, disorder_df, by = "GENE")
datatable = na.omit(datatable)
datatable = datatable %>% remove_rownames %>% column_to_rownames(var="GENE")
datatable = as.matrix(datatable)

##### Find region with highest PLS weight
big_weight=pls_weights_component_1[order(pls_weights_component_1$boot.weight, decreasing = T),]
big_weight = big_weight[1,big_weight$name]

##### To understand what directionality means we are plotting a region with high value against genetic effects of disorder
# Plot brain region with high value
correlation = cor.test(datatable[,paste(big_weight)], y=datatable[,paste("SCZ")])

p= ggplot(datatable, aes(x=datatable[,paste(big_weight)], y=datatable[,paste("SCZ")])) + 
      geom_point(size=0.5, alpha = 0.8)+
      geom_smooth(method=lm,colour = "red")+ 
      scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
      theme_bw()+
      theme(axis.text=element_text(size=16),axis.title=element_blank())

ggsave(p,filename = "~/output/plots/subset/SA_big_pls_weight.pdf",height = 2, width = 2.5)





