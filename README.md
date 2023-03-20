# Genetic-relationship-schizophrenia-and-brain-structure

Author: Eva-Maria Stauffer

This repository contains the code for the main analyses of the manuscript "The genetic relationships between brain structure and schizophrenia" by Eva-Maria Stauffer, Richard A.I. Bethlehem, Lena Dorfschmidt, Hyejung Won, Varun Warrier and Edward Bullmore.

Details can be found in the manuscript:

# Data
Summary statistics on cortical MRI phenotypes are available for download on https://portal.ide-cam.org.uk/overview/ and are described in detail in Warrier et al., 2023 (https://www.biorxiv.org/content/10.1101/2022.09.08.507084v1). Summary statistics for schizophrenia can be accessed from the Psychiatric Genomics Consortium https://pgc.unc.edu/for-researchers/download-results/ and are described in Trubetskoy et al., 2022 (https://www.nature.com/articles/s41586-022-04434-5). Imaging data used to construct structural covariance matrices may be requested through the
UK Biobank database https://www.ukbiobank.ac.uk/. Spatiotemporal gene expression data can be accessed from PsychENCODE http://development.psychencode.org/. Cell type specific expression data can be downloaded from http://solo.bmap.ucla.edu/shiny/webapp/. Information on constraint genes can be accessed from gnomAD https://gnomad.broadinstitute.org/downloads#v2-constraint.

# Software
Analyses were conducted using R version 3.5.2. We used H-MAGMA for SNP to gene mapping. Details on how to use H-MAGMA and required reference files can be found here https://github.com/thewonlab/H-MAGMA.

# Scripts

1) FB_hmagma.sh : sample code for snp-to-gene mapping. Shown for region one measured using surface area. Annotation files can be accessed here https://github.com/thewonlab/H-MAGMA/tree/master/Input_Files. 

2) genetic_data_correlation_matrices.R: this script calculates correlation matrices for genetic data (adult and fetal brain)

2) MRI_associated_genes.R: Identfies and characterises genes associated with MRI phenotypes (GO analysis and devvelopmental expresison profiles)

3) SCZ_associated_genes.R: Identfies and characterises genes associated with schizophrenia (GO analysis and devvelopmental expresison profiles)

4) Overlap_SCZ_MRI_genes.R: Extracts overlapping genes between MRI modality and SCZ, tests whether overlap is significant

5) PLS.R: Runs PLS on schizophrenia and brain structure

6) PLS_weight_interpretation.R: This script plots Z scores from SCZ and brain regions with highest PLS weight

7) PLS_weight: This script plots PLS weights on brain and extracts T and U scores for genes 

8) PLS_weights_enrichment_mesulam.R : This script performs enrichments of PLS weights in Mesulam class

9) correlation_sc_gs_heatmaps.R: correlates structural covariance and genetic similarity, creates heatmaps

10) cophenetic_correlations_within_phenotypes.R : caluclates cophenetic correlations 

11) louvain_modules.R : identfies Louvain modules

12) hubs.R: calculates hubness for imaging and genetic data, compares hubness to PLS weights

13) intra_inter_modular_degree.R: calculates inter- and intra-modular degree, comparison to PLS weights

14) PLS_gene_enrichment.R: This script plots T and U scores and performs gene enrichments 

15) LAVA.R: This script runs LAVA on positionally enriched regions (FUMA)


FUMA enrichments were analysed using the FUMA web interface https://fuma.ctglab.nl/
