###SNP relate###
#Making PCAs from SNP data#

#Install package#
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")


#Load packages#
library(SNPRelate)
library(ggplot2)
library(dplyr)


#Load VCF#
pruned_vcf <- "C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics_RNAseq/Final_VCF/PLINK_LD/vcf_pruned_by_plink.vcf"

#Load sample info#
sample_info <- read.csv("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics_RNAseq/Sample_info/trinity_efulgens_table.csv")


snpgdsVCF2GDS(pruned_vcf, "snps.gds")
geno_file <- snpgdsOpen("snps.gds")

#Do PCA#
pca <- snpgdsPCA(geno_file, num.thread=2, autosome.only = FALSE) #pca with no LD pruning
pca$sample.id <- gsub("_.*","", pca$sample.id)
