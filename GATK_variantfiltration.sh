#!/bin/bash

#VariantFiltration
#https://gatk.broadinstitute.org/hc/en-us/articles/13832750065947-VariantFiltration
#Choosing the appropriate filters in the VariantFiltration step of GATK for population genomics involves a careful consideration of various parameters to effectively separate true genetic variants from sequencing artifacts and noise

workdir=/home/jacqueline/data2/epidendrum_RNA


#Filter variants with Phred-scaled p-value (FS) > 30 and Variant Confidence/Quality by Depth (QD) < 2.0

#/home/ziv/gatk-4.1.8.1/gatk VariantFiltration \
#-R ${workdir}/files/yahs.out_scaffolds_final_p_ctg.fa \
#-V ${workdir}/final_vcf_file/epidendrum-all-samples-vcf.vcf.gz --window 35 --cluster 3 \
#-O ${workdir}/final_vcf_file/epidendrum-final-filtered-vcf.vcf.gz \
#--filter-name "FS" \
#--filter-expression "FS > 30.0" \
#--filter-name "QD" \
#--filter-expression "QD < 2.0" 


/home/ziv/gatk-4.1.8.1/gatk SelectVariants \
	-R ${workdir}/files/yahs.out_scaffolds_final_p_ctg.fa \
	-V ${workdir}/final_vcf_file/epidendrum-final-filtered-vcf.vcf.gz \
	--select-type-to-include SNP \
	--restrict-alleles-to BIALLELIC \
	--exclude-filtered \
	-O ${workdir}/final_vcf_file/epidendrum-final-variants-removed.vcf.gz


