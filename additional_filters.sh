#!/bin/bash 

###Additional filters after the Variant Filtration step from GATK### 

workdir=/home/jacqueline/data2/epidendrum_RNA

# minor allele frequency and missingness filters for vcftools
#MAF=0.05
#MISS=0.75

#vcftools --gzvcf ${workdir}/final_vcf_file/epidendrum-final-variants-removed.vcf.gz \
#	--maf $MAF \
#	--max-missing $MISS \
#	--recode --stdout | gzip -c \
#	> ${workdir}/final_vcf_file/epidendrum-final-variants-removed.maf_miss.vcf.gz




# trying to run with bcftools 
bcftools +fill-tags ${workdir}/final_vcf_file/epidendrum-final-variants-removed.vcf.gz -- -t MAF,F_MISSING | bcftools filter -i 'MAF>=0.05 & F_MISSING<0.25' -Oz -o ${workdir}/final_vcf_file/epidendrum-final-variants-removed_filtered_bcftools.vcf.gz











