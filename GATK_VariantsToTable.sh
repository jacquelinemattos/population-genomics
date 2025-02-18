#!/bin/bash


#VariantsToTable
workdir=/home/jacqueline/data2/epidendrum_RNA


/home/ziv/gatk-4.1.8.1/gatk VariantsToTable \
     -V ${workdir}/final_vcf_file/epidendrum-vcf_only_biallelic_SNPs.vcf \
     -F CHROM -F POS -F TYPE -GF AD \
     -O epidendrum-vcf-table.table
