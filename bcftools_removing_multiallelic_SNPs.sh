#!/bin/bash

#Removing multi-allelic SNPs with bcftools

vcf_file=/home/jacqueline/data2/epidendrum_RNA/final_vcf_file/epidendrum-all-samples-vcf.vcf.gz

bcftools view --types snps -m 2 -M 2 $vcf_file > epidendrum-all-samples-vcf-only-biallelic-snps.vcf.gz
