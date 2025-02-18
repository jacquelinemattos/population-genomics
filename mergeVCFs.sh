#!/bin/bash
#SBATCH -t 23:00:00 -c 1 -o ./merge_vcfs.out

java -jar picard.jar MergeVcfs \
I=/home/bssleal/Ipuberulus/08_gatk/5_genotypegvcf/filtered_scaffold_1.vcf \
I=/home/bssleal/Ipuberulus/08_gatk/5_genotypegvcf/filtered_scaffold_1a.vcf \
I=/home/bssleal/Ipuberulus/08_gatk/5_genotypegvcf/filtered_Scaffold_1b.vcf \
I=/home/bssleal/Ipuberulus/08_gatk/5_genotypegvcf/filtered_scaffold_2.vcf \
I=/home/bssleal/Ipuberulus/08_gatk/5_genotypegvcf/filtered_scaffold_3.vcf \
I=/home/bssleal/Ipuberulus/08_gatk/5_genotypegvcf/filtered_scaffold_4.vcf \
I=/home/bssleal/Ipuberulus/08_gatk/5_genotypegvcf/filtered_scaffold_5.vcf \
I=/home/bssleal/Ipuberulus/08_gatk/5_genotypegvcf/filtered_scaffold_6.vcf \
I=/home/bssleal/Ipuberulus/08_gatk/5_genotypegvcf/filtered_scaffold_7.vcf \
I=/home/bssleal/Ipuberulus/08_gatk/5_genotypegvcf/filtered_scaffold_8.vcf \
I=/home/bssleal/Ipuberulus/08_gatk/5_genotypegvcf/filtered_scaffold_9.vcf \
I=/home/bssleal/Ipuberulus/08_gatk/5_genotypegvcf/filtered_scaffold_10.vcf \
I=/home/bssleal/Ipuberulus/08_gatk/5_genotypegvcf/filtered_scaffold_11.vcf \
O=/home/bssleal/Ipuberulus/08_gatk/5_genotypegvcf/ipuberulus_prefilteredSNPs_scaffolds1-11.vcf

export PATH=/home/bssleal/htslib-1.13:$PATH
bgzip -c ipuberulus_prefilteredSNPs_scaffolds1-11.vcf > ipuberulus_prefilteredSNPs_scaffolds1-11.vcf.gz
