#!/bin/bash


#GenotypeGVCFs
#https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs


workdir=/home/jacqueline/data2/epidendrum_RNA

/home/ziv/gatk-4.1.8.1/gatk --java-options "-Xmx4g" GenotypeGVCFs \
        -R ${workdir}/files/yahs.out_scaffolds_final_p_ctg.fa \
        -V ${workdir}/cohort_variants_gvcf/epidendrum_all_GVCFs_cohort.g.vcf.gz \
   	-O epidendrum-all-samples-vcf.vcf.gz
