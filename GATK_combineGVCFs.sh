#!/bin/bash


workdir=/home/jacqueline/data2/epidendrum_RNA

#Using CombineGVCFs instead of GenomicsDBImport because you can combine multiple intervals at once without building a GenomicsDB.
#CombineGVCFs is slower than GenomicsDBImport though.




for gvcf in ${workdir}/gatk/haplotypecaller/*_variants.g.vcf.gz; do

	gvcf_args+="--variant $gvcf " 
done

		
echo running CombineGVCFs on $gvcf_args

/home/ziv/gatk-4.1.8.1/gatk CombineGVCFs \
-R ${workdir}/files/yahs.out_scaffolds_final_p_ctg.fa \
$gvcf_args \
-O epidendrum_all_GVCFs_cohort.g.vcf.gz







