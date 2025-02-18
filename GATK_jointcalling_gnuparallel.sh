#!/bin/bash


###Pass a list of sample names as argument to the script###

sample_names=$1
workdir=/home/jacqueline/data2/epidendrum_RNA

#Using CombineGVCFs instead of GenomicsDBImport because you can combine multiple intervals at once without building a GenomicsDB.
#CombineGVCFs is slower than GenomicsDBImport though.


combineGVCFs(){

	local SEQ_NAME=$1

	/home/ziv/gatk-4.1.8.1/gatk CombineGVCFs \ 
	-R ${workdir}/files/yahs.out_scaffolds_final_p_ctg.fa \ 
	--variant ${workdir}/gatk/haplotypecaller/${SEQ_NAME}_variants.g.vcf.gz \ 
	-O epidendrum_all_GVCFs_cohort.g.vcf.gz

}



genotypeGVCFs(){

	/home/ziv/gatk-4.1.8.1/gatk --java-options "-Xmx4g" GenotypeGVCFs \ 
	-R reference-genome.fa \ 
	-V gendb://test_workspace_SCAFFOLD \ 
	-O output.vcf

}




###Exporting the function
export -f genomicsDBImport



###Parallelizing the function to all the samples
parallel -j 12 genomicsDBImport :::: $sample_names



