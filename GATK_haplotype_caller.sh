#!/bin/bash


###Pass a list of sample names as argument to the script###

sample_names=$1


haplotypecaller(){

        local SEQ_NAME=$1

	#--java-options "-Xmx4g"
	/home/ziv/gatk-4.1.8.1/gatk --java-options "-Xmx4g" HaplotypeCaller \
	-R /home/jacqueline/data2/epidendrum_RNA/files/yahs.out_scaffolds_final_p_ctg.fa \
	-I /home/jacqueline/data2/epidendrum_RNA/picard/splitN/${SEQ_NAME}_duplMarked_rg_split.bam \
	-O /home/jacqueline/data2/epidendrum_RNA/gatk/haplotypecaller/${SEQ_NAME}_variants.g.vcf.gz \
	-ERC GVCF \
	--dont-use-soft-clipped-bases \
	--ploidy 2 \
	--standard-min-confidence-threshold-for-calling 20

}



###Exporting the function

export -f haplotypecaller



###Parallelizing the function to all the samples


parallel -j 12 haplotypecaller :::: $sample_names



#testing the function with --dryrun
#parallel --dryrun process_bams :::: $sample_names

