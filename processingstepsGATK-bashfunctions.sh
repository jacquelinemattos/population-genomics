#!/bin/bash

###Parallelizing the pre-processing commands before GATK HaplotypeCaller###


sample_names=$1


process_bams(){

	local SEQ_NAME=$1
	
	#common_dir=/home/jacqueline/data20/epidendrum_RNA/star_results
	#creating an index for the bams, which is required for GATK	
	#samtools index -b ${common_dir}/2nd_pass/${SEQ_NAME}_2ndpass_Aligned.sortedByCoord.out.bam	



	#marking duplicates with picard
	java -jar /home/superstar/rahul/anaconda3/envs/snpcall/share/picard-2.21.9-0/picard.jar MarkDuplicates I=/home/jacqueline/data2/epidendrum_RNA/star_2ndpass_bams/${SEQ_NAME}_2ndpass_Aligned.sortedByCoord.out.bam O=/home/jacqueline/data2/epidendrum_RNA/picard/${SEQ_NAME}_duplMarked.bam M=/home/jacqueline/data2/epidendrum_RNA/picard/${SEQ_NAME}_duplMarked.metrics


	#adding read groups using picard
	java -jar /home/superstar/rahul/anaconda3/envs/snpcall/share/picard-2.21.9-0/picard.jar AddOrReplaceReadGroups I=/home/jacqueline/data2/epidendrum_RNA/picard/${SEQ_NAME}_duplMarked.bam O=/home/jacqueline/data2/epidendrum_RNA/picard/${SEQ_NAME}_duplMarked_rg.bam RGLB=${SEQ_NAME}_lib RGPL=illumina RGPU=none RGSM=${SEQ_NAME}

	rm /home/jacqueline/data2/epidendrum_RNA/picard/${SEQ_NAME}_duplMarked.bam


	#SplitNCigarReads function that will split reads that contain Ns in their cigar string (e.g. spanning splicing events in RNAseq data).
        /home/superstar/rahul/anaconda3/envs/snpcall/share/gatk4-4.1.4.1-1/gatk SplitNCigarReads -R /home/jacqueline/data2/epidendrum_RNA/files/yahs.out_scaffolds_final_p_ctg.fa -I /home/jacqueline/data2/epidendrum_RNA/picard/${SEQ_NAME}_duplMarked_rg.bam -O /home/jacqueline/data2/epidendrum_RNA/picard/splitN/${SEQ_NAME}_duplMarked_rg_split.bam


	rm /home/jacqueline/data2/epidendrum_RNA/picard/${SEQ_NAME}_duplMarked_rg.bam


	#creating an index for the bams, which is required for GATK
        samtools index -b /home/jacqueline/data2/epidendrum_RNA/picard/splitN/${SEQ_NAME}_duplMarked_rg_split.bam

}


# Export the function

export -f process_bams


# Parallelizing the function to all the samples

parallel -j 20 process_bams :::: $sample_names



#testing the function with --dryrun
#parallel --dryrun process_bams :::: $sample_names





