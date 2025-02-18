/home/ziv/gatk-4.1.8.1/gatk --java-options "-Xmx4g" HaplotypeCaller \
        -R /home/jacqueline/data2/epidendrum_RNA/files/yahs.out_scaffolds_final_p_ctg.fa \
        -I /home/jacqueline/data2/epidendrum_RNA/picard/splitN/203_S76_duplMarked_rg_split.bam \
        -O /home/jacqueline/data2/epidendrum_RNA/gatk/haplotypecaller/203_S76_variants.g.vcf.gz \
        -ERC GVCF \
        --dont-use-soft-clipped-bases \
        --ploidy 2 \
        --standard-min-confidence-threshold-for-calling 20
