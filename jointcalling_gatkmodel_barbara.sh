#!/bin/bash
#SBATCH -t 23:00:00 -c 8 -o ./SCAFFOLD_jointcalling.out

export PATH="/home/bssleal/Epidendrum/14_variant_calling/gatk-4.2.0.0/:$PATH"
#Run GenomicsDBImport on GVCFs to consolidate
#parallelization: https://www.ibm.com/downloads/cas/ZJQD0QAL
gatk --java-options "-Djava.io.tmpdir=/store/bssleal/tmp/tmp_SCAFFOLD -XX:+UseParallelGC -XX:ParallelGCThreads=4" GenomicsDBImport \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C12-35_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C12-36_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C12-38_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C12-39_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C12-42_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C1-45_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C1-47_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C1-49_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C1-50_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C1-51_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C1-52_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C1-54_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C1-55_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C2-62_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C2-63_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C2-69_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C2-72_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C3-77_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C3-78_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C3-80_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C3-82_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C3-83_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C3-84_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C3-86_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C3-88_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C6-18_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C6-19_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C6-23_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C6-26_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C6-28_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C6-30_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C9-03_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C9-04_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C9-05_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C9-06_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C9-12_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C9-13_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C9-14_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C9-15_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C9-16_variants.g.vcf.gz \
-V /home/bssleal/Ipuberulus/08_gatk/4_haplotypecaller/C9-17_variants.g.vcf.gz \
--genomicsdb-workspace-path /home/bssleal/Ipuberulus/08_gatk/test_workspace_SCAFFOLD --intervals SCAFFOLD \
--tmp-dir /store/bssleal/tmp/tmp_SCAFFOLD

export PATH="/home/bssleal/Epidendrum/14_variant_calling/gatk-4.2.0.0/:$PATH"
#Run GenotypeGVCFs on the GDB workspace to produce final multisample VCF
gatk GenotypeGVCFs \
-R /store/bssleal/Ip_genome/som_cromossomos2/ISCPUB_v0.2_scaffolds.yahs_chroms_cor.fa \
-V gendb://test_workspace_SCAFFOLD \
-G StandardAnnotation \
-O /home/bssleal/Ipuberulus/08_gatk/5_genotypegvcf/SCAFFOLD.vcf


export PATH="/home/bssleal/Epidendrum/14_variant_calling/gatk-4.2.0.0/:$PATH"
#Filter variants with Phred-scaled p-value (FS) > 30 and Variant Confidence/Quality by Depth (QD) < 2.0
gatk VariantFiltration \
-R /store/bssleal/Ip_genome/som_cromossomos2/ISCPUB_v0.2_scaffolds.yahs_chroms_cor.fa \
-V /home/bssleal/Ipuberulus/08_gatk/5_genotypegvcf/SCAFFOLD.vcf -window 35 -cluster 3 \
--filter-name FS -filter-expression "FS > 30.0" \
--filter-name QD -filter-expression "QD < 2.0" \
-O /home/bssleal/Ipuberulus/08_gatk/5_genotypegvcf/filtered_SCAFFOLD.vcf

rm -r /store/bssleal/tmp/tmp_SCAFFOLD
rm -r /home/bssleal/Ipuberulus/08_gatk/test_workspace_SCAFFOLD




