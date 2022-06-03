# Somatic variant calling

Somatic variant calling is perfomed with Mutect2 on analysis-ready bam files (deduplicated and recalibrated). 


```bash
gatk Mutect2 -R ../Annotations/human_g1k_v37.fasta -I ../task2/Tumor.recal.bam \
-I ../task2/Control.recal.bam -L ../Captured_Regions.bed -O somatic_mt2.vcf
```






