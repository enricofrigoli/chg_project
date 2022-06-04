# Somatic variant calling

Somatic variant calling is perfomed with Mutect2 on analysis-ready bam files (deduplicated and recalibrated). 


```bash
gatk Mutect2 -R ../Annotations/human_g1k_v37.fasta -I ../task2/Tumor.recal.bam \
-I ../task2/Control.recal.bam -L ../Captured_Regions.bed -O somatic_unfiltered_mt2.vcf
```

A needed tool to filter the raw callset found by Mutect2 is `FilterMutectCalls`:

```bash
gatk FilterMutectCalls -R ../Annotations/human_g1k_v37.fasta -V somatic_unfiltered_mt2.vcf -O filtered_mt2.vcf
```


