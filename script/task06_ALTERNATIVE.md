# Somatic variant calling

Somatic variant calling is perfomed with Mutect2 on analysis-ready bam files (deduplicated and recalibrated). 

```bash
gatk Mutect2 -R ../Annotations/human_g1k_v37.fasta -I ../task2/Tumor.recal.bam \
-I ../task2/Control.recal.bam -normal TCGA-A7-A4SE-Control \
--germline-resource ../task3/filtered.vcf -L ../Captured_Regions.bed \
-O somatic_unfiltered_mt2.vcf
```

A needed tool to filter the raw callset found by Mutect2 is `FilterMutectCalls`:

```bash
gatk FilterMutectCalls -R ../Annotations/human_g1k_v37.fasta -V somatic_unfiltered_mt2.vcf -O filtered_mt2.vcf
```

Then use `SelectVariants` to select SNPs:

```bash
gatk SelectVariants -R ../Annotations/human_g1k_v37.fasta -V filtered_mt2.vcf \
-O snps_filtered_mt2.vcf --select-type-to-include SNP
```

Annotation is performed with snpEff (outside of the Docker container):

```bash
java -jar ~/snpEff/snpEff.jar -v hg19 tmp/snp_filtered_mt2_0806.vcf \
-s report_snpEff.html > ann_somatic_snps.vcf
```

The total number of point mutations can be found with:

```bash
cat snps_filtered_mt2.vcf | grep -v "^#" | wc -l
```

