# Task 10 - SPIA

Use _GATK ASEReadCounter_ with the following setup, the same as task 9 bit with different sites. In this case we considere all SNP sites stored in the hapmap file (`hapmap_3.3.b37.vcf`) and not only the ones found with the variant calling.

```bash
java -jar ../Tools/GenomeAnalysisTK.jar -T ASEReadCounter \
-R human_g1k_v37.fasta \
-o control.hap.csv \
-I Control.sorted.dedup.realigned.recal.bam \
-sites hapmap_3.3.b37.vcf  \
-U ALLOW_N_CIGAR_READS \
-minDepth 20 \
--minMappingQuality 20 \
--minBaseQuality 20
```

Same for the tumor sample:

```bash
java -jar ../Tools/GenomeAnalysisTK.jar -T ASEReadCounter \
-R human_g1k_v37.fasta \
-o tumor.hap.csv \
-I Tumor.sorted.dedup.realigned.recal.bam \
-sites hapmap_3.3.b37.vcf  \
-U ALLOW_N_CIGAR_READS \
-minDepth 20 \
--minMappingQuality 20 \
--minBaseQuality 20
```

## SPIA analysis

Rscript: I will upload it when it is ready