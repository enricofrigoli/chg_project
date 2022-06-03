# Task 9 - Purity and Ploidy estimation

_GATK ASEReadCounter_ calculates read counts per allele for allele-specific expression analysis of RNAseq (??) data. The `.vcf` file is needed to specify the positions to evaluate: in this case, the variants identified during task 3 on the control sample.
Only biallelic heterozygous SNPs are considered, so it is quicker to extract them before with:

```bash
grep -E "(^#|0/1)" Control.UniGen.vcf > control.het.vcf
```

ASEReadCounter was run with the following setup:

```bash
java -jar ../Tools/GenomeAnalysisTK.jar -T ASEReadCounter \
-R human_g1k_v37.fasta \
-o control.csv \
-I Control.sorted.dedup.realigned.recal.bam \
-sites control.het.vcf  \
-U ALLOW_N_CIGAR_READS \
-minDepth 20 \
--minMappingQuality 20 \
--minBaseQuality 20
```

Parameters explanation:
* `-U ALLOW_N_CIGAR_READS`: allows the `N` string in CIGAR (it is required for RNASeq so I do not know why it is here, it is used in lesson 10 of the lab)
* `-MinDepth`: minimum number of bases that pass filters
* `--minMappingQuality`: minimum read mapping quality
* `--minBaseQuality`: minimum base quality

Same for the tumor sample:

```bash
java -jar ../Tools/GenomeAnalysisTK.jar -T ASEReadCounter \
-R human_g1k_v37.fasta \
-o tumor.csv \
-I Tumor.sorted.dedup.realigned.recal.bam \
-sites control.het.vcf  \
-U ALLOW_N_CIGAR_READS \
-minDepth 20 \
--minMappingQuality 20 \
--minBaseQuality 20
```

__Segmentation__ file from the somatic copy number calling is also needed

## CLONET: CLONality Estimate in Tumors

Rscript: I will upload it when it is ready

## TPES

Rscript: I will upload it when it is ready