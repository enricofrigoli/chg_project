## Data cleaning

Starting from `Control.bam`, we need to sort and index it.

```bash
samtools sort Control.bam > Control.sorted.bam
samtools index Control.sorted.bam
```

### Deduplication

Deduplication is performed with MarkDuplicates (from GATK) before Variant Calling:

```bash
gatk MarkDuplicates -I Control.sorted.bam -O Control.dedup.bam -ASSUME_SORT_ORDER coordinate \
-REMOVE_DUPLICATES true -M Control_dedup_metrics.txt
```

Using HaplotypeCaller does not require realignement, so we can proceed directly with recalibration.

### Recalibration

We use BaseRecalibrator tool from GATK that builds a model of covariation based on the input data and a set of known variants, producing a recalibration file; then we use ApplyBSQR (replacing PrintReads) tool to adjust the base quality scores in the data based on the model. The known variants are used to mask out bases at sites of real (expected) variation, to avoid counting real variants as errors. Outside of the masked sites, every mismatch is counted as an error.

An additional step involves the building of a second model and the generation of before/after plots to visualize the effects of the recalibration process.

```bash
gatk BaseRecalibrator -I Control.sorted.dedup.bam \
-R ~/Documents/Project_2022/Annotations/human_g1k_v37.fasta \
--known-sites ~/Documents/Project_2022/Annotations/hapmap_3.3.b37.vcf \
-O recal.Control.table -L ../Captured_Regions.bed
```

```bash
gatk ApplyBQSR -R ~/Documents/Project_2022/Annotations/human_g1k_v37.fasta \
-I Control.sorted.dedup.bam --bqsr-recal-file recal.Control.table \
-O Contro.recal.bam -L ../Captured_Regions.bed \
--emit-original-quals true
```


## Variant Calling

To perform variant calling, we used HaplotypeCaller module from GATK.

```bash
gatk HaplotypeCaller -R ../Annotations/human_g1k_v37.fasta -I ../task2/Control.dedup.bam \
-O variants_HC.vcf -L Captured_Regions.bed
```

## Variant Annotation

Variant annotation was perfomed with VariantAnnotator from GATK, using clinvar_Pathogenic.vcf

```bash
gatk VariantAnnotator -R ../Annotations/human_g1k_v37.fasta -V variants_HC.vcf \
-O annotated_var_VA.vcf --resource ../Annotations/clinvar_Pathogenic.vcf
```

We then filtered for the identified SNPs that are in the clinvar_Pathogenic.vcf using  snpEff

```bash
cat annotated_var_VA.vcf | java -jar ~/snpEff/SnpSift.jar filter "(exists CLNSIG)"
```


