We propose here two paraller pipelines to process the data: the one we saw in class (using UnifiedGenotyper etc.) that makes use of tools that are no longer considered "best-practice", and another one that uses the most updated tools (HaplotypeCaller etc.). The following is the description of the second pipeline. 

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
gatk BaseRecalibrator -I Control.dedup.bam \
-R ../Annotations/human_g1k_v37.fasta \
--known-sites ../Annotations/hapmap_3.3.b37.vcf \
-O recal.Control.table -L ../Captured_Regions.bed
```

```bash
gatk ApplyBQSR -R ../Annotations/human_g1k_v37.fasta \
-I Control.dedup.bam --bqsr-recal-file recal.Control.table \
-O Contro.recal.bam -L ../Captured_Regions.bed \
--emit-original-quals true
```

Then, we need to perform a second step of recalibration on the recalibrated file for quality control purposes.

```bash
gatk BaseRecalibrator -R ../Annotations/human_g1k_v37.fasta -I Control.recal.bam \
-O after.recal.Control.table -L ../Captured_Regions.bed --known-sites ../Annotations/hapmap_3.3.b37.vcf
```

Finally, AnalyzeCovariates allows to generate before/after plots.

```bash
gatk AnalyzeCovariates -before recal.Control.table -after after.recal.Control.table \
-csv recal_Control_report.csv -plots recal_Control_plots.pdf
```


## Variant Calling

To perform variant calling, we used HaplotypeCaller module from GATK.

```bash
gatk HaplotypeCaller -R ../Annotations/human_g1k_v37.fasta -I ../task2/Control.recal.bam \
-O variants_HC.vcf -L ../Captured_Regions.bed
```

Following the GATK recommendation, the variant callset is refined using CNNScoreVariants, that annotates the VCF with scores from a CNN using a 2D model with pre-trained architecture. This tool annotates each variant with a score indicating the model's prediction of the quality of each variant.

```bash
gatk CNNScoreVariants -I ../task2/Control.recal.bam -V variants_HC.vcf \
-R ../Annotations/human_g1k_v37.fasta -O annotated.vcf -tensor-type read_tensor
```

Since we are interested only in SNPs and not indels, let's remove them with SelectVariants.

```bash
gatk SelectVariants -R ../Annotations/human_g1k_v37.fasta -V annotated.vcf \
--select-type-to-include SNP -O annotated_snp.vcf
```

The callset is further refined using FilterVariantTranches, that filter SNP calls using the scores computed in the previous step and the provided hapmap known variants to further refine the callset. Default SNP sensitivity is set at 99.95

```bash
gatk FilterVariantTranches -V annotated_snp.vcf --resource ../Annotations/hapmap_3.3.b37.vcf \
 --info-key CNN_2D -O filtered.vcf
```

Let's compute the number of heterozygous SNPs found by GATK.

```bash
cat filtered.vcf | grep "0/1" | wc -l
```

6937 heterozygous SNPs were found.


## Variant Annotation

Variant annotation was perfomed with SnpSift from snpEff, using clinvar_Pathogenic.vcf

```bash
sudo java -jar ~/snpEff/SnpSift.jar Annotate Annotations/clinvar_Pathogenic.vcf filtered.vcf > annotated_clinvar_snpsift.vcf
```

Using the same tool we can find the variants that has a clinical significance retrieved from clinvar_Pathogenic.vcf 

```bash
cat annotated_clinvar_snpsift.vcf | java -jar ~/snpEff/SnpSift.jar filter "(exists CLNSIG)"
```

As we can see, there is a SNP in the gene BRCA1 that is associated with hereditary breast and ovarian cancer syndrome.

