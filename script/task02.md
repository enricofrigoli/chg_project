## Deduplication

Deduplication allows to remove technical duplicates among the pool of raw aligned reads. First we sort and indec the BAMs.

```bash
samtools sort Control.bam > Control.sorted.bam
samtools index Control.sorted.bam
```

Then, we remove duplicates using `MarkDuplicates` from Picard, specifying the input BAM, the desired name of the output, the option to remove duplicated instead of just marking them, the folder used as *temp*, the request for a file containing metrics, and the sorted status as true.

```bash
java -jar ~/Documents/HumanGenomics/Tools/picard.jar MarkDuplicates \
I=Control.sorted.bam O=Control.sorted.dedup.bam REMOVE_DUPLICATES=true \
TMP_DIR=/tmp METRICS_FILE=Control.picard.log ASSUME_SORTED=true
samtools index Control.sorted.dedup.bam
```

We do the same for the Tumor sample.

```bash
samtools sort Tumor.bam > Tumor.sorted.bam
samtools index Tumor.sorted.bam
java -jar ~/Documents/HumanGenomics/Tools/picard.jar MarkDuplicates \
I=Tumor.sorted.bam O=Tumor.sorted.dedup.bam REMOVE_DUPLICATES=true \
TMP_DIR=/tmp METRICS_FILE=Control.picard.log ASSUME_SORTED=true
samtools index Tumor.sorted.dedup.bam
```


## Realignment

Realignment is required for variant calling using `UnifiedGenotyper`.
To perform the realignment, we run `RealignerTargetCreator` and then `IndelRealigner`, specifying
the regions we want to target with `Captured_Regions.bed` and the reference genome.

```bash
java -jar ~/Documents/HumanGenomics/Tools/GenomeAnalysisTK.jar \
-T RealignerTargetCreator -R ~/Documents/HumanGenomics/Annotations/human_g1k_v37.fasta \
-I Control.sorted.dedup.bam -o realigner.intervals -L ../Captured_Regions.bed
```

```bash
java -jar ~/Documents/HumanGenomics/Tools/GenomeAnalysisTK.jar \
-T IndelRealigner -R ~/Documents/HumanGenomics/Annotations/human_g1k_v37.fasta \
-I Control.sorted.dedup.bam -targetIntervals realigner.intervals \
-o Control.sorted.dedup.realigned.bam -L ../Captured_Regions.bed
```

To count the number of realigned reads run:

```bash
samtools view Control.sorted.dedup.realigned.bam | grep OC | wc -l
```

where "OC" stands for *Original CIGAR*.

The same was done for the the Tumor.bam file:

```bash
java -jar ~/Documents/HumanGenomics/Tools/GenomeAnalysisTK.jar \
-T RealignerTargetCreator -R ~/Documents/HumanGenomics/Annotations/human_g1k_v37.fasta \
-I Tumor.sorted.dedup.bam -o realigner.Tumor.intervals -L ../Captured_Regions.bed
```

```bash
java -jar ~/Documents/HumanGenomics/Tools/GenomeAnalysisTK.jar \
-T IndelRealigner -R ~/Documents/HumanGenomics/Annotations/human_g1k_v37.fasta \
-I Tumor.sorted.dedup.bam -targetIntervals realigner.Tumor.intervals \
-o Tumor.sorted.dedup.realigned.bam -L ../Captured_Regions.bed
```

```bash
samtools view Tumor.sorted.dedup.realigned.bam | grep OC | wc -l
```


## Recalibration

Recalibration is instead required for all variant calling tools, that must use analysis-ready reads. Here (setup 1), it was perfomed using `BaseRecalibrator` that makes a recalibration table, which is then used by `PrintReads` to actually recalibrate the reads. Another round of recalibration is performed on the recalibrated file as a quality control.

```bash
java -jar ~/Documents/HumanGenomics/Tools/GenomeAnalysisTK.jar -T BaseRecalibrator \
-R ~/Documents/HumanGenomics/Annotations/human_g1k_v37.fasta -I Control.sorted.dedup.realigned.bam \
-knownSites ~/Documents/HumanGenomics/Annotations/hapmap_3.3.b37.vcf \
-o recal.Control.table -L ../Captured_Regions.bed
```

```bash
java -jar ~/Documents/HumanGenomics/Tools/GenomeAnalysisTK.jar -T PrintReads \
-R ~/Documents/HumanGenomics/Annotations/human_g1k_v37.fasta \
-I Control.sorted.dedup.realigned.bam -BQSR recal.Control.table \
-o Control.sorted.dedup.realigned.recal.bam -L ../Captured_Regions.bed \
--emit_original_quals
```

```bash
java -jar ~/Documents/HumanGenomics/Tools/GenomeAnalysisTK.jar \
-T BaseRecalibrator -R ~/Documents/HumanGenomics/Annotations/human_g1k_v37.fasta \
-I Control.sorted.dedup.realigned.bam \
-knownSites ~/Documents/HumanGenomics/Annotations/hapmap_3.3.b37.vcf \
-BQSR recal.Control.table -o after_recal.Control.table -L ../Captured_Regions.bed
```

The recalibration tables retrieved from the first and the recond round of recalibration were used to compute recalibration plots with `AnalyzeCovariates`.

```bash
java -jar ~/Documents/HumanGenomics/Tools/GenomeAnalysisTK.jar \
-T AnalyzeCovariates -R ~/Documents/HumanGenomics/Annotations/human_g1k_v37.fasta \
-before recal.Control.table -after after_recal.Control.table \
-csv recal.Control.csv -plots recal.Control.pdf
```

The number of recalibrated reads can be computed grepping "OQ" that stands for *Original Quality*, that was emitted for every recalibrated reads using the argument `--emit_original_quals` in the previous command.

```bash
samtools view Control.sorted.dedup.realigned.recal.bam | grep OQ | wc -l
```


And then we do the same for the Tumor sample.

```bash
java -jar ~/Documents/HumanGenomics/Tools/GenomeAnalysisTK.jar -T BaseRecalibrator \
-R ~/Documents/HumanGenomics/Annotations/human_g1k_v37.fasta -I Tumor.sorted.dedup.realigned.bam \
-knownSites ~/Documents/HumanGenomics/Annotations/hapmap_3.3.b37.vcf \
-o recal.Tumor.table -L ../Captured_Regions.bed
```

```bash
java -jar ~/Documents/HumanGenomics/Tools/GenomeAnalysisTK.jar -T PrintReads \
-R ~/Documents/HumanGenomics/Annotations/human_g1k_v37.fasta \
-I Tumor.sorted.dedup.realigned.bam -BQSR recal.Tumor.table \
-o Tumor.sorted.dedup.realigned.recal.bam -L ../Captured_Regions.bed \
--emit_original_quals
```

```bash
java -jar ~/Documents/HumanGenomics/Tools/GenomeAnalysisTK.jar \
-T BaseRecalibrator -R ~/Documents/HumanGenomics/Annotations/human_g1k_v37.fasta \
-I Tumor.sorted.dedup.realigned.bam \
-knownSites ~/Documents/HumanGenomics/Annotations/hapmap_3.3.b37.vcf \
-BQSR recal.Tumor.table -o after_recal.Tumor.table -L ../Captured_Regions.bed
```

```bash
java -jar ~/Documents/HumanGenomics/Tools/GenomeAnalysisTK.jar \
-T AnalyzeCovariates -R ~/Documents/HumanGenomics/Annotations/human_g1k_v37.fasta \
-before recal.Tumor.table -after after_recal.Tumor.table \
-csv recal.Tumor.csv -plots recal.Tumor.pdf
```

```bash
samtools view Tumor.sorted.dedup.realigned.recal.bam | grep OQ | wc -l
```


```bash
gatk HaplotypeCaller -R ../Annotations/human_g1k_v37.fasta -I ../task2/Control.dedup.bam -O variants_HC.vcf -L Captured_Regions.bed
```



