## Deduplication

First, we remove duplicates using Picard.

```bash
samtools sort Control.bam > Control.sorted.bam
samtools index Control.sorted.bam
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

To perform the realignment, we run RealignerTargetCreator and then IndelRealigner, specifying
the regions we want to target with `Captured_Regions.bed`.

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

To count ...? run:

```bash
samtools view Control.sorted.dedup.realigned.bam | grep OC | wc -l
```

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

To count ...? run:

```bash
samtools view Tumor.sorted.dedup.realigned.bam | grep OC | wc -l
```

## Recalibration

Firstly, we run:

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

We use AnalyzeCovariates to perform a before/after comparison and get the differences in a .cvs and in a .pdf file.

```bash
java -jar ~/Documents/HumanGenomics/Tools/GenomeAnalysisTK.jar \
-T AnalyzeCovariates -R ~/Documents/HumanGenomics/Annotations/human_g1k_v37.fasta \
-before recal.Control.table -after after_recal.Control.table \
-csv recal.Control.csv -plots recal.Control.pdf
```

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



### temporary notes

In the GATK docker container I did:

```bash
gatk MarkDuplicates -I Control.sorted.bam -O Control.dedup.bam -ASSUME_SORT_ORDER coordinate \
-REMOVE_DUPLICATES true -M Control_dedup_metrics.txt
```

and the same for tumor



