First, we remove duplicates using Picard.

```bash
samtools sort Control.bam > Control.sorted.bam
samtools index Control.sorted.bam
java -jar ~/Documents/HumanGenomics/Tools/picard.jar MarkDuplicates \
I=Control.sorted.bam O=Control.sorted.dedup.bam REMOVE_DUPLICATES=true \
TMP_DIR=/tmp METRICS_FILE=Control.picard.log ASSUME_SORTED=true
samtools index Control.sorted.dedup.bam
```

```bash
samtools sort Tumor.bam > Tumor.sorted.bam
samtools index Tumor.sorted.bam
java -jar ~/Documents/HumanGenomics/Tools/picard.jar MarkDuplicates \
I=Tumor.sorted.bam O=Tumor.sorted.dedup.bam REMOVE_DUPLICATES=true \
TMP_DIR=/tmp METRICS_FILE=Control.picard.log ASSUME_SORTED=true
samtools index Tumor.sorted.dedup.bam

```


Then, we run RealignerTargetCreator and then IndelRealigner.

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

