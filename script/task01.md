# Task 1

The statistics of the raw aligned reads contained in the two BAM files can be retrieved using `samtools stats` or `samtools flagstat`, for detailed and general statistics respectively. 

```bash
samtools stats Control.bam > stats_control.txt
samtools stats Tumor.bam > stats_tumor.txt
```

```bash
samtools flagstat Control.bam
samtools flagstat Tumor.bam
```
