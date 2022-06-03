# Somatic point mutations calling

To call somatic point mutations, we use VarScan.

```bash
samtools mpileup -q 1 -f ~/Documents/HumanGenomics/Annotations/human_g1k_v37.fasta \
../task2/Control.sorted.dedup.realigned.recal.bam > Control.pileup
```

```bash
samtools mpileup -q 1 -f ~/Documents/HumanGenomics/Annotations/human_g1k_v37.fasta \
../task2/Tumor.sorted.dedup.realigned.recal.bam > Tumor.pileup
```
