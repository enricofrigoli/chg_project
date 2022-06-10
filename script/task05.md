# Somatic Copy Number Variant Calling

To call CNV we will use VarScan (v2.3.9). First we use the module copynumber on the output of mpileup.

```bash
samtools mpileup -q 1 -f ~/Documents/HumanGenomics/Annotations/human_g1k_v37.fasta \
../task2/Control.sorted.dedup.realigned.recal.bam ../task2/Tumor.sorted.dedup.realigned.recal.bam \
| java -jar ~/Documents/HumanGenomics/Tools/VarScan.v2.3.9.jar \
copynumber --output-file SCNA --mpileup 1
```

Then we use the module copyCaller that generates the file that will be analyzed with DNAcopy in R.

```bash
java -jar ~/Documents/HumanGenomics/Tools/VarScan.v2.3.9.jar copyCaller SCNA.copynumber \
--output-file SCNA.copynumber.called
```

The analysis continues on the R markdown file `task05_R.md`
