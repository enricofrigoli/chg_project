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

```bash
java -jar ~/Documents/HumanGenomics/Tools/VarScan.v2.3.9.jar somatic \
Control.pileup Tumor.pileup  --output-snp somatic.pm --output-vcf 1
```

```bash
java -Xmx8g -jar ~/Documents/HumanGenomics/Tools/snpEff/snpEff.jar \
-v hg19 somatic.pm.vcf -s somatic.pm.vcf.html > somatic.pm.ann.vcf
```

```bash
less somatic.pm.ann.vcf
```

Filter the vcf for somatic point mutations:

```bash
cat somatic.pm.ann.vcf | java -Xmx4g -jar ../../Tools/snpEff/SnpSift.jar filter "(SS = 2)"
```
Non so se questo funziona perché SS è un parametro di INFO nel .vfc

Alternativamente si possono contare con un grep per "SS=2"
