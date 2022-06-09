# Germline SNPs calling

To perform germline SNPs calling with the setup 1, `UnifiedGenotyper` was used.

```bash
java -jar ~/Documents/HumanGenomics/Tools/GenomeAnalysisTK.jar -T UnifiedGenotyper \
-R ~/Documents/HumanGenomics/Annotations/human_g1k_v37.fasta -I Control.sorted.dedup.realigned.recal.bam \
-o Control.UniGen.vcf -L ../Captured_Regions.bed
```

Then, variants were filtered to keep only SNPs and remove indels:

```bash
vcftools --minQ 20 --max-meanDP 200 --min-meanDP 5 --remove-indels \
--vcf Control.UniGen.vcf --out Control.UniGen --recode --recode-INFO-all
```

How many heterozygous SNPs have been called?

```bash
cat Control.UniGen.recode.vcf | grep "0/1" | wc -l
```


# SNPs annotation

To further annotate the called variants, we use the tool snpEff. In this way, we add annotation taken from the internal hg19kg (known genes) database, hapmap_3.3.b37.vcf file, and from clinvar_Pathogenic.vcf file. 

```bash
java -jar ~/Documents/HumanGenomics/Tools/snpEff/snpEff.jar -v hg19kg \
Control.UniGen.recode.vcf -s Control.ann_kg.html > Control.ann_kg.vcf
```

```bash
java -jar ~/Documents/HumanGenomics/Tools/snpEff/SnpSift.jar \
Annotate ~/Documents/HumanGenomics/Annotations/hapmap_3.3.b37.vcf \
Control.ann_kg.vcf > Control.ann_kg.hapmap.vcf
```

```bash
java -jar ~/Documents/HumanGenomics/Tools/snpEff/SnpSift.jar \
Annotate ~/Documents/HumanGenomics/Annotations/clinvar_Pathogenic.vcf \
Control.ann_kg.hapmap.vcf > Control.ann_kg.hapmap.clvar.vcf
```

Select heterozygous ones:

```bash
cat Control.ann_kg.hapmap.clvar.vcf | grep -E "(^#|0/1)" > het_Control_ann_kg.hapmap.clvar.vcf
```

Which identified SNPs are in the list `clinvar_Pathogenic.vcf`?

```bash 
bedtools intersect -a Control.ann_kg.hapmap.vcf -b ~/Documents/HumanGenomics/Annotations/clinvar_Pathogenic.vcf | less
```

Which SNPs have clinical significance among the annotated ones?

```bash
cat Control.ann_kg.hapmap.clvar.vcf | java -jar ~/Documents/HumanGenomics/Tools/snpEff/SnpSift.jar \
filter "(exists CLNSIG)"
```

