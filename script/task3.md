# SNPs calling

Task 3 is about heterozygous germline SNPs calling. To do this, we have three options:

* Use `bcftools`
* Use `UnifiedGenotyper` from GATK
* Use `HaplotypeCaller` from GATK

It has been requested to use GATK, so let's start with UnifiedGenotyper:

```bash
java -jar ~/Documents/HumanGenomics/Tools/GenomeAnalysisTK.jar -T UnifiedGenotyper \
-R ~/Documents/HumanGenomics/Annotations/human_g1k_v37.fasta -I Control.sorted.dedup.realigned.recal.bam \
-o Control.UniGen.vcf -L ../Captured_Regions.bed
```

Then we need to filter the variants to keep only SNPs and remove indels:

```bash
vcftools --minQ 20 --max-meanDP 200 --min-meanDP 5 --remove-indels \
--vcf Control.UniGen.vcf --out Control.UniGen --recode --recode-INFO-all
```

# SNPs annotation

To further annotate the called variants, we use the tool SnpSift. In this way, we add annotation taken from the hapmap_3.3.b37.vcf file, and from clinvar_Pathogenic.vcf file. 

```bash
java -Xmx4g -jar ~/Documents/HumanGenomics/Tools/snpEff/SnpSift.jar \
Annotate ~/Documents/HumanGenomics/Annotations/hapmap_3.3.b37.vcf \
Control.UniGen.recode.vcf > Control.UniGen.recode.ann_hapmap.vcf
```

```bash
java -Xmx4g -jar ~/Documents/HumanGenomics/Tools/snpEff/SnpSift.jar \
Annotate ~/Documents/HumanGenomics/Annotations/clinvar_Pathogenic.vcf \
Control.UniGen.recode.ann_hapmap.vcf > Control.UniGen.recode.ann_clinv.vcf
```

Finally, let's count how many SNPs of our callset are in clinvar dataset

```bash
cat Control.UniGen.recode.ann_clinv.vcf ~/Documents/HumanGenomics/Tools/snpEff/SnpSift.jar filter ("exists CLNSIG")
```

