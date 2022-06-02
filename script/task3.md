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

Then we need to filter the variants:

```bash
vcftools --minQ 20 --max-meanDP 200 --min-meanDP 5 --remove-indels \
--vcf Control.UniGen.vcf --out Control.UniGen --recode --recode-INFO-all
```

To find which identified SNPs are in the list clinvar_Pathogenic.vcf run:

```bash
vcftools --vcf Control.UniGen.vcf --diff ~/Documents/HumanGenomics/Annotations/clinvar_Pathogenic.vcf --diff-site
```
