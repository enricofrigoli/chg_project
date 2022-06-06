## Which DNA repair genes (use file DNA_Repair_Genes.bed) overlap both heterozygous deletions and heterozygous SNPs of the patient that are in Clinvar?

Let's start with the import and rename of the raw callset (at the beginning of task03 we selected only SNPs, removing indels).

```bash
cp ../task3/Control.UniGen.vcf ./raw_callset.vcf
```

Then select for heterozygous variants.

```bash
grep -E "(^#|0/1)" raw_callset.vcf > het_variants.vcf
```

Let's filter using vcftools using the same filters used in task03.

```bash
vcftools --minQ 20 --max-meanDP 200 --min-meanDP 5 \
--vcf het_variants.vcf --out het_variants --recode --recode-INFO-all
```

Since annotation was performed only on SNPs, let's annotate the new file (comprising indels) starting with `snpEff`.

```bash
java -jar ~/Documents/HumanGenomics/Tools/snpEff/snpEff.jar -v hg19kg \
het_variants.recode.vcf -s het_variants.ann_kg.html > het_variants.ann_kg.vcf
```

Then use `snpSift` to annotate variants using as resources the hapmap and clinvar files.

```bash
java -jar ~/Documents/HumanGenomics/Tools/snpEff/SnpSift.jar \
Annotate ~/Documents/HumanGenomics/Annotations/hapmap_3.3.b37.vcf \
het_variants.ann_kg.vcf > het_variants.ann_kg.hapmap.vcf
```

```bash
java -jar ~/Documents/HumanGenomics/Tools/snpEff/SnpSift.jar \
Annotate ~/Documents/HumanGenomics/Annotations/clinvar_Pathogenic.vcf \
het_variants.ann_kg.hapmap.vcf > het_variants.ann_kg.hapmap.clvar.vcf
```

Finally, let's intersect the clinvar file with our annotation.

```bash
cat het_variants.ann_kg.hapmap.clvar.vcf | java -jar ~/Documents/HumanGenomics/Tools/snpEff/SnpSift.jar \
filter "(exists CLNSIG)" > intersect_clinvar.vcf
```

Now let's compute the overlap between germline heterozygous indels and SNPs with `DNA_Repair_Genes` using bedtools.

```bash
bedtools intersect -a ../DNA_Repair_Genes.bed -b intersect_clinvar.vcf
```

We can find the the only variant linked to a clinical significance is the one in BRCA1.
If we instead extend the intersection at all the heterozygous variants called, we found many SNPs in DNA repair genes.

```bash
bedtools intersect -a ../DNA_Repair_Genes.bed -b het_variants.ann_kg.hapmap.clvar.vcf > overlap.tsv
```

```bash
cat overlap.tsv | cut -f4,4 | uniq -c 
```



