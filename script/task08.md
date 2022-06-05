## Determine which DNA repair genes overlap both heterozygous deletions and somatic point mutations of the patient

Somatic point mutations were called with VarScan in `task06.md`, let's call also indels.

```bash
java -jar ~/Documents/HumanGenomics/Tools/VarScan.v2.3.9.jar \
somatic ../task6/Control.pileup ../task6/Tumor.pileup \
--output-indel somatic.indel --output-vcf 1
```

Let's find the overlap between `DNA_Repair_Genes.bed` and the callset of somatic variants.

```bash
bedtools intersect -a ../DNA_Repair_Genes.bed -b somatic.indel.vcf > indel_overlap.tsv
```

```bash
bedtools intersect -a ../DNA_Repair_Genes.bed -b ../task6/somatic.pm.vcf > pm_overlap.tsv
```

```bash
cat indel_overlap.tsv | cut -f4,4 | uniq > indel_genes.txt
```

```bash
cat pm_overlap.tsv | cut -f4,4 | uniq > pm_genes.txt
```

```bash
comm -12 indel_genes.txt pm_genes.txt > comm_overlap.txt
```


