## Determine which DNA repair genes overlap both heterozygous deletions and somatic point mutations of the patient

Somatic point mutations were called with Mutect in `task06_ALTERNATIVE.md`.

Starting from the filtered callsets of point mutations and other events, let's find the overlap between `DNA_Repair_Genes.bed` and the callset of somatic variants.

```bash
bedtools intersect -a ../DNA_Repair_Genes.bed -b snp_filtered.vcf | cut -f4,4 | uniq -c > snp_overlap.txt
```

```bash
bedtools intersect -a ../DNA_Repair_Genes.bed -b non_snp_filtered.vcf | cut -f4,4 | uniq -c > indel_overlap.txt
```

These .txt files were imported in R to intersect the gene found with the list of genes that harbor copy number deletions.
