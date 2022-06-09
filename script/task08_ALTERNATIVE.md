## Determine which DNA repair genes overlap both heterozygous deletions and somatic point mutations of the patient [setup 2].

Somatic point mutations were called with Mutect2 in `task06_ALTERNATIVE.md`. The callset of somatic point found with VarScan mutations was discarded since only 2 of them were found DNA Repair Genes and noone of them were part of the subset of hemideleted ones.

Starting from the filtered callsets of point mutations, let's find the overlap between `DNA_Repair_Genes.bed` and the callset of somatic variants.

```bash
bedtools intersect -a ../DNA_Repair_Genes.bed -b snp_filtered.vcf | cut -f4,4 | uniq -c > snp_overlap.txt
```

```bash
bedtools intersect -a ../DNA_Repair_Genes.bed -b non_snp_filtered.vcf | cut -f4,4 | uniq -c > indel_overlap.txt
```

These .txt files were imported in R to intersect the gene found with the list of genes that harbor copy number deletions.
