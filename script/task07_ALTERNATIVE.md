# Which DNA repair genes (use file DNA_Repair_Genes.bed) overlap both heterozygous deletions and heterozygous SNPs of the patient that are in Clinvar?

We already know that the only germline pathogenic SNPs is the one in BRCA1, which also in an hemideleted fragment (in the tumor sample). Let's check if other germline variants (that are not necessarily linked to a disease are in fragments that harbor somatic heterozygous deletions.

Let's start by intersecting the bed files containing the heterozygous (copy number) deletions with the bed containing DNA Repair Genes.

```bash
bedtools intersect -a ../DNA_Repair_Genes.bed \
-b /gatk/my_downloads/heterozygous.deletions.bed > hemideleted_dna_repair_genes.bed
```

In this way we retrieve the the subset of DNA Repair Genes that harbor heterozygous deletions.
Then we intersect this new bed file with the heterozygous germline SNPs callset:

```bash
bedtools intersect -a hemideleted_dna_repair_genes.bed -b het_filtered_clinv_annotated.vcf \
| cut -f4,4 | uniq -c > snps_hemidel_dna_repair.txt
```

This .txt is then used to generate the plot.
