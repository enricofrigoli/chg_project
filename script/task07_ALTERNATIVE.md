# Which DNA repair genes (use file DNA_Repair_Genes.bed) overlap both heterozygous deletions and heterozygous SNPs of the patient that are in Clinvar? [setup 2]

We already know that the only germline pathogenic SNPs is the one in BRCA1, which is also in an hemideleted fragment (in the tumor sample). Let's check if other germline variants (that are not necessarily linked to a disease are in fragments that harbor somatic heterozygous deletions.

Let's start by intersecting the bed files containing the heterozygous (copy number) deletions with the bed containing DNA Repair Genes.

```bash
bedtools intersect -a ../DNA_Repair_Genes.bed \
-b /gatk/my_downloads/heterozygous.deletions.bed > hemideleted_dna_repair_genes.bed
```

In this way we retrieve a new bed file containing the subset of DNA Repair Genes that harbor heterozygous deletions.
Then we intersect this new bed file with the heterozygous germline SNPs callset (from HaplotypeCaller):

```bash
bedtools intersect -a hemideleted_dna_repair_genes.bed -b het_filtered_clinv_annotated.vcf \
| cut -f4,4 | uniq -c > snps_hemidel_dna_repair.txt
```

This .txt is then used to generate the plot.

Note: here we first intersected the bed containing hemideleted regions with the bed of DNA repair genes, and then we found the variants that were inside these intervals; for the task8, the vcf of somatic variants was firstly intersected with the bed of DNA repair genes, and then, with R, the variants in DNA repair genes were differentiated based on the presence or not of a hemizygous deletion in that segment.
