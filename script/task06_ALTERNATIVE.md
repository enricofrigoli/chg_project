Starting from the filtered callset from Mutect2 (somatic short variants), we used `VariantAnnotator` from GATK to annotate the callset with `clinvar_Pathogenic.vcf` file.

```bash
gatk VariantAnnotator -R ../Annotations/human_g1k_v37.fasta \
--resource ../Annotations/clinvar_Pathogenic.vcf -V filtered_mt2.vcf  \
-O filtered_ann.vcf
```





