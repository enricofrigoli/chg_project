Task 10
================

Use *GATK ASEReadCounter* with the following setup, the same as task 9
bit with different sites. In this case we considere all SNP sites stored
in the hapmap file (`hapmap_3.3.b37.vcf`) and not only the ones found
with the variant calling.

It is better to filter out first the heterozygous SNPs to avoid warn
messages from ASEReadCounter.

``` bash
grep -E "(^#|0/1)" hapmap_3.3.b37.vcf > hapmap.het.vcf
```

``` bash
java -jar ../Tools/GenomeAnalysisTK.jar -T ASEReadCounter \
-R human_g1k_v37.fasta \
-o control.hap.csv \
-I Control.sorted.dedup.realigned.recal.bam \
-sites hapmap.het.vcf  \
-U ALLOW_N_CIGAR_READS \
-minDepth 20 \
--minMappingQuality 20 \
--minBaseQuality 20
```

Same for the tumor sample:

``` bash
java -jar ../Tools/GenomeAnalysisTK.jar -T ASEReadCounter \
-R human_g1k_v37.fasta \
-o tumor.hap.csv \
-I Tumor.sorted.dedup.realigned.recal.bam \
-sites hapmap.het.vcf  \
-U ALLOW_N_CIGAR_READS \
-minDepth 20 \
--minMappingQuality 20 \
--minBaseQuality 20
```

## SPIA analysis

\[SPIA documentation\]
(<https://cran.r-project.org/web/packages/SPIAssay/SPIAssay.pdf>)

``` r
library(data.table)
library(SPIAssay)
```

Load the tables created with *GATK ASEReadCounter* and use them to
compute the **allelic fractions**. The interesting columns of the
`control` and `tumor` tables are the `altCount` (number of occurrences
of alternative base) and `totalCount` (number of occurrences of the
SNP).

``` r
control <-  fread("control.hap.csv",data.table=F)
control_af <-  data.frame(control = control$altCount/control$totalCount)
rownames(control_af) <- control$variantID

tumor <- fread("tumor.hap.csv",data.table=F)
tumor_af <- data.frame(tumor = tumor$altCount/tumor$totalCount)
rownames(tumor_af) <- tumor$variantID
```

Encode the data into SPIA format based on allelic fractions: SPIA uses 0
for AA, 1 for BB, 2 for AB.

``` r
genotype <- merge(control_af, tumor_af, by=0)
genotype$control[-which(genotype$control %in% c(0,1))] <- 2
genotype$tumor[-which(genotype$tumor %in% c(0,1))] <- 2
colnames(genotype) <- c("variantID", "control", "tumor")
```

Set the SPIA parameters (still to be defined)

``` r
spia_parameters <- list(Pmm = 0.1, nsigma = 2, Pmm_nonM = 0.6, nsigma_nonM = 3, PercValidCall=0.7)
```

Perform the SPIA test

``` r
SPIA_commmon <-  SPIATest(x = genotype, row.names = TRUE, test.param = spia_parameters) 
```

    ## .

``` r
SPIAPlot(SPIA_commmon)
```

![](SPIA_task10_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
SPIA_commmon["SPIAresult"]
```

    ## $SPIAresult
    ##      Cell_1    Cell_2  Distance  SPIA_Score  SNP_available Total_SNP One_SNP_NA
    ## [1,] "control" "tumor" "0.18809" "Uncertain" "13196"       "13196"   "0"       
    ##      Bot_SNP_NA Diff_AvsB_or_BvsA Diff_AorBvsAB_or_vic DiffABvsAorB
    ## [1,] "0"        "0"               "2482"               "947"       
    ##      counterBothHomoz counterBothHeter
    ## [1,] "6744"           "3970"
