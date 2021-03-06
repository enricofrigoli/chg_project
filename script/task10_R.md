Task 10 - Similarity between Control and Tumor samples
================

Use *GATK ASEReadCounter* with the following setup, the same as task 9
but with different sites. In this case we consider all SNP sites stored
in the hap map file (`hapmap_3.3.b37.vcf`) and not only the ones found
with the variant calling.

It is better to filter out first the heterozygous SNPs to avoid warn
messages from ASEReadCounter.

``` bash
grep -E "(^#|0/1)" hapmap_3.3.b37.vcf > ./task10/hapmap.het.vcf
```

``` bash
java -jar ../Tools/GenomeAnalysisTK.jar -T ASEReadCounter \
-R human_g1k_v37.fasta \
-o ./task10/control.hap.csv \
-I Control.sorted.dedup.realigned.recal.bam \
-sites ./task10/hapmap.het.vcf  \
-U ALLOW_N_CIGAR_READS \
-minDepth 20 \
--minMappingQuality 20 \
--minBaseQuality 20
```

Same for the tumor sample:

``` bash
java -jar ../Tools/GenomeAnalysisTK.jar -T ASEReadCounter \
-R human_g1k_v37.fasta \
-o ./task10/tumor.hap.csv \
-I Tumor.sorted.dedup.realigned.recal.bam \
-sites ./task10/hapmap.het.vcf  \
-U ALLOW_N_CIGAR_READS \
-minDepth 20 \
--minMappingQuality 20 \
--minBaseQuality 20
```

## SPIA analysis

[SPIA
documentation](https://cran.r-project.org/web/packages/SPIAssay/SPIAssay.pdf)

``` r
library(data.table)
library(SPIAssay)
library(ggplot2)
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
for AA, 1 for BB, 2 for AB. Only the SNPs present in both samples are
considered.

``` r
genotype <- merge(control_af, tumor_af, by=0)
genotype$control[-which(genotype$control %in% c(0,1))] <- 2
genotype$tumor[-which(genotype$tumor %in% c(0,1))] <- 2
colnames(genotype) <- c("variantID", "control", "tumor")
```

Set the SPIA parameters: they determine how the distance score is
considered. It can be different, similar or uncertain.

``` r
spia_parameters <- list(Pmm = 0.1, nsigma = 2, Pmm_nonM = 0.6, nsigma_nonM = 3, PercValidCall=1)
```

Perform the SPIA test

``` r
SPIA_commmon <-  SPIATest(x = genotype, row.names = TRUE, test.param = spia_parameters) 
```

    ## .

``` r
result_common <- t(as.data.frame(SPIA_commmon["SPIAresult"]))
result_common
```

    ##                                 [,1]       
    ## SPIAresult.Cell_1               "control"  
    ## SPIAresult.Cell_2               "tumor"    
    ## SPIAresult.Distance             "0.18809"  
    ## SPIAresult.SPIA_Score           "Uncertain"
    ## SPIAresult.SNP_available        "13196"    
    ## SPIAresult.Total_SNP            "13196"    
    ## SPIAresult.One_SNP_NA           "0"        
    ## SPIAresult.Bot_SNP_NA           "0"        
    ## SPIAresult.Diff_AvsB_or_BvsA    "0"        
    ## SPIAresult.Diff_AorBvsAB_or_vic "2482"     
    ## SPIAresult.DiffABvsAorB         "947"      
    ## SPIAresult.counterBothHomoz     "6744"     
    ## SPIAresult.counterBothHeter     "3970"

``` r
res2plot <- data.frame(counts = as.numeric(result_common[c(9,10,12,13),]))
rownames(res2plot) <- c("AA<->BB", "AA/BB<->AB", "Hom-Hom", "Het-Het")

ggplot(data=res2plot, aes(x=rownames(res2plot), y=counts)) +
  geom_bar(stat="identity", width = 0.5, fill="#ffb703")+
  xlab(NULL)
```

![](task10_R_files/figure-gfm/spia1_plot-1.png)<!-- -->

The results give a SPIA score of **0.18809**, marked as uncertain for
the parameters given. Since we only considered SNPs in common between
the samples, no NAs were present. The analysis also detected that there
were no changes between two homozygous variants (AA\<-\>BB) but only
changes from homozygous to heterozygous genotype and the opposite
(AA/BB\<-\>AB)(2482)
