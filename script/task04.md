# Task 4

Run Ancestry Analysis with the EthSEQ R package and the following setup:

```R
library(EthSEQ)

ethseq.Analysis(
  bam.list = "./BAMs_list.txt",
  out.dir = "./",
  model.gds = "./SS2.Light.Model.gds",
  cores=1,
  mbq = 20,
  mrq = 1,
  mdc = 10,
  run.genotype = TRUE,
  verbose=TRUE,
  composite.model.call.rate = 1,
  space="3D")
```
input BAM file: `Control.sorted.dedup.realigned.recal.bam`

input model: `SS2.Light.Model.gds`--> reference model

Other parameters:
* `cores`: number of parallel cores to be used
* `mbq`: minimum base quality used for the pileup
* `mrq`: minimum read quality used in the piluep 
* `mdc`: minimum read count accettable for genotype inference
* `run.genotype`: to run the genotype (from ASEQ package)
* `composite.model.call.rate`: SNP call rate for PCA
* `space`: number od principal component dimensions to be used


The BAM file is genotyped at _reference model_ positions and a _target model_ is created. 
PCA is performed on the aggregates target and reference models. Only SNPs that satisfy the requirements set as input are retained.
As output, some plots are generated with the first three principal components. 
The reference are depicted as polygons, the sample is a point falling into the Africa region.

A multi-step refinement analysis could be possibly performed but it is not clear how and if there is enough data in the reference model (GDS)

More information:
* [documentation](https://cran.r-project.org/web/packages/EthSEQ/EthSEQ.pdf)
* [vignette](https://cran.r-project.org/web/packages/EthSEQ/vignettes/EthSEQ.html)
* [paper](https://academic.oup.com/bioinformatics/article/33/15/2402/3091083)
