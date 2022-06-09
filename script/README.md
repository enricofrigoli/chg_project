# Description

These markdown files contain the command used in the analysis, and some comment. 
As mentioned in the report, two setups were used:

* the first one (setup **1**) involves the use of a provided Virtual Machine that already included GATK (v3.8-1-0), Picard (v2.22.3), VarScan (v2.3.9), and snpEff (v4.3t)

* a second one (setup **2**), consisting of GATK (v4.2.6.1) run in a Docker container with Picard (v2.27.1), and SnpEff (5.1d).

Both setups share custom R scripts, the provided raw BAM files, few files containing the reference genome, the .bed files containing the captured regions and the .vcf files for variant annotation. When the second setup was used, the name of the markdown with the commands contains "ALTERNATIVE".
