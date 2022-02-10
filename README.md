# Genome Analysis Toolkit (GATK)
The GATK is the industry standard for identifying SNPs and indels in germline DNA and RNAseq data. Its scope is now expanding to include somatic short variant calling, and to tackle copy number (CNV) and structural variation (SV). In addition to the variant callers themselves, the GATK also includes many utilities to perform related tasks such as processing and quality control of high-throughput sequencing data, and bundles the popular Picard toolkit.

It is a collection of command-line tools for analyzing high-throughput seqiencing data with a primary focus on variant discovery.Tools can be used individually or chained together into complete workflows.

This page contains the GATK best practice pipelines.

## Reference data pre-processing
Download Hg38 Reference Genome and unzip it
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```
BWA requires an idenxed reference genome. This step will produce 5 files (amb, ann, bwt, pac and sa) needed for the alignment of the reads to the genome.
-a (algorithm) use bwtsw which works for most genomes unless they are under a certain size.
```
bwa index -a bwtsw hg38.fa or bwa index hg38.fa
```
Create sequence dictionary using the Picard tool CreateSequenceDictionary.jar.
The tool will create and name the output as hg38.dict.
```
java -jar CreateSequenceDictionary.jar \
REFERENCE=hg38.fa
```
Index the reference genome using samtools
```
samtools faidx hg38.fa
```

Order of process
1. [Map and clean up short read sequence data effectively](https://github.com/marianakhoul/GATK/tree/main/Map%20and%20clean%20up%20short%20read%20sequence%20data%20efficiently)
2. [Data pre-processing for variant discovery](https://github.com/marianakhoul/GATK/tree/main/Data%20pre-processing%20for%20variant%20discovery)

Downloading [https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4]

## Optimizations for code --java-options
Java heap size
```
 -Xms - set initial Java heap size
 -Xmx - set maximum Java heap size
 -Xss - set java thread stack size
```
-Xms4000m means the initial heap size will be 4Gb and this helps to not readjust later on.



## References
[https://gatk.broadinstitute.org/hc/en-us/categories/360002302312]
