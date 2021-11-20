# (How to) Map and clean up short read sequence data efficiently

<img width="790" alt="Screen Shot 2021-11-20 at 5 55 45 PM" src="https://user-images.githubusercontent.com/31465978/142743126-ce947ded-7025-4289-a007-0fec21f27ff3.png">

## Step 1. Generate an unmapped BAM from FASTQ or aligned BAM

If you have raw read data in BAM format with appropriate assigned read group fileds, you can start with Step 2. Read group ID should differentiate factors contributing to technical batch effects. If no read group, will need to reassign the read group fields. 

### Check for read group information

```
samtools view -H sample.bam | grep '^@RG'
```
These tags, when assigned appropriately, allow to differentiate not only samples, but also various technical features that are associated with artifacts. With this information in hand, it can mitigate the effects of those artifacts during the duplicate marking and base recalibration steps. The GATK requires several read group fields to be present in input files and will fail with errors if this requirement is not satisfied. More on the fields can be found [here][1]


### Convert FASTQ to uBAM and add read group information using FastqToSAM

Picard's FastqToSam transforms a FASTQ file to an unmapped BAM. Requires 2 read group fields and makes optional specification of other read group fields. Below command only has the GATK best practice workflow parameters.
```
java -Xmx8G -jar picard.jar FastqToSam \
    FASTQ=6484_snippet_1.fastq \ #first read file of pair
    FASTQ2=6484_snippet_2.fastq \ #second read file of pair
    OUTPUT=6484_snippet_fastqtosam.bam \
    READ_GROUP_NAME=H0164.2 \ #required; changed from default of A
    SAMPLE_NAME=NA12878 \ #required
    LIBRARY_NAME=Solexa-272222 \ #required 
    PLATFORM_UNIT=H0164ALXX140820.2 \ 
    PLATFORM=illumina \ #recommended
    SEQUENCING_CENTER=BI \ 
    RUN_DATE=2014-08-20T00:00:00-0400
```
For paired reads, specify each FASTQ file with FASTQ and FASTQ2. Tool assumes identical ordering for pairs -- records in file must be queryname sorted. For single end reads, just specify FASTQ.
QUALITY_FORMAT is detected automatically if unspecified.
SORT_ORDER by default is queryname.

### Convert aligned BAM to uBAM and discard problematic records using RevertSam

Removes alignment information and generate an unmapped BAM.
```
java -Xmx8G -jar /path/picard.jar RevertSam \
    I=6484_snippet.bam \
    O=6484_snippet_revertsam.bam \
    SANITIZE=true \ 
    MAX_DISCARD_FRACTION=0.005 \ #informational; does not affect processing
    ATTRIBUTE_TO_CLEAR=XT \
    ATTRIBUTE_TO_CLEAR=XN \
    ATTRIBUTE_TO_CLEAR=AS \ #Picard release of 9/2015 clears AS by default
    ATTRIBUTE_TO_CLEAR=OC \
    ATTRIBUTE_TO_CLEAR=OP \
    SORT_ORDER=queryname \ #default
    RESTORE_ORIGINAL_QUALITIES=true \ #default
    REMOVE_DUPLICATE_INFORMATION=true \ #default
    REMOVE_ALIGNMENT_INFORMATION=true #default
```
Use the below command to see the tags in a BAM to be removed:
```
samtools view input.bam | cut -f 12- | tr '\t' '\n' | cut -d ':' -f 1 | awk '
```

## Step 2. Mark adapter sequences using MarkIlluminaAdapters

MarkIlluminaAdapters adds the XT tag to a record to mark the 5'start position of the adapter sequence and produces a metrics file. Tools use the XT tag to effectively remove adapter sequence contributions to read alignment and alignment scoring metrics.
```
java -Xmx8G -jar /path/picard.jar MarkIlluminaAdapters \
I=6484_snippet_revertsam.bam \
O=6484_snippet_markilluminaadapters.bam \
M=6484_snippet_markilluminaadapters_metrics.txt \ #naming required
TMP_DIR=/path/shlee #optional to process large files
```
BAM files are identical, but the adapters are marked in the new BAM.




### Reference
1. https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently

2. https://sites.google.com/a/broadinstitute.org/legacy-gatk-forum-discussions/tutorials/6484-how-to-generate-an-unmapped-bam-from-fastq-or-aligned-bam


[1]: https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
