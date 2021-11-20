# (How to) Map and clean up short read sequence data efficiently

<img width="790" alt="Screen Shot 2021-11-20 at 5 55 45 PM" src="https://user-images.githubusercontent.com/31465978/142743126-ce947ded-7025-4289-a007-0fec21f27ff3.png">

## Step 1. Generate an unmapped BAM from FASTQ or aligned BAM

If you have raw read data in BAM format with appropriate assigned read group fileds, you can start with Step 2. Read group ID should differentiate factors contributing to technical batch effects. If no read group, will need to reassign the read group fields. 

### Check for read group information

These tags, when assigned appropriately, allow to differentiate not only samples, but also various technical features that are associated with artifacts. With this information in hand, it can mitigate the effects of those artifacts during the duplicate marking and base recalibration steps. The GATK requires several read group fields to be present in input files and will fail with errors if this requirement is not satisfied. More on the fields can be found here [1]
```
samtools view -H sample.bam | grep '^@RG'
```


### Convert FASTQ to uBAM and add read group information using FastqToSAM

### Convert aligned BAM to uBAM and discard problematic records using RevertSam

### Reference
1. https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently

2. https://sites.google.com/a/broadinstitute.org/legacy-gatk-forum-discussions/tutorials/6484-how-to-generate-an-unmapped-bam-from-fastq-or-aligned-bam


[1] https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
