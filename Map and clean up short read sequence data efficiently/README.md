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

## Step 3. Align reads with BWA-MEM and merge with uBAM using MergeBamAlignment

Here we pipe 3 tools SamToFastq, BWA-MEM and MergeBamAlignment. By piping these we bypass storage of larger intermediate FASTQ and SAM files. We additionally save time by eliminating the need for the processor to read in and write out data for two of the processes, as piping retains data in the processor's input-output (I/O) device for the next process

### Step 3A. Convert BAM to FASTQ and discount adapter sequences using SamToFastq

We use additional options to effectively remove previously marked adapter sequences, in this example marked with an XT tag. By specifying CLIPPING_ATTRIBUTE=XT and CLIPPING_ACTION=2, SamToFastq changes the quality scores of bases marked by XT to two--a rather low score in the Phred scale. This effectively removes the adapter portion of sequences from contributing to downstream read alignment and alignment scoring metrics.

```
java -Xmx8G -jar /path/picard.jar SamToFastq \
I=6484_snippet_markilluminaadapters.bam \
FASTQ=6484_snippet_samtofastq_interleaved.fq \
CLIPPING_ATTRIBUTE=XT \
CLIPPING_ACTION=2 \
INTERLEAVE=true \ --PE reads
TMP_DIR=/path/shlee #optional to process large files         
```
### Step 3B. Align reads and flag secondary hits using BWA-MEM

BWA-MEM is suitable for aligning high-quality long reads ranging from 70 bp to 1 Mbp against a large reference genome such as the human genome. BWA alignment requires an indexed reference genome file. Indexing is specific to algorithms. To index the human genome for BWA, we apply BWA's index function on the reference genome file, e.g. human_g1k_v37_decoy.fasta. This produces five index files with the extensions amb, ann, bwt, pac and sa.

```
bwa index -a bwtsw human_g1k_v37_decoy.fasta -- creates the fasta

## Creates the alignment SAM file; Interleaved FASTQ has both the PE reads from the SamToFastq

/path/bwa mem -M -t 7 -p /path/human_g1k_v37_decoy.fasta \ 
6484_snippet_samtofastq_interleaved.fq > 6484_snippet_bwa_mem.sam
```

### Step 3C. Restore altered data and apply & adjust meta information using MergeBamAlignment

Broadly, the tool merges defined information from the unmapped BAM (uBAM, step 1) with that of the aligned BAM (step 3) to conserve read data, e.g. original read information and base quality scores. The tool also generates additional meta information based on the information generated by the aligner, which may alter aligner-generated designations, e.g. mate information and secondary alignment flags. The tool then makes adjustments so that all meta information is congruent, e.g. read and mate strand information based on proper mate designations. We ascribe the resulting BAM as clean.

Specifically, the aligned BAM generated in step 3 lacks read group information and certain tags--the UQ (Phred likelihood of the segment), MC (CIGAR string for mate) and MQ (mapping quality of mate) tags. It has hard-clipped sequences from split reads and altered base qualities. The reads also have what some call mapping artifacts but what are really just features we should not expect from our aligner. For example, the meta information so far does not consider whether pairs are optimally mapped and whether a mate is unmapped (in reality or for accounting purposes). Depending on these assignments, MergeBamAlignment adjusts the read and read mate strand orientations for reads in a proper pair. Finally, the alignment records are sorted by query name. We would like to fix all of these issues before taking our data to a variant discovery workflow.

```
java -Xmx16G -jar /path/picard.jar MergeBamAlignment \
R=/path/Homo_sapiens_assembly19.fasta \ 
UNMAPPED_BAM=6384_snippet_revertsam.bam \ 
ALIGNED_BAM=6484_snippet_bwa_mem.sam \ #accepts either SAM or BAM
O=6484_snippet_mergebamalignment.bam \
CREATE_INDEX=true \ #standard Picard option for coordinate-sorted outputs
ADD_MATE_CIGAR=true \ #default; adds MC tag
CLIP_ADAPTERS=false \ #changed from default
CLIP_OVERLAPPING_READS=true \ #default; soft-clips ends so mates do not extend past each other
INCLUDE_SECONDARY_ALIGNMENTS=true \ #default
MAX_INSERTIONS_OR_DELETIONS=-1 \ #changed to allow any number of insertions or deletions
PRIMARY_ALIGNMENT_STRATEGY=MostDistant \ #changed from default BestMapq
ATTRIBUTES_TO_RETAIN=XS \ #specify multiple times to retain tags starting with X, Y, or Z 
TMP_DIR=/path/shlee #optional to process large files
```
This generates the coordinate-sorted and clean BAM with it's index file .bai.

BAM file is ready for the analysis workflows that start with MarkDuplicates.

### Reference
1. https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently

2. https://sites.google.com/a/broadinstitute.org/legacy-gatk-forum-discussions/tutorials/6484-how-to-generate-an-unmapped-bam-from-fastq-or-aligned-bam


[1]: https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
