# Data pre-processing for variant discovery

<img width="389" alt="Screen Shot 2021-11-20 at 6 55 29 PM" src="https://user-images.githubusercontent.com/31465978/142744295-159aaabf-c881-431f-8fdd-272d4a2ccd35.png">

The first 3 steps can be found [here][1].

## Mark Duplicates
SortSam + MarkDuplicates.

Use SortSam to Sort the BAM file
+ INPUT: The BAM or SAM file to sort. Required.
+ OUTPUT: The sorted BAM or SAM output file. Required.
+ SORT_ORDER: Sort order of output file Required. Possible values: {unsorted, queryname, coordinate, duplicate}
+ VALIDATION_STRINGENCY: Validation stringency for all SAM files read by this program. Possible values: {STRICT, LENIENT, SILENT}
```
java -Xmx8G -jar $PICARD/picard-2.8.0.jar SortSam \
INPUT= aln-pe.bam \ -- SAM or BAM file
OUTPUT=aln-pe_sorted.bam \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=SILENT
```

Picard MarkDupicates function to mark the duplicate reads
+ INPUT: The sorted BAM or SAM file to sort. Required.
+ OUTPUT: The BAM or SAM output file. Required.
+ METRICS_FILE: File to write duplication metrics to Required.
+ ASSUME_SORTED: If true, assume that the input file is coordinate sorted even if the header says otherwise. Default value: false. Possible values: {true, false}
+ VALIDATION_STRINGENCY: Validation stringency for all SAM files read by this program. Default value: STRICT. Possible values: {STRICT, LENIENT, SILENT}
```
java -Xmx8G -jar $PICARD/picard-2.8.0.jar MarkDuplicates \
INPUT= aln-pe_sorted.bam \
OUTPUT= aln-pe_sorted_marked.bam \ -- Outputs a sorted and marked BAM file
METRICS_FILE=metrics.txt \
ASSUME_SORTED=true \
VALIDATION_STRINGENCY=SILENT
```

## Recalibrate Base Quality Scores

BQSR stands for Base Quality Score Recalibration. In a nutshell, it is a data pre-processing step that detects systematic errors made by the sequencing machine when it estimates the accuracy of each base call.

Tools involved are BaseRecalibrator, Apply Recalibration, AnalyzeCovariates(optional). AnalyzeCovariates creates plots and analysis.

Workflow
```
gatk BaseRecalibrator \
   -I my_reads.bam \
   -R reference.fasta \
   --known-sites sites_of_variation.vcf \
   --known-sites another/optional/setOfSitesToMask.vcf \
   -O recal_data.table
   
 gatk ApplyBQSR \
   -R reference.fasta \
   -I input.bam \
   --bqsr-recal-file recal_data.table \
   -O output.bam
```
Output bam is ready for analysis.


[1]: https://github.com/marianakhoul/GATK/tree/main/Map%20and%20clean%20up%20short%20read%20sequence%20data%20efficiently

### Referencess
1. https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-
