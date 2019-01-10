# kf-rnaseq-workflow
RNA-Seq workflow for Kids-First DRC

## Introduction
According to GDC's RNA-Seq workflow, this workflow is composed by several steps.

### STAR
RNA-Seq raw data alignment.
### RSEM
Calculation of gene expression.
### STAR-Fusion
Fusion detection for `STAR` chimeric reads.
### HTSeq
htseq-count tool for gene expression count using genomic bam file from `STAR`.
### RNA-SeQC
### Kallisto
Fusion detection using genomic bam file from `STAR`.
### Pizzly
Further fusion analysis tool and generate fusion fasta file from `Kallisto` result.


## Usage
Detail information can be found on Cavatica.
https://cavatica.sbgenomics.com/u/zhangb1/kf-rnaseq-workflow-test/apps/#zhangb1/kf-rnaseq-workflow-test/kf-rnaseq-wf

