# D3b rMATS Workflow

## Introduction

The rMATS workflow can also be run as a standalone workflow. In this workflow, rMATS is run on the input BAM files to generate 5 junction files: `[alternative_3_prime_splice_sites_jc, alternative_5_prime_splice_sites_jc, mutually_exclusive_exons_jc, retained_introns_jc, skipped_exons_jc]`. The workflow next grabs the sample information from the `sample_1_bams` by parsing the read group information from the BAM header for use in the output names. Each of the five junction files then undergo a simple filtering process where calls that have junction counts less than 10 are removed. These filtered junction files are returned as the final outputs.

## Usage

### Inputs

 - `gtf_annotation`: Input gtf annotation file
 - `sample_1_bams:`: Input sample 1 bam files
 - `sample_2_bams:`: Input sample 2 bam files
 - `read_length:`: Input read length for sample reads
 - `variable_read_length`: Allow reads with lengths that differ from --readLength to be processed. --readLength will still be used to determine IncFormLen and SkipFormLen
 - `read_type`: Select one option for input read type either paired or single. Tool default: paired
 - `strandedness`: Select one option for input strandedness. Tool default: fr-unstranded
 - `novel_splice_sites:`: Select for novel splice site detection or unannotated splice sites. 'true' to detect or add this parameter, 'false' to disable denovo detection. Tool Default: false
 - `stat_off:`: Select to skip statistical analysis, either between two groups or on single sample group. 'true' to add this parameter. Tool default: false
 - `allow_clipping:`: Allow alignments with soft or hard clipping to be used
 - `output_basename:`: String to use as basename for output files
 - `rmats_threads:`: Threads to allocate to RMATs
 - `rmats_ram:`: GB of RAM to allocate to RMATs

### Outputs

 - `filtered_alternative_3_prime_splice_sites_jc`: File extension `filtered.A3SS.MATS.JC.txt`. Alternative 3 prime splice sites JC.txt output from RMATs containing only those calls with 10 or more junction spanning read counts of support
 - `filtered_alternative_5_prime_splice_sites_jc`: File extension `filtered.A5SS.MATS.JC.txt`. Alternative 5 prime splice sites JC.txt output from RMATs containing only those calls with 10 or more junction spanning read counts of support
 - `filtered_mutually_exclusive_exons_jc`: File extension `filtered.MXE.MATS.JC.txt`. Mutually exclusive exons JC.txt output from RMATs containing only those calls with 10 or more junction spanning read counts of support
 - `filtered_retained_introns_jc`: File extension `filtered.RI.MATS.JC.txt`. Retained introns JC.txt output from RMATs containing only those calls with 10 or more junction spanning read counts of support
 - `filtered_skipped_exons_jc`: File extension `filtered.SE.MATS.JC.txt`. Skipped exons JC.txt output from RMATs containing only those calls with 10 or more junction spanning read counts of support
