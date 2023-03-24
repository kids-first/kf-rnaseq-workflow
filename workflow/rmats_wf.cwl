cwlVersion: v1.2
class: Workflow
id: rmats_wf
label: "rMATS Turbo"
doc: |
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
requirements:
- class: ScatterFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
inputs:
  gtf_annotation: {type: 'File', doc: "Input gtf annotation file."}
  sample_1_bams: {type: 'File[]', doc: "Input sample 1 bam files."}
  sample_2_bams: {type: 'File[]?', doc: "Input sample 2 bam files."}
  read_length: {type: 'int?', doc: "Input read length for sample reads."}
  variable_read_length: {type: 'boolean?', doc: "Allow reads with lengths that differ\
      \ from --readLength to be processed. --readLength will still be used to determine\
      \ IncFormLen and SkipFormLen."}
  read_type:
    type:
    - "null"
    - type: enum
      symbols:
      - paired
      - single
      name: read_type
    doc: "Select one option for input read type either paired or single. Tool default:\
      \ paired"
  strandedness:
    type:
    - "null"
    - type: enum
      symbols:
      - fr-unstranded
      - fr-firststrand
      - fr-secondstrand
      name: strandedness
    doc: "Select one option for input strandedness. Tool default: fr-unstranded"
  novel_splice_sites: {type: 'boolean?', doc: "Select for novel splice site detection\
      \ or unannotated splice sites. 'true' to detect or add this parameter, 'false'\
      \ to disable denovo detection. Tool Default: false"}
  stat_off: {type: 'boolean?', doc: "Select to skip statistical analysis, either between\
      \ two groups or on single sample group. 'true' to add this parameter. Tool default:\
      \ false"}
  allow_clipping: {type: 'boolean?', doc: "Allow alignments with soft or hard clipping\
      \ to be used."}
  output_basename: {type: 'string', doc: "String to use as basename for output files"}
  rmats_threads: {type: 'int?', doc: "Threads to allocate to RMATs."}
  rmats_ram: {type: 'int?', doc: "GB of RAM to allocate to RMATs."}
outputs:
  filtered_alternative_3_prime_splice_sites_jc: {type: 'File', outputSource: filter_alt_3_prime/output,
    doc: "Alternative 3 prime splice sites JC.txt output from RMATs containing only\
      \ those calls with 10 or more read counts of support"}
  filtered_alternative_5_prime_splice_sites_jc: {type: 'File', outputSource: filter_alt_5_prime/output,
    doc: "Alternative 5 prime splice sites JC.txt output from RMATs containing only\
      \ those calls with 10 or more read counts of support"}
  filtered_mutually_exclusive_exons_jc: {type: 'File', outputSource: filter_me_exons/output,
    doc: "Mutually exclusive exons JC.txt output from RMATs containing only those\
      \ calls with 10 or more read counts of support"}
  filtered_retained_introns_jc: {type: 'File', outputSource: filter_retained_introns/output,
    doc: "Retained introns JC.txt output from RMATs containing only those calls with\
      \ 10 or more read counts of support"}
  filtered_skipped_exons_jc: {type: 'File', outputSource: filter_skipped_exons/output,
    doc: "Skipped exons JC.txt output from RMATs containing only those calls with\
      \ 10 or more read counts of support"}
steps:
  samtools_readlength_bam:
    run: ../tools/samtools_readlength_bam.cwl
    in:
      input_bam:
        source: sample_1_bams
        valueFrom: |
          $(self[0])
    out: [output, top_readlength, variable_readlength]
  rmats_both_bam:
    run: ../tools/rmats_both_bam.cwl
    in:
      gtf_annotation: gtf_annotation
      sample_1: sample_1_bams
      sample_2: sample_2_bams
      read_length:
        source: [read_length, samtools_readlength_bam/top_readlength]
        pickValue: first_non_null
      variable_read_length:
        source: [variable_read_length, samtools_readlength_bam/variable_readlength]
        pickValue: first_non_null
      read_type: read_type
      strandedness: strandedness
      allow_clipping: allow_clipping
      novel_splice_sites: novel_splice_sites
      stat_off: stat_off
      output_directory: output_basename
      threads: rmats_threads
      ram: rmats_ram
    out: [alternative_3_prime_splice_sites_jc, alternative_5_prime_splice_sites_jc,
      mutually_exclusive_exons_jc, retained_introns_jc, skipped_exons_jc, temp_read_outcomes,
      summary_file]
  filter_alt_3_prime:
    run: ../tools/awk_junction_filtering.cwl
    in:
      input_jc_file: rmats_both_bam/alternative_3_prime_splice_sites_jc
    out: [output]
  filter_alt_5_prime:
    run: ../tools/awk_junction_filtering.cwl
    in:
      input_jc_file: rmats_both_bam/alternative_5_prime_splice_sites_jc
    out: [output]
  filter_me_exons:
    run: ../tools/awk_junction_filtering.cwl
    in:
      input_jc_file: rmats_both_bam/mutually_exclusive_exons_jc
    out: [output]
  filter_retained_introns:
    run: ../tools/awk_junction_filtering.cwl
    in:
      input_jc_file: rmats_both_bam/retained_introns_jc
    out: [output]
  filter_skipped_exons:
    run: ../tools/awk_junction_filtering.cwl
    in:
      input_jc_file: rmats_both_bam/skipped_exons_jc
    out: [output]
sbg:license: Apache License 2.0
sbg:publisher: KFDRC
