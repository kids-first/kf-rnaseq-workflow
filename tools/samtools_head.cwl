class: CommandLineTool
cwlVersion: v1.2
id: samtools_head
doc: |-
  By default, prints all headers from the specified input file to standard
  output in SAM format. The input alignment file may be in SAM, BAM, or CRAM
  format; if no FILE is specified, standard input will be read. With appropriate
  options, only some of the headers and/or additionally some of the alignment
  records will be printed.
  
  The samtools head command outputs SAM headers exactly as they appear in the
  input file; in particular, it never adds an @PG header itself. (Other samtools
  commands add such @PG headers to facilitate provenance tracking in analysis
  pipelines, but because samtools head never outputs more than a handful of
  alignment records it is unsuitable for use in such contexts anyway.) 

  Tool can also has the option to run an additional grep filter on the output
  of samtools head.
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'staphb/samtools:1.20'
baseCommand: [samtools, head]
arguments:
  - position: 5
    shellQuote: false
    valueFrom: >
      > header.txt
inputs:
  input_bam: { type: File, inputBinding: { position: 3 }, doc: "Input bam file" }
  num_records: { type: 'int?', inputBinding: { position: 2, prefix: "--records" }, doc: "Number of record lines to output" }
  num_headers: { type: 'int?', inputBinding: { position: 2, prefix: "--headers" }, doc: "Number of header lines to output" }
  line_filter: { type: 'string?', inputBinding: { position: 4, shellQuote: false, prefix: "| grep" }, doc: "Additional grep filter for samtools head output" } 
outputs:
  header_file:
    type: File
    outputBinding:
      glob: header.txt

$namespaces:
  sbg: https://sevenbridges.com
