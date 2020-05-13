cwlVersion: v1.0
class: CommandLineTool
id: kallisto
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/uros_sipetic/kallisto:0.43.1'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 10000

baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      kallisto quant
      -i $(inputs.transcript_idx.path)
      -o output
      -b 10
      -t 8
      $(inputs.strand ? inputs.strand == "default" ? "" : "--"+inputs.strand : "")
      $(inputs.reads2 ? inputs.reads1.path+" "+inputs.reads2.path : "--single -l "+inputs.avg_frag_len+" -s "+inputs.std_dev+" "+inputs.reads1.path) &&
      mv output/abundance.tsv $(inputs.SampleID).kallisto.abundance.tsv &&
      gzip $(inputs.SampleID).kallisto.abundance.tsv

inputs:
  transcript_idx: File
  strand: {type: ['null', string], doc: "input none if unstranded, otherwise rf-stranded or fr-stranded"}
  reads1: File
  reads2: File?
  SampleID: string
  std_dev: long?
  avg_frag_len: int?

outputs:
  abundance_out:
    type: File
    outputBinding:
      glob: '*.abundance.tsv.gz'
