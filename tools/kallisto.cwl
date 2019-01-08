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

baseCommand: [kallisto]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      quant -i $(inputs.transcript_idx.path)
      -o output --fusion
      -b 10 -t $(inputs.runThreadN)
      $(inputs.reads1.path)
      $(inputs.reads2.path) &&
      mv output/abundance.tsv $(inputs.SampleID).abundance.tsv &&
      mv output/fusion.txt $(inputs.SampleID).fusion.txt

inputs:
  transcript_idx: File
  GTF: File
  reads1: File
  reads2: File
  runThreadN: int
  SampleID: string

outputs:
  abundance_out:
    type: File
    outputBinding:
      glob: '*.abundance.tsv'

  fusion_out:
    type: File
    outputBinding:
      glob: '*.fusion.txt'

