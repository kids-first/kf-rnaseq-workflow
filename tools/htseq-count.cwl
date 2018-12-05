cwlVersion: v1.0
class: CommandLineTool
id: htseq_count
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'dmccloskey/htseq-count'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 10000

baseCommand: [htseq-count]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -f bam
      $(inputs.Aligned_bam.path)
      $(inputs.GTF.path)
      > $(inputs.SampleID).HTSeq-count.txt

inputs:
  Aligned_bam: File
  GTF: File
  SampleID: string

outputs:
  HTSeq_count:
    type: File
    outputBinding:
      glob: '*.HTSeq-count.txt'

