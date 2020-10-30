cwlVersion: v1.0
class: CommandLineTool
id: samtools_flagstat
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 4
    ramMin: 4000

baseCommand: [samtools]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      flagstat $(inputs.input_bam.path)
      -@ 4
      > $(inputs.input_bam.nameroot).flagstat

inputs:
  input_bam: File

outputs:
  flagstat:
    type: File
    outputBinding:
      glob: '*.flagstat'
