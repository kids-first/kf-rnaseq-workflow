cwlVersion: v1.0
class: CommandLineTool
id: samtools_sort
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/samtools:1.9'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 24000

baseCommand: [samtools]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      sort $(inputs.unsorted_bam.path)
      -@ 8
      -m 2G
      -O bam
      > $(inputs.unsorted_bam.nameroot).sorted.bam

inputs:
  unsorted_bam: File

outputs:
  sorted_bam:
    type: File
    outputBinding:
      glob: '*.sorted.bam'
