cwlVersion: v1.0
class: CommandLineTool
id: samtools_sort
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/samtools:1.9'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 16
    ramMin: 24000

baseCommand: [samtools]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      sort $(inputs.unsorted_bam.path)
      -@ 16
      -m 1G
      -O bam
      > $(inputs.unsorted_bam.nameroot).sorted.bam &&
      samtools
      index
      -@ 16
      $(inputs.unsorted_bam.nameroot).sorted.bam
      $(inputs.unsorted_bam.nameroot).sorted.bai

inputs:
  unsorted_bam: File

outputs:
  sorted_bam:
    type: File
    outputBinding:
      glob: '*.sorted.bam'
  sorted_bai:
    type: File
    outputBinding:
      glob: '*.sorted.bai'
