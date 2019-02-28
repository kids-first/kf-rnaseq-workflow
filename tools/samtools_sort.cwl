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
      $(inputs.unsorted_bam.nameroot).sorted.bai &&
      samtools view
      -bh
      -@ 16
      $(inputs.chimeric_sam_out.path)
      -o $(inputs.chimeric_sam_out.nameroot).bam

inputs:
  unsorted_bam: File
  chimeric_sam_out: File

outputs:
  sorted_bam:
    type: File
    outputBinding:
      glob: '*.sorted.bam'
  sorted_bai:
    type: File
    outputBinding:
      glob: '*.sorted.bai'
  chimeric_bam_out:
    type: File
    outputBinding:
      glob: "$(inputs.chimeric_sam_out.nameroot).bam"
