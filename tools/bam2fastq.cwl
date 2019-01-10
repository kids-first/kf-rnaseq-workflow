cwlVersion: v1.0
class: CommandLineTool
id: bam2fastq
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/samtools:1.9'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 36
    ramMin: 30000

baseCommand: [samtools, sort]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -n
      -@ $(inputs.runThreadN)
      -m 1G
      -O SAM
      $(inputs.input_bam.path) |
      samtools
      fastq
      -1 $(inputs.SampleID).converted_1.fastq -2 $(inputs.SampleID).converted_2.fastq -@ $(inputs.runThreadN) - &&
      ls ./*.fastq | xargs -IFN -P 2 gzip FN

inputs:
  input_bam: File
  SampleID: string
  runThreadN: int

outputs:
  fq1:
    type: File
    outputBinding:
      glob: '*.converted_1.fastq.gz'

  fq2:
    type: File
    outputBinding:
      glob: '*.converted_2.fastq.gz'

