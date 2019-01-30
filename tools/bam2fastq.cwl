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
      -m 1G -n -o -O SAM -@ $(inputs.runThreadN) $(inputs.input_bam.path) &&
      samtools fastq -1 $(inputs.SampleID).converted_1.fastq.gz -2 $(inputs.SampleID).converted_2.fastq.gz -@ $(inputs.runThreadN) -

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
