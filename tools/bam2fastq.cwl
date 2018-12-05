cwlVersion: v1.0
class: CommandLineTool
id: bam2fastq
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/picard:2.18.2-dev'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 10000

baseCommand: [java, -Dsamjdk.compression_level=2, -Xms4000m, -jar, /picard.jar, SamToFastq]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      I=$(inputs.Aligned_out.path)
      F=$(inputs.SampleID).kallisto_1.fastq.gz
      F2=$(inputs.SampleID).kallisto_2.fastq.gz

inputs:
  Aligned_out: File
  SampleID: string

outputs:
  kallisto_fq1:
    type: File
    outputBinding:
      glob: '*.kallisto_1.fastq.gz'

  kallisto_fq2:
    type: File
    outputBinding:
      glob: '*.kallisto_2.fastq.gz'

