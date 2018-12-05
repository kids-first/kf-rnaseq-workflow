cwlVersion: v1.0
class: CommandLineTool
id: rsem-calculate-expression
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'jinh2/star-rsem'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 10000

baseCommand: [tar]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -zxf $(inputs.genomeDir.path) &&
      /usr/local/RSEM-1.3.1/rsem-calculate-expression
      --paired-end
      --alignments
      --no-bam-output
      -p $(inputs.runThreadN)
      $(inputs.bam.path)
      ./STAR-RSEM-Ref/Homo_sapiens_assembly38
      $(inputs.outFileNamePrefix)

inputs:
  bam: File
  genomeDir: File
  runThreadN: int
  outFileNamePrefix: string

outputs:
  gene_out:
    type: File
    outputBinding: 
      glob: '*genes.results'

  isoform_out:
    type: File
    outputBinding:
      glob: '*isoforms.results'

