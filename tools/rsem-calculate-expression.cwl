cwlVersion: v1.0
class: CommandLineTool
id: rsem-calculate-expression
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/uros_sipetic/rsem:1.3.1'
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
      rsem-calculate-expression
      --paired-end
      --alignments
      --append-names
      --no-bam-output
      -p 8
      $(inputs.bam.path)
      ./$(inputs.genomeDir.nameroot.split('.')[0])/$(inputs.genomeDir.nameroot.split('.')[0])
      $(inputs.outFileNamePrefix)

inputs:
  bam: File
  genomeDir: File
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

