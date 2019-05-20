cwlVersion: v1.0
class: CommandLineTool
id: arriba_fusion
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/arriba:1.1.0'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 5000

baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      /arriba_v1.1.0/draw_fusions.R
      --annotation=$(inputs.gtf_anno.path)
      --fusions=$(inputs.fusion_tsv.path)
      --alignments=$(inputs.genome_aligned_bam.path)
      --cytobands=/arriba_v1.1.0/database/cytobands_hg38_GRCh38_2018-02-23.tsv
      --proteinDomains=/arriba_v1.1.0/database/protein_domains_hg38_GRCh38_2018-03-06.gff3
      --output=$(inputs.outFileNamePrefix).arriba.fusions.pdf

inputs:
  genome_aligned_bam:
    type: File
    secondaryFiles: ^.bai
  fusion_tsv: File
  gtf_anno: File
  outFileNamePrefix: string

outputs:

  arriba_pdf:
    type: File
    outputBinding:
      glob: "$(inputs.outFileNamePrefix).arriba.fusions.pdf"
