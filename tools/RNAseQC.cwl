cwlVersion: v1.0
class: CommandLineTool
id: rnaseqc
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'gcr.io/broad-cga-aarong-gtex/rnaseqc:latest'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 10000

baseCommand: [rnaseqc]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.GTF.path)
      $(inputs.Aligned_bam.path)
      output/
      --legacy
      --stranded=rf
      --rpkm

inputs:
  Aligned_bam: File
  GTF: File
  GenomeRef:
    type: File
    secondaryFiles: [.fai, ^.dict]
  SampleID: string

outputs:
  Metrics:
    type: File
    outputBinding:
      glob: '*.metrics.tsv'

  Gene_TPM:
    type: File
    outputBinding:
      glob: '*.gene_tpm.gct'

  Gene_count:
    type: File
    outputBinding:
      glob: '*.gene_reads.gct'

  Exon_count:
    type: File
    outputBinding:
      glob: '*.exon_reads.gct'

