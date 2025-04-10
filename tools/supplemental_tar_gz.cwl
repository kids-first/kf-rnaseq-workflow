cwlVersion: v1.0
class: CommandLineTool
id: supplemental_tar_gz
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/ubuntu:18.04'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 4
    ramMin: 8000

baseCommand: [mkdir]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.outFileNamePrefix)_RNASeQC_counts

      cp $(inputs.Gene_TPM.path)
      $(inputs.Gene_count.path)
      $(inputs.Exon_count.path)
      $(inputs.outFileNamePrefix)_RNASeQC_counts

      tar -czf
      $(inputs.outFileNamePrefix).RNASeQC.counts.tar.gz
      $(inputs.outFileNamePrefix)_RNASeQC_counts

inputs:
  outFileNamePrefix: string
  Gene_TPM: File
  Gene_count: File
  Exon_count: File

outputs:
  RNASeQC_counts:
    type: File
    outputBinding:
      glob: "$(inputs.outFileNamePrefix).RNASeQC.counts.tar.gz"
