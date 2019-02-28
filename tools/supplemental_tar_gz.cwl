cwlVersion: v1.0
class: CommandLineTool
id: supplemental_tar_gz
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'ubuntu:18.04'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 4
    ramMin: 1600

baseCommand: [mkdir]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.outFileNamePrefix)_STAR_supplemental

      cp $(inputs.log_final_out.path)
      $(inputs.junctions_out.path)
      $(inputs.chimeric_junctions.path)
      $(inputs.gene_counts.path)
      $(inputs.outFileNamePrefix)_STAR_supplemental

      tar -czf
      $(inputs.outFileNamePrefix).STAR.supplemental.tar.gz
      $(inputs.outFileNamePrefix)_STAR_supplemental

      mkdir $(inputs.outFileNamePrefix)_RNASeQC_supplemental

      cp $(inputs.Gene_TPM.path)
      $(inputs.Gene_count.path)
      $(inputs.Exon_count.path)
      $(inputs.outFileNamePrefix)_RNASeQC_supplemental

      tar -czf
      $(inputs.outFileNamePrefix).RNASeQC.supplemental.tar.gz
      $(inputs.outFileNamePrefix)_RNASeQC_supplemental

inputs:
  outFileNamePrefix: string
  log_final_out: File
  junctions_out: File
  chimeric_junctions: File
  gene_counts: File
  Gene_TPM: File
  Gene_count: File
  Exon_count: File

outputs:
  STAR_supplemental:
    type: File
    outputBinding:
      glob: "$(inputs.outFileNamePrefix).STAR.supplemental.tar.gz"
  RNASeQC_supplemental:
    type: File
    outputBinding:
      glob: "$(inputs.outFileNamePrefix).RNASeQC.supplemental.tar.gz"
