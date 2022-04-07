cwlVersion: v1.0
class: CommandLineTool
id: rnaseqc
label: RNA-SeQC 2.4.2
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'gcr.io/broad-cga-aarong-gtex/rnaseqc:2.4.2'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 10000

baseCommand: [rnaseqc]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.collapsed_gtf.path)
      $(inputs.Aligned_sorted_bam.path)

inputs:
  Aligned_sorted_bam: File
  collapsed_gtf: File
  unpaired: { type: 'boolean?', doc: "If single-end, set to true", default: false, inputBinding: { position: 1, prefix: "--unpaired"} }
  legacy: { type: 'boolean?', doc: "If true, will output format compatible with version 1.1.9", default: true, inputBinding: { position: 1, prefix: "--legacy" } }
  stranded: { type: [ 'null', {type: enum, name: rnaseqc_std, symbols: ["rf", "fr"]}],
  doc: "If stranded, specify", default: "rf", inputBinding: { position: 1, prefix: "--stranded=", separate: false} }

outputs:
  Metrics:
    type: File
    outputBinding:
      glob: 'output/*.metrics.tsv'
  Gene_TPM:
    type: File
    outputBinding:
      glob: 'output/*.gene_tpm.gct'
  Gene_count:
    type: File
    outputBinding:
      glob: 'output/*.gene_reads.gct'
  Exon_count:
    type: File
    outputBinding:
      glob: 'output/*.exon_reads.gct'
