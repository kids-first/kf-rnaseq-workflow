cwlVersion: v1.0
class: CommandLineTool
id: rnaseqc
label: RNA-SeQC 2.4.2
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/rnaseqc:v2.4.2'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 16000

baseCommand: [rnaseqc]

inputs:
  collapsed_gtf: { type: File, doc: "Collapsed GTF file",
    inputBinding: { position: 0 } }
  aligned_sorted_reads: { type: File, doc: "Aligned and sorted BAM or CRAM file",
    inputBinding: { position: 1} }
  output_dirname: { type: 'string?', default: "output", doc: "output dirname",
    inputBinding: { position: 2} }
  stranded: { type: [ 'null', {type: enum, name: rnaseqc_std, symbols: ["rf", "fr"]}], doc: "If stranded, specify", default: "rf",
    inputBinding: { position: 3, prefix: "--stranded=", separate: false} }
  unpaired: { type: 'boolean?', doc: "If single-end, set to true", default: false,
    inputBinding: { position: 3, prefix: "--unpaired"} }
  legacy: { type: 'boolean?', doc: "If true, will output format compatible with version 1.1.9", default: true,
    inputBinding: { position: 3, prefix: "--legacy" } }
  bed: { type: 'File?', doc: "BED file with intervals for estimating insert size distribution",
    inputBinding: { position: 3, prefix: "--bed" } }
  fasta: { type: 'File?', secondaryFiles: ['.fai'], doc: "If input is CRAM, provide reference",
    inputBinding: { position: 3, prefix: "--fasta"} }
  logging: { type: 'boolean?', doc: "Turn on this flag to get progress updates",
    inputBinding: { position: 3, prefix: "-v"} }

outputs:
  Metrics:
    type: File
    outputBinding:
      glob: '$(inputs.output_dirname)/*.metrics.tsv'
  Gene_TPM:
    type: File
    outputBinding:
      glob: '$(inputs.output_dirname)/*.gene_tpm.gct'
  Gene_count:
    type: File
    outputBinding:
      glob: '$(inputs.output_dirname)/*.gene_reads.gct'
  Exon_count:
    type: File
    outputBinding:
      glob: '$(inputs.output_dirname)/*.exon_reads.gct'
