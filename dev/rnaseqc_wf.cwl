cwlVersion: v1.2
class: Workflow
id: rnaseqc-wf
doc: "Temp workflow to run just RNAseQC and compress outputs"
requirements:
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement

inputs:
  collapsed_gtf: { type: File, doc: "Collapsed GTF file" }
  aligned_sorted_reads: { type: File, doc: "Aligned and sorted BAM or CRAM file" }
  outFileNamePrefix: { type: string, doc: "For naming tar gz output of counts" }
  stranded: { type: [ 'null', {type: enum, name: rnaseqc_std, symbols: ["rf", "fr"]}], doc: "If stranded, specify", default: "rf" }
  unpaired: { type: 'boolean?', doc: "If single-end, set to true", default: false }
  legacy: { type: 'boolean?', doc: "If true, will output format compatible with version 1.1.9", default: false }
  bed: { type: 'File?', doc: "BED file with intervals for estimating insert size distribution" }
  fasta: { type: 'File?', secondaryFiles: ['.fai'], doc: "If input is CRAM, provide reference" }
  logging: { type: 'boolean?', doc: "Turn on this flag to get progress updates" }
outputs:
  RNASeQC_Metrics: {type: 'File', outputSource: rna_seqc/Metrics, doc: "Metrics on
      mapping, intronic, exonic rates, count information, etc"}
  RNASeQC_counts: {type: 'File', outputSource: supplemental/RNASeQC_counts, doc: "Contains
      gene tpm, gene read, and exon counts" }
steps:
  rna_seqc:
    run: ../tools/rnaseqc_2.4.2.cwl
    in:
      aligned_sorted_reads: aligned_sorted_reads
      collapsed_gtf: collapsed_gtf
      stranded: stranded
      unpaired: unpaired
      bed: bed
      fasta: fasta
      legacy: legacy
      logging: logging
    out: [Metrics, Gene_TPM, Gene_count, Exon_count]
  supplemental:
    run: ../tools/supplemental_tar_gz.cwl
    in:
      outFileNamePrefix: outFileNamePrefix
      Gene_TPM: rna_seqc/Gene_TPM
      Gene_count: rna_seqc/Gene_count
      Exon_count: rna_seqc/Exon_count
    out: [RNASeQC_counts]
