cwlVersion: v1.0
class: Workflow
id: temp_strand_ct_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  sample_name: string
  tx_bam: File
  gbam: {type: File, secondaryFiles: ['^.bai']}
  RSEMgenome: File
  wf_strand_param: {type: ['null', string], doc: "use 'default' for unstranded/auto, rf_stranded if read1 in the fastq read pairs is reverse complement to the transcript, fr-stranded if read1 same sense as transcript"}
  RNAseQC_GTF: File

outputs:
  RSEM_isoform: {type: File, outputSource: rsem/isoform_out}
  RSEM_gene: {type: File, outputSource: rsem/gene_out}
  RNASeQC_Metrics: {type: File, outputSource: rna_seqc/Metrics}
  RNASeQC_counts: {type: File, outputSource: supplemental/RNASeQC_counts}

steps:
  strand_parse:
    run: ../../tools/expression_parse_strand_param.cwl
    in:
      strand: wf_strand_param
    out:
      [
        rsem_std,
        kallisto_std,
        rnaseqc_std,
        arriba_std
      ]

  rsem:
    run: ../../tools/rsem_calc_expression.cwl
    in:
      bam: tx_bam
      genomeDir: RSEMgenome
      outFileNamePrefix: sample_name
      strandedness: strand_parse/rsem_std
    out: [
      gene_out,
      isoform_out
    ]

  rna_seqc:
    run: ../../tools/rnaseqc.cwl
    in:
      Aligned_sorted_bam: gbam
      collapsed_gtf: RNAseQC_GTF
      strand: strand_parse/rnaseqc_std
    out: [
      Metrics,
      Gene_TPM,
      Gene_count,
      Exon_count
    ]

  supplemental:
    run: ../../tools/supplemental_tar_gz.cwl
    in:
      outFileNamePrefix: sample_name
      Gene_TPM: rna_seqc/Gene_TPM
      Gene_count: rna_seqc/Gene_count
      Exon_count: rna_seqc/Exon_count
    out: [
      RNASeQC_counts
    ]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 3