cwlVersion: v1.0
class: Workflow
id: temp_rnaseqc_only_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  sample_name: string
  wf_strand_param: {type: ['null', string], doc: "use 'default' for unstranded/auto, rf_stranded if read1 in the fastq read pairs is reverse complement to the transcript, fr-stranded if read1 same sense as transcript"}
  RNAseQC_GTF: File
  genome_aligned_bam: {type: File, secondaryFiles: ['^.bai']}

outputs:
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
  rna_seqc:
    run: ../../tools/rnaseqc.cwl
    in:
      Aligned_sorted_bam: genome_aligned_bam
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
# $namespaces:
#   sbg: https://sevenbridges.com
# hints:
#   - class: 'sbg:AWSInstanceType'
#     value: r4.8xlarge;ebs-gp2;200
