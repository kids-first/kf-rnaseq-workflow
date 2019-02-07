cwlVersion: v1.0
class: Workflow
id: kfdrc_rnaseq_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  sample_name: string
  r1_adapter: {type: ['null', string]}
  r2_adapter: {type: ['null', string]}
  reads1: File
  reads2: File
  STARgenome: File
  RSEMgenome: File
  FusionGenome: File
  runThread: int
  STAR_outSAMattrRGline: string
  RNAseQC_GTF: File
  GTF_Anno: File
  kallisto_idx: File
  pizzly_transcript_ref: File


outputs:
  cutadapt_stats: {type: File, outputSource: cutadapt/cutadapt_stats}
  STAR_transcriptome_bam: {type: File, outputSource: star/transcriptome_bam_out}
  STAR_junctions: {type: File, outputSource: star/junctions_out}
  STAR_genomic_bam: {type: File, outputSource: star/genomic_bam_out}
  STAR_gene_counts: {type: File, outputSource: star/gene_counts}
  STAR_chimeric_junctions: {type: File, outputSource: star/chimeric_junctions}
  STAR_chimeric_sam: {type: File, outputSource: star/chimeric_sam_out}
  STAR_Fusion: {type: File, outputSource: star_fusion/fusion_out}
  RSEM_isoform: {type: File, outputSource: rsem/isoform_out}
  RSEM_gene: {type: File, outputSource: rsem/gene_out}
  RNASeQC_Metrics: {type: File, outputSource: rna_seqc/Metrics}
  RNASeQC_Gene_TPM: {type: File, outputSource: rna_seqc/Gene_TPM}
  RNASeQC_Gene_count: {type: File, outputSource: rna_seqc/Gene_count}
  RNASeQC_Exon_count: {type: File, outputSource: rna_seqc/Exon_count}
  kallisto_Abundance: {type: File, outputSource: kallisto/abundance_out}
  kallisto_Fusion: {type: File, outputSource: kallisto/fusion_out}
  Pizzly_Fasta: {type: File, outputSource: pizzly/fusions_fasta}
  Pizzly_unfiltered_Fasta: {type: File, outputSource: pizzly/unfiltered_fusion_fasta}

steps:
  cutadapt:
    run: ../tools/cutadapter.cwl
    in:
      readFilesIn1: reads1
      readFilesIn2: reads2
      r1_adapter: r1_adapter
      r2_adapter: r2_adapter
      sample_name: sample_name
    out: [
    trimmedReadsR1,
    trimmedReadsR2,
    cutadapt_stats
    ]

  star:
    run: ../tools/star_align.cwl
    in:
      outSAMattrRGline: STAR_outSAMattrRGline
      readFilesIn1: cutadapt/trimmedReadsR1
      readFilesIn2: cutadapt/trimmedReadsR2
      genomeDir: STARgenome
      runThreadN: runThread
      outFileNamePrefix: sample_name
      r1_adapter: r1_adapter
      r2_adapter: r2_adapter

    out: [
      chimeric_junctions,
      chimeric_sam_out,
      gene_counts,
      genomic_bam_out,
      junctions_out,
      log_final_out,
      log_out,
      log_progress_out,
      transcriptome_bam_out
    ]

  star_fusion:
    run: ../tools/STAR-Fusion.cwl
    in:
      Chimeric: star/chimeric_junctions
      genomeDir: FusionGenome
      runThreadN: runThread
      SampleID: sample_name
    out:
      [fusion_out]

  rsem:
    run: ../tools/rsem-calculate-expression.cwl
    in:
      bam: star/transcriptome_bam_out
      genomeDir: RSEMgenome
      outFileNamePrefix: sample_name
    out: [
      gene_out,
      isoform_out
    ]

  samtools_sort:
    run: ../tools/samtools_sort.cwl
    in:
      unsorted_bam: star/genomic_bam_out
    out:
      [sorted_bam]

  rna_seqc:
    run: ../tools/RNAseQC.cwl
    in:
      Aligned_sorted_bam: samtools_sort/sorted_bam
      collapsed_gtf: RNAseQC_GTF
    out: [
      Metrics,
      Gene_TPM,
      Gene_count,
      Exon_count
    ]

  kallisto:
    run: ../tools/kallisto.cwl
    in:
      transcript_idx: kallisto_idx
      reads1: cutadapt/trimmedReadsR1
      reads2: cutadapt/trimmedReadsR2
      runThreadN: runThread
      SampleID: sample_name
    out: [
      abundance_out,
      fusion_out
    ]

  pizzly:
    run: ../tools/pizzly.cwl
    in:
      transcript_fa: pizzly_transcript_ref
      GTF: GTF_Anno
      kallisto_fusion: kallisto/fusion_out
      SampleID: sample_name
    out: [
      fusions_fasta,
      unfiltered_fusion_fasta
    ]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4