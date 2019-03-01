cwlVersion: v1.0
class: Workflow
id: kfdrc_rnaseq_wf_bam_in
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  sample_name: string
  r1_adapter: {type: ['null', string]}
  r2_adapter: {type: ['null', string]}
  STAR_outSAMattrRGline: string
  STARgenome: File
  RSEMgenome: File
  reference_fasta: File
  gtf_anno: File
  wf_strand_param: {type: ['null', string], doc: "use 'default' or leave blank for unstranded/auto, rf_stranded if read1 in the fastq read pairs is reverse complement to the transcript, fr-stranded if read1 same sense as transcript"}
  FusionGenome: File
  runThread: int
  input_bam: File
  RNAseQC_GTF: File
  kallisto_idx: File
  pizzly_transcript_ref: File

outputs:
  cutadapt_stats: {type: File, outputSource: cutadapt/cutadapt_stats}
  STAR_transcriptome_bam: {type: File, outputSource: star/transcriptome_bam_out}
  STAR_sorted_genomic_bam: {type: File, outputSource: samtools_sort/sorted_bam}
  STAR_sorted_genomic_bai: {type: File, outputSource: samtools_sort/sorted_bai}
  STAR_chimeric_bam_out: {type: File, outputSource: samtools_sort/chimeric_bam_out}
  STAR_chimeric_junctions: {type: File, outputSource: star_fusion/chimeric_junction_compressed}
  STAR_gene_count: {type: File, outputSource: star/gene_counts}
  STAR_junctions_out: {type: File, outputSource: star/junctions_out}
  STAR_final_log: {type: File, outputSource: star/log_final_out}
  STAR-Fusion_results: {type: File, outputSource: star_fusion/abridged_coding}
  pizzly_fusion_results: {type: File, outputSource: pizzly/fusions_flattened}
  arriba_fusion_results: {type: File, outputSource: arriba_fusion/arriba_fusions}
  arriba_fusion_viz: {type: File, outputSource: arriba_fusion/arriba_pdf}
  RSEM_isoform: {type: File, outputSource: rsem/isoform_out}
  RSEM_gene: {type: File, outputSource: rsem/gene_out}
  RNASeQC_Metrics: {type: File, outputSource: rna_seqc/Metrics}
  RNASeQC_counts: {type: File, outputSource: supplemental/RNASeQC_counts}
  kallisto_Abundance: {type: File, outputSource: kallisto/abundance_out}

steps:

  bam2fastq:
    run: ../tools/bam2fastq.cwl
    in:
      input_bam: input_bam
      SampleID: sample_name
      runThreadN: runThread
    out:
      [fq1, fq2]

  cutadapt:
    run: ../tools/cutadapter.cwl
    in:
      readFilesIn1: bam2fastq/fq1
      readFilesIn2: bam2fastq/fq2
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

  samtools_sort:
    run: ../tools/samtools_sort.cwl
    in:
      unsorted_bam: star/genomic_bam_out
      chimeric_sam_out: star/chimeric_sam_out
    out:
      [sorted_bam, sorted_bai, chimeric_bam_out]

  strand_parse:
    run: ../tools/expression_parse_strand_param.cwl
    in:
      strand: wf_strand_param
    out:
      [
        rsem_std,
        kallisto_std,
        rnaseqc_std,
        arriba_std
      ]

  star_fusion:
    run: ../tools/STAR-Fusion.cwl
    in:
      readFilesIn1: cutadapt/trimmedReadsR1
      readFilesIn2: cutadapt/trimmedReadsR2
      Chimeric_junction: star/chimeric_junctions
      genomeDir: FusionGenome
      SampleID: sample_name
    out:
      [abridged_coding, chimeric_junction_compressed]

  arriba_fusion:
    run: ../tools/arriba.cwl
    in:
      genome_aligned_bam: samtools_sort/sorted_bam
      genome_aligned_bai: samtools_sort/sorted_bai
      chimeric_sam_out: star/chimeric_sam_out
      reference_fasta: reference_fasta
      gtf_anno: gtf_anno
      outFileNamePrefix: sample_name
      arriba_strand_flag: strand_parse/arriba_std
    out:
      [
        arriba_fusions,
        arriba_pdf
      ]

  rsem:
    run: ../tools/rsem-calculate-expression.cwl
    in:
      bam: star/transcriptome_bam_out
      genomeDir: RSEMgenome
      outFileNamePrefix: sample_name
      forward_prob: strand_parse/rsem_std
    out: [
      gene_out,
      isoform_out
    ]

  rna_seqc:
    run: ../tools/RNAseQC.cwl
    in:
      Aligned_sorted_bam: samtools_sort/sorted_bam
      collapsed_gtf: RNAseQC_GTF
      strand: strand_parse/rnaseqc_std
    out: [
      Metrics,
      Gene_TPM,
      Gene_count,
      Exon_count
    ]

  supplemental:
    run: ../tools/supplemental_tar_gz.cwl
    in:
      outFileNamePrefix: sample_name
      Gene_TPM: rna_seqc/Gene_TPM
      Gene_count: rna_seqc/Gene_count
      Exon_count: rna_seqc/Exon_count
    out: [
      RNASeQC_counts
    ]

  kallisto:
    run: ../tools/kallisto.cwl
    in:
      transcript_idx: kallisto_idx
      strand: strand_parse/kallisto_std
      reads1: cutadapt/trimmedReadsR1
      reads2: cutadapt/trimmedReadsR2
      SampleID: sample_name
    out: [
      abundance_out,
      fusion_out
    ]

  pizzly:
    run: ../tools/pizzly.cwl
    in:
      transcript_fa: pizzly_transcript_ref
      GTF: gtf_anno
      kallisto_fusion: kallisto/fusion_out
      SampleID: sample_name
    out: [fusions_flattened]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 3