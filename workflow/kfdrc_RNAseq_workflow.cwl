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
  reference_genome: File
  gtf_anno: File
  wf_strand_param: {type: string, doc: "use 'default' for unstranded/auto, rf_stranded if read1 in the fastq read pairs is reverse complement to the transcript, fr-stranded if read1 same sense as transcript"}
  FusionGenome: File
  runThread: int
  STAR_outSAMattrRGline: string
  RNAseQC_GTF: File
  kallisto_idx: File
  pizzly_transcript_ref: File

outputs:
  cutadapt_stats: {type: File, outputSource: cutadapt/cutadapt_stats}
  STAR_transcriptome_bam: {type: File, outputSource: star/transcriptome_bam_out}
  STAR_junctions: {type: File, outputSource: star/junctions_out}
  STAR_sorted_genomic_bam: {type: File, outputSource: samtools_sort/sorted_bam}
  STAR_sorted_genomic_bai: {type: File, outputSource: samtools_sort/sorted_bai}
  STAR_gene_counts: {type: File, outputSource: star/gene_counts}
  STAR_chimeric_junctions: {type: File, outputSource: star/chimeric_junctions}
  STAR_chimeric_sam: {type: File, outputSource: star/chimeric_sam_out}
  STAR-Fusion_results: {type: File, outputSource: star_fusion/abridged_coding}
  pizzly_fusion_results: {type: File, outputSource: pizzly/fusions_flattened}
  arriba_fusion_results: {type: File, outputSource: arriba_fusion/arriba_fusions}
  RSEM_isoform: {type: File, outputSource: rsem/isoform_out}
  RSEM_gene: {type: File, outputSource: rsem/gene_out}
  RNASeQC_Metrics: {type: File, outputSource: rna_seqc/Metrics}
  RNASeQC_Gene_TPM: {type: File, outputSource: rna_seqc/Gene_TPM}
  RNASeQC_Gene_count: {type: File, outputSource: rna_seqc/Gene_count}
  RNASeQC_Exon_count: {type: File, outputSource: rna_seqc/Exon_count}
  kallisto_Abundance: {type: File, outputSource: kallisto/abundance_out}
  kallisto_Fusion: {type: File, outputSource: kallisto/fusion_out}

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

  pizzly:
    run: ../tools/pizzly.cwl
    in:
      transcript_fa: pizzly_transcript_ref
      GTF: gtf_anno
      kallisto_fusion: kallisto/fusion_out
      SampleID: sample_name
    out: [fusions_flattened]

  star_fusion:
    run: ../tools/STAR-Fusion.cwl
    in:
      readFilesIn1: cutadapt/trimmedReadsR1
      readFilesIn2: cutadapt/trimmedReadsR2
      Chimeric_junction: star/chimeric_junctions
      genomeDir: FusionGenome
      SampleID: sample_name
    out:
      [abridged_coding]

  arriba_fusion:
    run: ../tools/arriba.cwl
    in:
      genome_aligned_bam: star/genomic_bam_out
      reference_fasta: reference_genome
      gtf_anno: gtf_anno
      outFileNamePrefix: sample_name
      arriba_strand_flag: strand_parse/arriba_std
    out:
      [arriba_fusions]

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

  samtools_sort:
    run: ../tools/samtools_sort.cwl
    in:
      unsorted_bam: star/genomic_bam_out
    out:
      [sorted_bam, sorted_bai]

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


$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 2