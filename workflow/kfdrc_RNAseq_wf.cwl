cwlVersion: v1.0
class: Workflow
id: kf_rnaseq_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  sample_name: string
  reads1: File
  reads2: File
  STARgenome: File
  RSEMgenome: File
  FusionGenome: File
  runThread: int
  STAR_outSAMattrRGline: string

outputs:
  STAR_transcriptome_bam: {type: File, outputSource: star/transcriptome_bam_out}
  STAR_junctions: {type: File, outputSource: star/junctions_out}
  STAR_genomic_bam: {type: File, outputSource: star/genomic_bam_out}
  STAR_gene_counts: {type: File, outputSource: star/gene_counts}
  STAR_chimeric_junctions: {type: File, outputSource: star/chimeric_junctions}
  STAR_chimeric_sam: {type: File, outputSource: star/chimeric_sam_out}
  STAR_Fusion: {type: File, outputSource: star_fusion/fusion_out}
  RSEM_isoform: {type: File, outputSource: rsem/isoform_out}
  RSEM_gene: {type: File, outputSource: rsem/gene_out}

steps:
  star:
    run: ../tools/star_align.d3b.2.cwl
    in:
      outSAMattrRGline: STAR_outSAMattrRGline
      readFilesIn1: reads1
      readFilesIn2: reads2
      genomeDir: STARgenome
      runThreadN: runThread
      outFileNamePrefix: sample_name
    out: [
      chimeric_junctions,
      chimeric_sam_out,
      gene_counts,
      genomic_bam_out,
      junctions_out,
    #   log_final_out,
    #   log_out,
    #   log_progress_out,
    #   star_pass1,
    #   star_pass1_genome,
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
      runThreadN: runThread
      outFileNamePrefix: sample_name
    out: [
      gene_out,
      isoform_out
    ]


$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:AWSInstanceType'
    value: c4.8xlarge;ebs-gp2;850
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4