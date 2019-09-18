cwlVersion: v1.0
class: Workflow
doc: "Kids First Data Resource Center RNAseq Workflow (fastq input). This workflow follows the GDC outline for analysis in [mRNA Analysis Pipeline.](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/) The related Github branch for this app is located [here.](https://github.com/kids-first/kf-rnaseq-workflow/tree/be-publish) This workflow takes fastq input and uses STAR for alignment to the reference fasta version GRCh38. STAR gives aligned genomic, chimeric, transcriptomic, and junction aligned read outputs. Arriba and STAR fusion mode are run for fusion estimations on STAR alignment chimeric output. RSEM is run for gene expression estimations on STAR transcriptomic output. RNAseQC is run to generate a number of metrics including mapping rates, transcript counts, gene counts, sense and antisense mapping along with others. Gene abundance estimations are gathered from Kallisto run on trimmed reads of the raw read input."
id: kfdrc_rnaseq_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  sample_name: {type: string, doc: "Sample name used for file base name of all outputs"}
  r1_adapter: {type: ['null', string], doc: "Optional input. If the input reads have already been trimmed, leave these as null. If they do need trimming, supply the adapters."}
  r2_adapter: {type: ['null', string], doc: "Optional input. If the input reads have already been trimmed, leave these as null. If they do need trimming, supply the adapters."}
  reads1: {type: File, doc: "Input reads 1 file"}
  reads2: {type: File, doc: "Input reads 2 file"}
  STARgenome: {type: File, doc: "STAR_GENCODE27.tar.gz"}
  RSEMgenome: {type: File, doc: "RSEM_GENCODE27.tar.gz"}
  reference_fasta: {type: File, doc: "GRCh38.primary_assembly.genome.fa"}
  gtf_anno: {type: File, doc: "gencode.v27.primary_assembly.annotation.gtf"}
  wf_strand_param: {type: ['null', string], doc: "use 'default' for unstranded/auto, rf_stranded if read1 in the fastq read pairs is reverse complement to the transcript, fr-stranded if read1 same sense as transcript"}
  FusionGenome: {type: File, doc: "GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz"}
  runThread: {type: int, doc: "Amount of threads for analysis."}
  STAR_outSAMattrRGline: {type: string, doc: "Suggested setting, with TABS SEPARATING THE TAGS, format is: ID:sample_name LB:aliquot_id PL:platform SM:BSID for example ID:7316-242 LB:750189 PL:ILLUMINA SM:BS_W72364MN"}
  RNAseQC_GTF: {type: File, doc: "gencode.v27.primary_assembly.RNAseQC.gtf"}
  kallisto_idx: {type: File, doc: "gencode.v27.kallisto.index"}

outputs:
  cutadapt_stats: {type: File, outputSource: cutadapt/cutadapt_stats, doc: "Cutadapt stats output, only if adapter is supplied."}
  STAR_transcriptome_bam: {type: File, outputSource: star/transcriptome_bam_out, doc: "STAR bam of transcriptome reads"}
  STAR_sorted_genomic_bam: {type: File, outputSource: samtools_sort/sorted_bam, doc: "STAR sorted alignment bam"}
  STAR_sorted_genomic_bai: {type: File, outputSource: samtools_sort/sorted_bai, doc: "STAR index for sorted aligned bam"}
  STAR_chimeric_bam_out: {type: File, outputSource: samtools_sort/chimeric_bam_out, doc: "STAR bam output of chimeric reads"}
  STAR_chimeric_junctions: {type: File, outputSource: star_fusion/chimeric_junction_compressed, doc: "STAR chimeric junctions"}
  STAR_gene_count: {type: File, outputSource: star/gene_counts, doc: "STAR gene counts"}
  STAR_junctions_out: {type: File, outputSource: star/junctions_out, doc: "STAR junction reads"}
  STAR_final_log: {type: File, outputSource: star/log_final_out, doc: "STAR metrics log file of unique, multi-mapping, unmapped, and chimeric reads"}
  STAR-Fusion_results: {type: File, outputSource: star_fusion/abridged_coding, doc: "STAR fusion detection from chimeric reads"}
  arriba_fusion_results: {type: File, outputSource: arriba_fusion/arriba_fusions, doc: "Fusion output from Arriba"}
  arriba_fusion_viz: {type: File, outputSource: arriba_fusion/arriba_pdf, doc: "pdf output from Arriba"}
  RSEM_isoform: {type: File, outputSource: rsem/isoform_out, doc: "RSEM isoform expression estimates"}
  RSEM_gene: {type: File, outputSource: rsem/gene_out, doc: "RSEM gene expression estimates"}
  RNASeQC_Metrics: {type: File, outputSource: rna_seqc/Metrics, doc: "Metrics on mapping, intronic, exonic rates, count information, etc"}
  RNASeQC_counts: {type: File, outputSource: supplemental/RNASeQC_counts, doc: "Contains gene tpm, gene read, and exon counts"}
  kallisto_Abundance: {type: File, outputSource: kallisto/abundance_out, doc: "Gene abundance output from STAR genomic bam file"}

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

  samtools_sort:
    run: ../tools/samtools_sort.cwl
    in:
      unsorted_bam: star/genomic_bam_out
      chimeric_sam_out: star/chimeric_sam_out
    out:
      [sorted_bam, sorted_bai]

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
    run: ../tools/star_fusion.cwl
    in:
      Chimeric_junction: star/chimeric_junctions
      genomeDir: FusionGenome
      SampleID: sample_name
    out:
      [abridged_coding, chimeric_junction_compressed]

  arriba_fusion:
    run: ../tools/arriba_fusion.cwl
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
    run: ../tools/rsem_calc_expression.cwl
    in:
      bam: star/transcriptome_bam_out
      genomeDir: RSEMgenome
      outFileNamePrefix: sample_name
      strandedness: strand_parse/rsem_std
    out: [
      gene_out,
      isoform_out
    ]

  rna_seqc:
    run: ../tools/rnaseqc.cwl
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
    run: ../tools/kallisto_calc_expression.cwl
    in:
      transcript_idx: kallisto_idx
      strand: strand_parse/kallisto_std
      reads1: cutadapt/trimmedReadsR1
      reads2: cutadapt/trimmedReadsR2
      SampleID: sample_name
    out: [
      abundance_out
    ]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 3