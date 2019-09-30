cwlVersion: v1.0
class: Workflow
doc: >-
  ![data service logo](../docs/dataservice.png) 
  
  
  This is the Kids First Data Resource Center RNA-Seq Workflow, which includes fusion and expression detection. 
  This workflow takes either BAM or FASTQ input and uses STAR for alignment to the reference. 
  STAR gives aligned genomic, chimeric, transcriptomic, and junction read outputs. 
  Arriba and STAR fusion mode are run for fusion estimations on STAR alignment chimeric output. 
  RSEM is run for gene expression estimations on STAR transcriptomic output.
  Additional gene abundance estimations are gathered from Kallisto which is run on trimmed reads of the raw read input.
  RNAseQC is run to generate a number of metrics including mapping rates, transcript counts, gene counts, sense and antisense mapping along with others. 
  
  
  ### Gene Expression Abundance Estimation:
  
  [STAR v2.6.1d](https://doi.org/f4h523) is used to align paired-end RNA-seq reads.
  This output is used for all subsequent RNA analysis. The reference we used, and is recommended for use, was that of ENSEMBL's [GENCODE 27](https://www.gencodegenes.org/human/release_27.html), "Comprehensive gene annotation."
  [RSEM v1.3.1](https://doi:10/cwg8n5) is used for transcript- and gene-level quantification.
  A second method of quantification was added using [Kallisto v0.43.1](https://doi:10.1038/nbt.3519).
  This method differs in that it uses pseudoaligments using fastq reads directly to the aforementioned GENCODE 27 reference.

  ### RNA Fusion Calling:
  
  [Arriba v1.1.0](https://github.com/suhrig/arriba/) and [STAR-Fusion 1.5.0](https://doi:10.1101/120295) fusion detection tools are set up for fusion calling.
  For both of these tools, aligned BAM and chimeric SAM files are used from STAR as inputs and `GRCh38_gencode_v27` GTF for gene annotation is recommended.
  STAR-Fusion is set with default parameters and we recommend to annotate all fusion calls with `GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz` provided in the STAR-fusion release.
  For Arriba, we used a blacklist file `blacklist_hg38_GRCh38_2018-11-04.tsv.gz` from the Arriba release tarballs to remove recurrent fusion artifacts and transcripts present in healthy tissue.
  
  ### Tips To Run:

  1) For fastq input, please enter the reads 1 file in reads1 and the reads 2 file in reads2. For bam input, please enter the reads file in reads1 and leave reads2 empty as it is optional.
  
  
  2) r1_adapter and r2_adapter are OPTIONAL.  If the input reads have already been trimmed, leave these as null and
  cutadapt step will simple pass on the fastq files to STAR.  If they do need trimming, supply the adapters and the cutadapt step will trim, and pass trimmed fastqs along
  
  
  3) `wf_strand_param` is a workflow convenience param so that, if you input the following, the equivalent will propagate
  to the four tools that use that parameter:
  
    - `default`: 'rsem_std': null, 'kallisto_std': null, 'rnaseqc_std': null, 'arriba_std': null. This means unstranded or auto in the case of arriba.
    - `rf_stranded`: 'rsem_std': 0, 'kallisto_std': 'rf-stranded', 'rnaseqc_std': 'rf', 'arriba_std': 'reverse'.  This means if read1 in the input fastq/bam is reverse complement to the transcript that it maps to.
    - `fr-stranded`: 'rsem_std': 1, 'kallisto_std': 'fr-stranded', 'rnaseqc_std': 'fr', 'arriba_std': 'yes'. This means if read1 in the input fastq/bam is the same sense (maps 5' to 3') to the transcript that it maps to.

  4) Suggested `STAR_outSAMattrRGline`, with **TABS SEPARATING THE TAGS**,  format is:

    - `ID:sample_name LB:aliquot_id   PL:platform SM:BSID`
    - for example: `ID:7316-242   LB:750189 PL:ILLUMINA SM:BS_W72364MN`

  5) Suggested inputs are:
  
    ```
    FusionGenome: GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz
    gtf_anno: gencode.v27.primary_assembly.annotation.gtf
    RNAseQC_GTF: gencode.v27.primary_assembly.RNAseQC.gtf
    RSEMgenome: RSEM_GENCODE27.tar.gz
    STARgenome: STAR_GENCODE27.tar.gz
    reference_fasta: GRCh38.primary_assembly.genome.fa
    kallisto_idx: gencode.v27.kallisto.index
    ```
    
  ### Links/Resources:
  
  The related Github branch for this app is located [here](https://github.com/kids-first/kf-rnaseq-workflow/tree/be-publish).
    
id: kfdrc-rnaseq-wf
label: Kids First DRC RNA-Seq Workflow
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  sample_name: {type: string, doc: "Sample name used for file base name of all outputs"}
  r1_adapter: {type: ['null', string], doc: "Optional input. If the input reads have already been trimmed, leave these as null. If they do need trimming, supply the adapters."}
  r2_adapter: {type: ['null', string], doc: "Optional input. If the input reads have already been trimmed, leave these as null. If they do need trimming, supply the adapters."}
  reads1: {type: File, doc: "For FASTQ input, please enter reads 1 here. For BAM input, please enter reads here."}
  reads2: {type: File?, doc: "For FASTQ input, please enter reads 2 here. For BAM input, leave empty."}
  STARgenome: {type: File, doc: "STAR_GENCODE27.tar.gz", sbg:suggestedValue: {class: 'File', path: '5d8bb21fe4b0950c4028f853', name: 'STAR_GENCODE27.tar.gz'}}
  RSEMgenome: {type: File, doc: "RSEM_GENCODE27.tar.gz", sbg:suggestedValue: {class: 'File', path: '5d8bb21fe4b0950c4028f851', name: 'RSEM_GENCODE27.tar.gz'}}
  reference_fasta: {type: File, doc: "GRCh38.primary_assembly.genome.fa", sbg:suggestedValue: {class: 'File', path: '5d8bb21fe4b0950c4028f855', name: 'GRCh38.primary_assembly.genome.fa'}}
  gtf_anno: {type: File, doc: "gencode.v27.primary_assembly.annotation.gtf", sbg:suggestedValue: {class: 'File', path: '5d8bb21fe4b0950c4028f84f', name: 'gencode.v27.primary_assembly.annotation.gtf'}}
  FusionGenome: {type: File, doc: "GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz", sbg:suggestedValue: {class: 'File', path: '5d8bb21fe4b0950c4028f854', name: 'GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz'}}
  runThread: {type: int, doc: "Amount of threads for analysis."}
  STAR_outSAMattrRGline: {type: string, doc: "Suggested setting, with TABS SEPARATING THE TAGS, format is: ID:sample_name LB:aliquot_id PL:platform SM:BSID for example ID:7316-242 LB:750189 PL:ILLUMINA SM:BS_W72364MN"}
  RNAseQC_GTF: {type: File, doc: "gencode.v27.primary_assembly.RNAseQC.gtf", sbg:suggestedValue: {class: 'File', path: '5d8bb21fe4b0950c4028f852', name: 'gencode.v27.primary_assembly.RNAseQC.gtf'}}
  kallisto_idx: {type: File, doc: "gencode.v27.kallisto.index", sbg:suggestedValue: {class: 'File', path: '5d8bb21fe4b0950c4028f850', name: 'gencode.v27.kallisto.index'}}
  wf_strand_param: {type: [{type: enum, name: wf_strand_param, symbols: ["default", "rf-stranded", "fr-stranded"]}], doc: "use 'default' for unstranded/auto, 'rf-stranded' if read1 in the fastq read pairs is reverse complement to the transcript, 'fr-stranded' if read1 same sense as transcript"}
  input_type: {type: [{type: enum, name: input_type, symbols: ["BAM", "FASTQ"]}], doc: "Please select one option for input file type, BAM or FASTQ."}

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

  bam2fastq:
    run: ../tools/samtools_fastq.cwl
    in:
      input_reads_1: reads1
      input_reads_2: reads2
      SampleID: sample_name
      runThreadN: runThread
      input_type: input_type
    out: [
      fq1,
      fq2
    ]

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
      wf_strand_param: wf_strand_param
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