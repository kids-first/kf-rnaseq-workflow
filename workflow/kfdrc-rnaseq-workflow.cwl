cwlVersion: v1.2
class: Workflow
id: kfdrc-rnaseq-workflow
label: Kids First DRC RNAseq Workflow
doc: |
  # Kids First RNA-Seq Workflow

  This is the Kids First RNA-Seq pipeline, which includes fusion and expression detection.

  ![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

  ## Introduction
  This pipeline utilizes cutadapt to trim adapters from the raw reads, if necessary, and passes the reads to STAR for alignment.
  The alignment output is used by RSEM for gene expression abundance estimation and rMATS for differential alternative splicing events detection.
  Additionally, Kallisto is used for quantification, but uses pseudoalignments to estimate the gene abundance from the raw data.
  Fusion calling is performed using Arriba and STAR-Fusion detection tools on the STAR alignment outputs.
  Filtering and prioritization of fusion calls is done by annoFuse.
  Metrics for the workflow are generated by RNA-SeQC.
  Junction files for the workflow are generated by rMATS.

  If you would like to run this workflow using the cavatica public app, a basic primer on running public apps can be found [here](https://www.notion.so/d3b/Starting-From-Scratch-Running-Cavatica-af5ebb78c38a4f3190e32e67b4ce12bb).
  Alternatively, if you'd like to run it locally using `cwltool`, a basic primer on that can be found [here](https://www.notion.so/d3b/Starting-From-Scratch-Running-CWLtool-b8dbbde2dc7742e4aff290b0a878344d) and combined with app-specific info from the readme below.
  This workflow is the current production workflow, equivalent to this [Cavatica public app](https://cavatica.sbgenomics.com/public/apps#cavatica/apps-publisher/kfdrc-rnaseq-workflow).

  ### Cutadapt
  [Cutadapt v2.5](https://github.com/marcelm/cutadapt) Cut adapter sequences from raw reads if needed.
  ### STAR
  [STAR v2.6.1d](https://doi.org/f4h523) RNA-Seq raw data alignment.
  ### RSEM
  [RSEM v1.3.1](https://doi:10/cwg8n5) Calculation of gene expression.
  ### Kallisto
  [Kallisto v0.43.1](https://doi:10.1038/nbt.3519) Raw data pseudoalignment to estimate gene abundance.
  ### STAR-Fusion
  [STAR-Fusion v1.5.0](https://doi:10.1101/120295) Fusion detection for `STAR` chimeric reads.
  ### Arriba
  [Arriba v1.1.0](https://github.com/suhrig/arriba/) Fusion caller that uses `STAR` aligned reads and chimeric reads output.
  ### annoFuse
  [annoFuse 0.90.0](https://github.com/d3b-center/annoFuse/releases/tag/v0.90.0) Filter and prioritize fusion calls. For more information, please see the following [paper](https://www.biorxiv.org/content/10.1101/839738v3).
  ### RNA-SeQC
  [RNA-SeQC v2.3.4](https://github.com/broadinstitute/rnaseqc) Generate metrics such as gene and transcript counts, sense/antisene mapping, mapping rates, etc
  ### rMATS
  [rMATS turbo v4.1.2](https://github.com/Xinglab/rmats-turbo) Computational tool to detect differential alternative splicing events from RNA-Seq data

  ## Usage

  ### Runtime Estimates:
  - 8 GB single end FASTQ input: 66 Minutes & $2.00
  - 17 GB single end FASTQ input: 58 Minutes & $2.00

  ### Inputs common:
  ```yaml
  inputs:
    sample_name: string
    r1_adapter: {type: ['null', string]}
    r2_adapter: {type: ['null', string]}
    STAR_outSAMattrRGline: string
    STARgenome: File
    RSEMgenome: File
    reference_fasta: File
    gtf_anno: File
    FusionGenome: File
    runThread: int
    RNAseQC_GTF: File
    kallisto_idx: File
    wf_strand_param: {type: [{type: enum, name: wf_strand_param, symbols: ["default", "rf-stranded", "fr-stranded"]}], doc: "use 'default' for unstranded/auto, 'rf-stranded' if read1 in the fastq read pairs is reverse complement to the transcript, 'fr-stranded' if read1 same sense as transcript"}
    input_type: {type: [{type: enum, name: input_type, symbols: ["BAM", "FASTQ"]}], doc: "Please select one option for input file type, BAM or FASTQ."}
    rmats_read_length: {type: 'int', doc: "Input read length for sample reads."}
    rmats_variable_read_length: {type: 'boolean?', doc: "Allow reads with lengths that differ from --readLength to be processed. --readLength will still be used to determine IncFormLen and SkipFormLen."}
    rmats_novel_splice_sites: {type: 'boolean?', doc: "Select for novel splice site detection or unannotated splice sites. 'true' to detect or add this parameter, 'false' to disable denovo detection. Tool Default: false"}
    rmats_stat_off: {type: 'boolean?', doc: "Select to skip statistical analysis, either between two groups or on single sample group. 'true' to add this parameter. Tool default: false"}
    rmats_allow_clipping: {type: 'boolean?', doc: "Allow alignments with soft or hard clipping to be used."}
    rmats_threads: {type: 'int?', doc: "Threads to allocate to RMATs."}
    rmats_ram: {type: 'int?', doc: "GB of RAM to allocate to RMATs."}
  ```

  ### Bam input-specific:
  ```yaml
  inputs:
    reads1: File
  ```

  ### PE Fastq input-specific:
  ```yaml
  inputs:
    reads1: File
    reads2: File
  ```

  ### SE Fastq input-specific:
  ```yaml
  inputs:
    reads1: File
  ```

  ### Run:

  1) For fastq or bam input, run `kfdrc-rnaseq-wf` as this can accept both file types.
  For PE fastq input, please enter the reads 1 file in `reads1` and the reads 2 file in `reads2`.
  For SE fastq input, enter the single ends reads file in `reads1` and leave `reads2` empty as it is optional.
  For bam input, please enter the reads file in `reads1` and leave `reads2` empty as it is optional.

  2) `r1_adapter` and `r2_adapter` are OPTIONAL.
  If the input reads have already been trimmed, leave these as null and cutadapt step will simple pass on the fastq files to STAR.
  If they do need trimming, supply the adapters and the cutadapt step will trim, and pass trimmed fastqs along.

  3) `wf_strand_param` is a workflow convenience param so that, if you input the following, the equivalent will propagate to the four tools that use that parameter:
      - `default`: 'rsem_std': null, 'kallisto_std': null, 'rnaseqc_std': null, 'arriba_std': null. This means unstranded or auto in the case of arriba.
      - `rf-stranded`: 'rsem_std': 0, 'kallisto_std': 'rf-stranded', 'rnaseqc_std': 'rf', 'arriba_std': 'reverse'.  This means if read1 in the input fastq/bam is reverse complement to the transcript that it maps to.
      - `fr-stranded`: 'rsem_std': 1, 'kallisto_std': 'fr-stranded', 'rnaseqc_std': 'fr', 'arriba_std': 'yes'. This means if read1 in the input fastq/bam is the same sense (maps 5' to 3') to the transcript that it maps to.

  4) Suggested `STAR_outSAMattrRGline`, with **TABS SEPARATING THE TAGS**,  format is:

      `ID:sample_name LB:aliquot_id   PL:platform SM:BSID` for example `ID:7316-242   LB:750189 PL:ILLUMINA SM:BS_W72364MN`
  5) Suggested inputs are:

      - `FusionGenome`: [GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz](https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.3/GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz)
      - `gtf_anno`: gencode.v27.primary_assembly.annotation.gtf, location: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.primary_assembly.annotation.gtf.gz, will need to unzip
      - `RNAseQC_GTF`: gencode.v27.primary_assembly.RNAseQC.gtf, built using `gtf_anno` and following build instructions [here](https://github.com/broadinstitute/rnaseqc#usage)
      - `RSEMgenome`: RSEM_GENCODE27.tar.gz, built using the `reference_fasta` and `gtf_anno`, following `GENCODE` instructions from [here](https://deweylab.github.io/RSEM/README.html), then creating a tar ball of the results.
      - `STARgenome`: STAR_GENCODE27.tar.gz, created using the star_genomegenerate.cwl tool, using the `reference_fasta`, `gtf_anno`, and setting `sjdbOverhang` to 100
      - `reference_fasta`: [GRCh38.primary_assembly.genome.fa](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh38.primary_assembly.genome.fa.gz), will need to unzip
      - `kallisto_idx`: gencode.v27.kallisto.index, built from gencode 27 trascript fasta: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.transcripts.fa.gz, following instructions from [here](https://pachterlab.github.io/kallisto/manual)

  6) rMATS requires you provide the length of the reads in the sample. If you are unsure of the length, you can set `rmats_variable_read_length` to true which will allow reads with a length other than the value you provided to be processed.

  ### Outputs:
  ```yaml
  outputs:
    cutadapt_stats: {type: 'File?', outputSource: cutadapt/cutadapt_stats, doc: "Cutadapt stats output, only if adapter is supplied."}
    STAR_transcriptome_bam: {type: 'File', outputSource: star/transcriptome_bam_out, doc: "STAR bam of transcriptome reads"}
    STAR_sorted_genomic_bam: {type: 'File', outputSource: samtools_sort/sorted_bam, doc: "STAR sorted alignment bam"}
    STAR_sorted_genomic_bai: {type: 'File', outputSource: samtools_sort/sorted_bai, doc: "STAR index for sorted aligned bam"}
    STAR_chimeric_bam_out: {type: 'File', outputSource: samtools_sort/chimeric_bam_out, doc: "STAR bam output of chimeric reads"}
    STAR_chimeric_junctions: {type: 'File', outputSource: star_fusion/chimeric_junction_compressed, doc: "STAR chimeric junctions"}
    STAR_gene_count: {type: 'File', outputSource: star/gene_counts, doc: "STAR gene counts"}
    STAR_junctions_out: {type: 'File', outputSource: star/junctions_out, doc: "STAR junction reads"}
    STAR_final_log: {type: 'File', outputSource: star/log_final_out, doc: "STAR metrics log file of unique, multi-mapping, unmapped, and chimeric reads"}
    STAR-Fusion_results: {type: 'File', outputSource: star_fusion/abridged_coding, doc: "STAR fusion detection from chimeric reads"}
    arriba_fusion_results: {type: 'File', outputSource: arriba_fusion/arriba_fusions, doc: "Fusion output from Arriba"}
    arriba_fusion_viz: {type: 'File', outputSource: arriba_fusion/arriba_pdf, doc: "pdf output from Arriba"}
    RSEM_isoform: {type: 'File', outputSource: rsem/isoform_out, doc: "RSEM isoform expression estimates"}
    RSEM_gene: {type: 'File', outputSource: rsem/gene_out, doc: "RSEM gene expression estimates"}
    RNASeQC_Metrics: {type: 'File', outputSource: rna_seqc/Metrics, doc: "Metrics on mapping, intronic, exonic rates, count information, etc"}
    RNASeQC_counts: {type: 'File', outputSource: supplemental/RNASeQC_counts, doc: "Contains gene tpm, gene read, and exon counts"}
    kallisto_Abundance: {type: 'File', outputSource: kallisto/abundance_out, doc: "Gene abundance output from STAR genomic bam file"}
    annofuse_filtered_fusions_tsv: {type: 'File?', outputSource: annofuse/annofuse_filtered_fusions_tsv, doc: "Filtered fusions called by annoFuse."}
    rmats_filtered_alternative_3_prime_splice_sites_jc: {type: 'File', outputSource: rmats/filtered_alternative_3_prime_splice_sites_jc, doc: "Alternative 3 prime splice sites JC.txt output from RMATs containing only those calls with 10 or more read counts of support" }
    rmats_filtered_alternative_5_prime_splice_sites_jc: {type: 'File', outputSource: rmats/filtered_alternative_5_prime_splice_sites_jc, doc: "Alternative 5 prime splice sites JC.txt output from RMATs containing only those calls with 10 or more read counts of support" }
    rmats_filtered_mutually_exclusive_exons_jc: {type: 'File', outputSource: rmats/filtered_mutually_exclusive_exons_jc, doc: "Mutually exclusive exons JC.txt output from RMATs containing only those calls with 10 or more read counts of support" }
    rmats_filtered_retained_introns_jc: {type: 'File', outputSource: rmats/filtered_retained_introns_jc, doc: "Retained introns JC.txt output from RMATs containing only those calls with 10 or more read counts of support" }
    rmats_filtered_skipped_exons_jc: {type: 'File', outputSource: rmats/filtered_skipped_exons_jc, doc: "Skipped exons JC.txt output from RMATs containing only those calls with 10 or more read counts of support" }
  ```

  ![pipeline flowchart](https://github.com/kids-first/kf-rnaseq-workflow/blob/master/docs/kfdrc-rnaseq-workflow.png?raw=true)

  # D3b annoFuse Workflow

  ## Introduction

  In this workflow, annoFuse performs standardization of StarFusion and arriba output files to retain information regarding fused genes, breakpoints, reading frame information as well as annotation from FusionAnnotator, output format description [here](https://github.com/d3b-center/annoFuse/wiki#1-standardize-calls-from-fusion-callers-to-retain-information-regarding-fused-genesbreakpoints-reading-frame-information-as-well-as-annotation-from-fusionannotator). Basic artifact filtering to remove fusions among gene paralogs, conjoined genes and fused genes found in normal samples is also performed by filtering fusions annotated by [FusionAnnotator](https://github.com/d3b-center/FusionAnnotator) with "GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap|ConjoinG". Each fusion call needs at least one junction reads support to be retained as true call. Additionally, if a fusion call has large number of spanning fragment reads compared to junction reads (spanning fragment minus junction read greater than ten), we remove these calls as potential false positives. An expression based filter is also applied, requiring a min FPKM value of 1 for the fusion genes in question.

  Please refer to [annoFuse](https://github.com/d3b-center/annoFuse) R package for additional applications like putative oncogene annotations.

  ## Usage

  ### Inputs

  ```yaml
  inputs:
    sample_name: { type: 'string', doc: "Sample name used for file base name of all outputs" }
    FusionGenome: { type: 'File', doc: "GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz", sbg:suggestedValue: { class: 'File', path: '5d9c8d04e4b0950cce147f94', name: 'GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz' }}
    genome_untar_path: { type: 'string?', doc: "This is what the path will be when genome_tar is unpackaged", default: "GRCh38_v27_CTAT_lib_Feb092018/ctat_genome_lib_build_dir" }
    rsem_expr_file: { type: 'File', doc: "gzipped rsem gene expression file" }
    arriba_output_file: { type: 'File', doc: "Output from arriba, usually extension arriba.fusions.tsv" }
    col_num: { type: 'int?', doc: "column number in file of fusion name." }
    star_fusion_output_file: { type: 'File', doc: "Output from arriba, usually extension STAR.fusion_predictions.abridged.coding_effect.tsv" }
    output_basename: { type: 'string', doc: "String to use as basename for outputs" }
  ```

  ### Run

  1) Outputs from the arriba and STAR Fusion runs are required ahead of time (main RNAseq worflow output)
  2) Gzipped rsem counts file, also generated in main RNAseq workflow
  3) `FusionGenome` should match what was used to run STAR Fusion

  ### Outputs

  ```yaml
  outputs:
    annofuse_filtered_fusions_tsv: {type: File, outputSource: annoFuse_filter/filtered_fusions_tsv, doc: "Filtred output of formatted and annotated Star Fusion and arriba results"}
  ```

  # D3b rMATS Workflow

  ## Introduction

  The rMATS workflow can also be run as a standalone workflow. In this workflow, rMATS is run on the input BAM files to generate 5 junction files: `[alternative_3_prime_splice_sites_jc, alternative_5_prime_splice_sites_jc, mutually_exclusive_exons_jc, retained_introns_jc, skipped_exons_jc]`. The workflow next grabs the sample information from the `sample_1_bams` by parsing the read group information from the BAM header for use in the output names. Each of the five junction files then undergo a simple filtering process where calls that have junction counts less than 10 are removed. These filtered junction files are returned as the final outputs.

  ## Usage

  ### Inputs

  ```yaml
  inputs:
    gtf_annotation: { type: 'File', doc: "Input gtf annotation file." }
    sample_1_bams: { type: 'File[]', doc: "Input sample 1 bam files." }
    sample_2_bams: { type: 'File[]?', doc: "Input sample 2 bam files." }
    read_length: { type: 'int', doc: "Input read length for sample reads." }
    variable_read_length: { type: 'boolean?', doc: "Allow reads with lengths that differ from --readLength to be processed. --readLength will still be used to determine IncFormLen and SkipFormLen." }
    read_type:
      type:
      - "null"
      - type: enum
        symbols:
        - paired
        - single
        name: read_type
      doc: "Select one option for input read type either paired or single. Tool default: paired"
    strandedness:
      type:
      - "null"
      - type: enum
        symbols:
          - fr-unstranded
          - fr-firststrand
          - fr-secondstrand
        name: strandedness
      doc: "Select one option for input strandedness. Tool default: fr-unstranded"
    novel_splice_sites: { type: 'boolean?', doc: "Select for novel splice site detection or unannotated splice sites. 'true' to detect or add this parameter, 'false' to disable denovo detection. Tool Default: false" }
    stat_off: { type: 'boolean?', doc: "Select to skip statistical analysis, either between two groups or on single sample group. 'true' to add this parameter. Tool default: false" }
    allow_clipping: { type: 'boolean?', doc: "Allow alignments with soft or hard clipping to be used." }
    output_basename: { type: 'string', doc: "String to use as basename for output files" }
    rmats_threads: { type: 'int?', doc: "Threads to allocate to RMATs." }
    rmats_ram: { type: 'int?', doc: "GB of RAM to allocate to RMATs." }
  ```

  ### Outputs

  ```yaml
  outputs:
    filtered_alternative_3_prime_splice_sites_jc: {type: 'File', outputSource: filter_alt_3_prime/output, doc: "Alternative 3 prime splice sites JC.txt output from RMATs containing only those calls with 10 or more read counts of support" }
    filtered_alternative_5_prime_splice_sites_jc: {type: 'File', outputSource: filter_alt_5_prime/output, doc: "Alternative 5 prime splice sites JC.txt output from RMATs containing only those calls with 10 or more read counts of support" }
    filtered_mutually_exclusive_exons_jc: {type: 'File', outputSource: filter_me_exons/output, doc: "Mutually exclusive exons JC.txt output from RMATs containing only those calls with 10 or more read counts of support" }
    filtered_retained_introns_jc: {type: 'File', outputSource: filter_retained_introns/output, doc: "Retained introns JC.txt output from RMATs containing only those calls with 10 or more read counts of support" }
    filtered_skipped_exons_jc: {type: 'File', outputSource: filter_skipped_exons/output, doc: "Skipped exons JC.txt output from RMATs containing only those calls with 10 or more read counts of support" }
  ```

requirements:
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: SubworkflowFeatureRequirement

inputs:
  sample_name: {type: 'string', doc: "Sample ID of the input reads"}
  output_basename: {type: 'string', doc: "String to use as basename for outputs"}
  r1_adapter: {type: 'string?', doc: "Optional input. If the input reads have already\
      \ been trimmed, leave these as null. If they do need trimming, supply the adapters."}
  r2_adapter: {type: 'string?', doc: "Optional input. If the input reads have already\
      \ been trimmed, leave these as null. If they do need trimming, supply the adapters."}
  reads1: {type: 'File', doc: "For FASTQ input, please enter reads 1 here. For BAM\
      \ input, please enter reads here."}
  reads2: {type: 'File?', doc: "For FASTQ input, please enter reads 2 here. For BAM\
      \ input, leave empty."}
  STARgenome: {type: 'File', doc: "STAR_GENCODE27.tar.gz", "sbg:suggestedValue": {
      class: File, path: 5f5001a6e4b054958bc8d2ec, name: STAR_GENCODE27.tar.gz}}
  RSEMgenome: {type: 'File', doc: "RSEM_GENCODE27.tar.gz", "sbg:suggestedValue": {
      class: File, path: 5f500135e4b0370371c051be, name: RSEM_GENCODE27.tar.gz}}
  reference_fasta: {type: 'File', doc: "GRCh38.primary_assembly.genome.fa", "sbg:suggestedValue": {
      class: File, path: 5f500135e4b0370371c051b4, name: GRCh38.primary_assembly.genome.fa}}
  gtf_anno: {type: 'File', doc: "gencode.v27.primary_assembly.annotation.gtf", "sbg:suggestedValue": {
      class: File, path: 5f500135e4b0370371c051c3, name: gencode.v27.primary_assembly.annotation.gtf}}
  FusionGenome: {type: 'File', doc: "GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz",
    "sbg:suggestedValue": {class: File, path: 5f500135e4b0370371c051b0, name: GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz}}
  runThread: {type: 'int', doc: "Amount of threads for analysis."}
  STAR_outSAMattrRGline: {type: 'string', doc: "Suggested setting, with TABS SEPARATING\
      \ THE TAGS, format is: ID:sample_name LB:aliquot_id PL:platform SM:BSID for\
      \ example ID:7316-242 LB:750189 PL:ILLUMINA SM:BS_W72364MN"}
  RNAseQC_GTF: {type: 'File', doc: "gencode.v27.primary_assembly.RNAseQC.gtf", "sbg:suggestedValue": {
      class: File, path: 5f500135e4b0370371c051c8, name: gencode.v27.primary_assembly.RNAseQC.gtf}}
  kallisto_idx: {type: 'File', doc: "gencode.v27.kallisto.index", "sbg:suggestedValue": {
      class: File, path: 5f500135e4b0370371c051bd, name: gencode.v27.kallisto.index}}
  wf_strand_param: {type: [{type: 'enum', name: wf_strand_param, symbols: ["default",
          "rf-stranded", "fr-stranded"]}], doc: "use 'default' for unstranded/auto,\
      \ 'rf-stranded' if read1 in the fastq read pairs is reverse complement to the\
      \ transcript, 'fr-stranded' if read1 same sense as transcript"}
  input_type: {type: [{type: 'enum', name: input_type, symbols: ["PEBAM", "SEBAM",
          "FASTQ"]}], doc: "Please select one option for input file type, PEBAM (paired-end\
      \ BAM), SEBAM (single-end BAM) or FASTQ."}
  kallisto_avg_frag_len: {type: 'int?', doc: "Optional input. Average fragment length\
      \ for Kallisto only if single end input."}
  kallisto_std_dev: {type: 'long?', doc: "Optional input. Standard Deviation of the\
      \ average fragment length for Kallisto only needed if single end input."}
  annofuse_genome_untar_path: {type: 'string?', doc: "This is what the path will be\
      \ when genome_tar is unpackaged", default: "GRCh38_v27_CTAT_lib_Feb092018/ctat_genome_lib_build_dir"}
  annofuse_col_num: {type: 'int?', doc: "column number in file of fusion name."}
  rmats_read_length: {type: 'int', doc: "Input read length for sample reads."}
  rmats_variable_read_length: {type: 'boolean?', doc: "Allow reads with lengths that\
      \ differ from --readLength to be processed. --readLength will still be used\
      \ to determine IncFormLen and SkipFormLen."}
  rmats_novel_splice_sites: {type: 'boolean?', doc: "Select for novel splice site\
      \ detection or unannotated splice sites. 'true' to detect or add this parameter,\
      \ 'false' to disable denovo detection. Tool Default: false"}
  rmats_stat_off: {type: 'boolean?', doc: "Select to skip statistical analysis, either\
      \ between two groups or on single sample group. 'true' to add this parameter.\
      \ Tool default: false"}
  rmats_allow_clipping: {type: 'boolean?', doc: "Allow alignments with soft or hard\
      \ clipping to be used."}
  rmats_threads: {type: 'int?', doc: "Threads to allocate to RMATs."}
  rmats_ram: {type: 'int?', doc: "GB of RAM to allocate to RMATs."}

outputs:
  cutadapt_stats: {type: 'File?', outputSource: cutadapt/cutadapt_stats, doc: "Cutadapt\
      \ stats output, only if adapter is supplied."}
  STAR_transcriptome_bam: {type: 'File', outputSource: star/transcriptome_bam_out,
    doc: "STAR bam of transcriptome reads"}
  STAR_sorted_genomic_bam: {type: 'File', outputSource: samtools_sort/sorted_bam,
    doc: "STAR sorted alignment bam"}
  STAR_sorted_genomic_bai: {type: 'File', outputSource: samtools_sort/sorted_bai,
    doc: "STAR index for sorted aligned bam"}
  STAR_chimeric_bam_out: {type: 'File', outputSource: samtools_sort/chimeric_bam_out,
    doc: "STAR bam output of chimeric reads"}
  STAR_chimeric_junctions: {type: 'File', outputSource: star_fusion/chimeric_junction_compressed,
    doc: "STAR chimeric junctions"}
  STAR_gene_count: {type: 'File', outputSource: star/gene_counts, doc: "STAR gene\
      \ counts"}
  STAR_junctions_out: {type: 'File', outputSource: star/junctions_out, doc: "STAR\
      \ junction reads"}
  STAR_final_log: {type: 'File', outputSource: star/log_final_out, doc: "STAR metrics\
      \ log file of unique, multi-mapping, unmapped, and chimeric reads"}
  STAR-Fusion_results: {type: 'File', outputSource: star_fusion/abridged_coding, doc: "STAR\
      \ fusion detection from chimeric reads"}
  arriba_fusion_results: {type: 'File', outputSource: arriba_fusion/arriba_fusions,
    doc: "Fusion output from Arriba"}
  arriba_fusion_viz: {type: 'File', outputSource: arriba_fusion/arriba_pdf, doc: "pdf\
      \ output from Arriba"}
  RSEM_isoform: {type: 'File', outputSource: rsem/isoform_out, doc: "RSEM isoform\
      \ expression estimates"}
  RSEM_gene: {type: 'File', outputSource: rsem/gene_out, doc: "RSEM gene expression\
      \ estimates"}
  RNASeQC_Metrics: {type: 'File', outputSource: rna_seqc/Metrics, doc: "Metrics on\
      \ mapping, intronic, exonic rates, count information, etc"}
  RNASeQC_counts: {type: 'File', outputSource: supplemental/RNASeQC_counts, doc: "Contains\
      \ gene tpm, gene read, and exon counts"}
  kallisto_Abundance: {type: 'File', outputSource: kallisto/abundance_out, doc: "Gene\
      \ abundance output from STAR genomic bam file"}
  annofuse_filtered_fusions_tsv: {type: 'File?', outputSource: annofuse/annofuse_filtered_fusions_tsv,
    doc: "Filtered fusions called by annoFuse."}
  rmats_filtered_alternative_3_prime_splice_sites_jc: {type: 'File', outputSource: rmats/filtered_alternative_3_prime_splice_sites_jc,
    doc: "Alternative 3 prime splice sites JC.txt output from RMATs containing only\
      \ those calls with 10 or more read counts of support"}
  rmats_filtered_alternative_5_prime_splice_sites_jc: {type: 'File', outputSource: rmats/filtered_alternative_5_prime_splice_sites_jc,
    doc: "Alternative 5 prime splice sites JC.txt output from RMATs containing only\
      \ those calls with 10 or more read counts of support"}
  rmats_filtered_mutually_exclusive_exons_jc: {type: 'File', outputSource: rmats/filtered_mutually_exclusive_exons_jc,
    doc: "Mutually exclusive exons JC.txt output from RMATs containing only those\
      \ calls with 10 or more read counts of support"}
  rmats_filtered_retained_introns_jc: {type: 'File', outputSource: rmats/filtered_retained_introns_jc,
    doc: "Retained introns JC.txt output from RMATs containing only those calls with\
      \ 10 or more read counts of support"}
  rmats_filtered_skipped_exons_jc: {type: 'File', outputSource: rmats/filtered_skipped_exons_jc,
    doc: "Skipped exons JC.txt output from RMATs containing only those calls with\
      \ 10 or more read counts of support"}

steps:

  bam2fastq:
    run: ../tools/samtools_fastq.cwl
    in:
      input_reads_1: reads1
      input_reads_2: reads2
      SampleID: output_basename
      runThreadN: runThread
      input_type: input_type
    out: [fq1, fq2]

  cutadapt:
    run: ../tools/cutadapter.cwl
    in:
      readFilesIn1: bam2fastq/fq1
      readFilesIn2: bam2fastq/fq2
      r1_adapter: r1_adapter
      r2_adapter: r2_adapter
      sample_name: output_basename
    out: [trimmedReadsR1, trimmedReadsR2, cutadapt_stats]

  star:
    run: ../tools/star_align.cwl
    in:
      outSAMattrRGline: STAR_outSAMattrRGline
      readFilesIn1: cutadapt/trimmedReadsR1
      readFilesIn2: cutadapt/trimmedReadsR2
      genomeDir: STARgenome
      runThreadN: runThread
      outFileNamePrefix: output_basename
    out: [chimeric_junctions, chimeric_sam_out, gene_counts, genomic_bam_out, junctions_out,
      log_final_out, log_out, log_progress_out, transcriptome_bam_out]

  samtools_sort:
    run: ../tools/samtools_sort.cwl
    in:
      unsorted_bam: star/genomic_bam_out
      chimeric_sam_out: star/chimeric_sam_out
    out: [sorted_bam, sorted_bai, chimeric_bam_out]

  rmats:
    run: ../workflow/rmats_wf.cwl
    in:
      gtf_annotation: gtf_anno
      sample_1_bams:
        source: samtools_sort/sorted_bam
        valueFrom: |
          $([self])
      read_length: rmats_read_length
      variable_read_length: rmats_variable_read_length
      read_type:
        source: bam2fastq/fq2
        valueFrom: |
          $(self == null ? "single" : "paired")
      strandedness:
        source: wf_strand_param
        valueFrom: |
          $(self == "rf-stranded" ? "fr-firststrand" : self == "fr-stranded" ? "fr-secondstrand" : "fr-unstranded")
      novel_splice_sites: rmats_novel_splice_sites
      stat_off: rmats_stat_off
      allow_clipping: rmats_allow_clipping
      output_basename: output_basename
      rmats_threads: rmats_threads
      rmats_ram: rmats_ram
    out: [filtered_alternative_3_prime_splice_sites_jc, filtered_alternative_5_prime_splice_sites_jc,
      filtered_mutually_exclusive_exons_jc, filtered_retained_introns_jc, filtered_skipped_exons_jc]

  strand_parse:
    run: ../tools/expression_parse_strand_param.cwl
    in:
      wf_strand_param: wf_strand_param
    out: [rsem_std, kallisto_std, rnaseqc_std, arriba_std]

  star_fusion:
    run: ../tools/star_fusion.cwl
    in:
      Chimeric_junction: star/chimeric_junctions
      genome_tar: FusionGenome
      SampleID: output_basename
    out: [abridged_coding, chimeric_junction_compressed]

  arriba_fusion:
    run: ../tools/arriba_fusion.cwl
    in:
      genome_aligned_bam: samtools_sort/sorted_bam
      genome_aligned_bai: samtools_sort/sorted_bai
      chimeric_sam_out: star/chimeric_sam_out
      reference_fasta: reference_fasta
      gtf_anno: gtf_anno
      outFileNamePrefix: output_basename
      arriba_strand_flag: strand_parse/arriba_std
    out: [arriba_fusions, arriba_pdf]

  rsem:
    run: ../tools/rsem_calc_expression.cwl
    in:
      bam: star/transcriptome_bam_out
      input_reads_2: cutadapt/trimmedReadsR2
      genomeDir: RSEMgenome
      outFileNamePrefix: output_basename
      strandedness: strand_parse/rsem_std
    out: [gene_out, isoform_out]

  rna_seqc:
    run: ../tools/RNAseQC.cwl
    in:
      Aligned_sorted_bam: samtools_sort/sorted_bam
      collapsed_gtf: RNAseQC_GTF
      strand: strand_parse/rnaseqc_std
      input_reads2: cutadapt/trimmedReadsR2
    out: [Metrics, Gene_TPM, Gene_count, Exon_count]

  supplemental:
    run: ../tools/supplemental_tar_gz.cwl
    in:
      outFileNamePrefix: output_basename
      Gene_TPM: rna_seqc/Gene_TPM
      Gene_count: rna_seqc/Gene_count
      Exon_count: rna_seqc/Exon_count
    out: [RNASeQC_counts]

  kallisto:
    run: ../tools/kallisto_calc_expression.cwl
    in:
      transcript_idx: kallisto_idx
      strand: strand_parse/kallisto_std
      reads1: cutadapt/trimmedReadsR1
      reads2: cutadapt/trimmedReadsR2
      SampleID: output_basename
      avg_frag_len: kallisto_avg_frag_len
      std_dev: kallisto_std_dev
    out: [abundance_out]

  annofuse:
    run: ../workflow/kfdrc_annoFuse_wf.cwl
    in:
      sample_name: sample_name
      FusionGenome: FusionGenome
      genome_untar_path: annofuse_genome_untar_path
      rsem_expr_file: rsem/gene_out
      arriba_output_file: arriba_fusion/arriba_fusions
      star_fusion_output_file: star_fusion/abridged_coding
      col_num: annofuse_col_num
      output_basename: output_basename
    out: [annofuse_filtered_fusions_tsv]

$namespaces:
  sbg: https://sevenbridges.com
hints:
- class: "sbg:maxNumberOfParallelInstances"
  value: 3
"sbg:license": Apache License 2.0
"sbg:publisher": KFDRC
"sbg:categories":
- ALIGNMENT
- ANNOFUSE
- ARRIBA
- BAM
- CUTADAPT
- FASTQ
- KALLISTO
- PE
- RNASEQ
- RNASEQC
- RSEM
- SE
- STAR
"sbg:links":
- id: 'https://github.com/kids-first/kf-rnaseq-workflow/releases/tag/v3.0.1'
  label: github-release
