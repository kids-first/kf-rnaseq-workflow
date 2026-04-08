cwlVersion: v1.2
class: Workflow
id: prepare_gencode_refs
label: GENCODE Ref Builder
doc: |
  # GENCODE Reference Builder

  ## Introduction

  GENCODE is a project that aims to identify and classify all gene features in the human genome.
  The tools of the Kids First RNAseq workflow makes extensive use of GENCODE-based references.
  These references provide the necessary information to perform read assignment, quanfitifaciton,
  and annoatation of RNAseq data.

  GENCODE frequently updates their annotations. Over the past five years, GENCODE has released a
  new version 2-3 times a year. While version-to-version changes tend to be minor, over the course
  of multiple years, the differences can become significant.

  For large projects, such as Kids First, it's important to have comparable data generated over large
  spans of time. To achieve this we lock our version of GENCODE and generate stable references that
  will be used for many years.

  The pipeline boils down the creation of the Kids First GENCODE-based references to a series of
  version declarations. The pipeline will download the appropriate reference file for each version
  provided and perform the necessary processes with the necessary flags to build the KF references.

  ### Runtime and Instance Recommendation

  If time is of the essence, we recommend using direct instances to limit retries of longer tasks.

  Expected runtime for tools:
  - Under 1 minute: `GTEX_COLLAPSE_ANNOTATION`, `T1K_BUILD`
  - Under 10 minutes: `DOWNLOAD_GENCODE`, `KALLISO_INDEX`
  - Under 2 hours: `RSEM_GENERATE_GENOME`, `STAR_GENERATE_GENOME`
  - 24+ hours: `STAR_FUSION_GENOME_GENERATE`

  `STAR_FUSION_GENOME_GENERATE` is an incredibly long running tool. Unfortunately, it seems like
  there's little motivation to improve that.

  ### Software Changes

  This pipeline is intended for changes of only the versions of the reference files. Changes in the
  software used to build the GENCODE-based references is not possible in CWL. Currently, the pipeline
  is using the versions of software used in the KF RNAseq pipeline. See the [dockers](./dockers_rnaseq.md).
  If updating the Dockers of this pipeline, make sure the command line parameters have not changed.

  ## Usage

  ### Inputs

  - `gencode_version`: Version of GENCODE to download from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/
  - `ctat_resource_version`: Version of CTAT resource SOURCE file to download from https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/
  - `ctat_fusion_version`: Version of CTAT fusion dat.gz to download from https://github.com/FusionAnnotator/CTAT_HumanFusionLib/releases/
  - `hla_version`: Version of HLA to download from https://github.com/ANHIG/IMGTHLA/releases/
  - `dfam_version`: Version of DFAM to download from https://dfam.org/releases/
  - `pfam_version`: Version of PFAM to download from https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/
  - `date_today`: Todays date in the following format: Mar012021 (MonthDayYear). Used to name output

  ### Outputs

  - `resource_manifest`: This file contains the URLs of all files downloaded. Mostly for recordkeeping purposes, it also ends up in the `star_fusion_genome` TAR
  - `gencode_genome`: GENCODE `*primary_assembly.fa` associated with the `gencode_version`
  - `gencode_annotation`: GENCODE `*primary_assembly.annotation.gtf.gz` associated with the `gencode_version`
  - `gtex_collapsed_annotation`: Collapsed (i.e., combining all isoforms of a gene into a single transcript) version of `gencode_annotation` created using https://github.com/broadinstitute/gtex-pipeline/tree/master/gene_model with the `--stranded` parameter set
  - `rsem_genome`: TAR of transcript references for RSEM created using the `gencode_genome` and `gencode_annotation`
  - `kallisto_idx`: IDX of the transcriptome FASTA within `rsem_genome` for Kallisto
  - `star_genome`: TAR of binary genome sequence, suffix arrays, text chromosome names/lengths, splice junctions coordinates, and transcripts/genes information made using STAR in `--runMode genomeGenerate` from `gencode_genome` and `gencode_annotation`
  - `star_fusion_genome`: STAR fusion genome TAR made from `gencode_genome`, `gencode_annotation`, `ctat_resource_version` SOURCE TAR + filtering rules, `ctat_fusion_version` TAR, `dfam_version` DB, and `pfam_version` DB
  - `star_fusion_annot`: Pared version of `star_fusion_genome` containing only `fusion_annot_lib.idx` and `blast_pairs.idx`. Only used to reduce I/O for fusion annotation steps.
  - `hla_rna_ref_seqs`: T1k HLA reference sequences made using `gencode_annotation` and `hla_version` release DAT file
  - `hla_rna_gene_coords`: T1k HLA reference coordinates made using `gencode_annotation` and `hla_version` release DAT file
requirements:
- class: ScatterFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
inputs:
  gencode_version: {type: 'string', doc: "Version of GENCODE to download from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"}
  ctat_resource_version: {type: 'string', doc: "Version of CTAT resource SOURCE file to download from https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/"}
  ctat_fusion_version: {type: 'string', doc: "Version of CTAT fusion dat.gz to download from https://github.com/FusionAnnotator/CTAT_HumanFusionLib/releases"}
  hla_version: {type: 'string', doc: "Version of HLA to download from https://github.com/ANHIG/IMGTHLA/releases"}
  dfam_version: {type: 'string', doc: "Version of DFAM to download from https://dfam.org/releases/"}
  pfam_version: {type: 'string', doc: "Version of PFAM to download from https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/"}
  date_today: {type: 'string', doc: "Todays date in the following format: Mar012021 (MonthDayYear)"}
outputs:
  resource_manifest: {type: 'File', outputSource: download_gencode/manifest}
  gencode_genome: {type: 'File', outputSource: download_gencode/gencode_genome}
  gencode_annotation: {type: 'File', outputSource: download_gencode/gencode_annotation}
  gtex_collapsed_annotation: {type: 'File', outputSource: gtex_collapse_annotation/collapsed_gtf}
  rsem_genome: {type: 'File', outputSource: rsem_generate_genome/genome_tar}
  kallisto_idx: {type: 'File', outputSource: kallisto_index/index}
  star_genome: {type: 'File', outputSource: star_genome_generate/star_ref}
  star_fusion_genome: {type: 'File', outputSource: star_fusion_genome_generate/star_fusion_reference}
  star_fusion_annot: {type: 'File', outputSource: generate_fusion_annot/fusion_annot_ref}
  hla_rna_ref_seqs: {type: 'File', outputSource: t1k_build/rna_ref_seqs}
  hla_rna_gene_coords: {type: 'File', outputSource: t1k_build/rna_gene_coords}
steps:
  download_gencode:
    run: ../tools/download_gencode.cwl
    in:
      gencode_version: gencode_version
      ctat_resource_version: ctat_resource_version
      ctat_fusion_version: ctat_fusion_version
      hla_version: hla_version
      dfam_version: dfam_version
      pfam_version: pfam_version
    out: [manifest, annot_filter_rule, gencode_genome, gencode_annotation, ctat_resource, ctat_fusion, hla, pfam, dfam]
  rsem_generate_genome:
    run: ../tools/rsem_prepare_reference.cwl
    in:
      gencode_version: gencode_version
      gtf: download_gencode/gencode_annotation
      reference: download_gencode/gencode_genome
      output_prefix:
        valueFrom: |
          $("RSEM_GENCODE" + inputs.gencode_version)
    out: [genome_tar, transcripts_fasta]
  gtex_collapse_annotation:
    run: ../tools/gtex_collapse_annotation.cwl
    in:
      gencode_version: gencode_version
      gtf: download_gencode/gencode_annotation
      output_filename:
        valueFrom: |
          $("gencode.v" + inputs.gencode_version + ".primary_assembly.rnaseqc.stranded.gtf")
      stranded:
        valueFrom: |
          $(true)
    out: [collapsed_gtf]
  star_genome_generate:
    run: ../tools/star_2.7.10a_genome_generate.cwl
    in:
      gencode_version: gencode_version
      genome_fa: download_gencode/gencode_genome
      gtf: download_gencode/gencode_annotation
      genomeDir:
        valueFrom: |
          $("STAR_2.7.10a_GENCODE" + inputs.gencode_version)
    out: [star_ref]
  kallisto_index:
    run: ../tools/kallisto_index.cwl
    in:
      gencode_version: gencode_version
      transcripts_fasta: rsem_generate_genome/transcripts_fasta
      output_filename:
        valueFrom: |
          $("RSEM_GENCODE" + inputs.gencode_version + ".transcripts.kallisto.idx")
    out: [index]
  t1k_build:
    run: ../tools/t1k_build.cwl
    in:
      gencode_version: gencode_version
      hla_version: hla_version
      dat: download_gencode/hla
      genome_annot: download_gencode/gencode_annotation
      prefix:
        valueFrom: |
          $("hla_v" + inputs.hla_version + "_gencode_v" + inputs.gencode_version)
    out: [dna_gene_coords, dna_ref_seqs, rna_gene_coords, rna_ref_seqs]
  star_fusion_genome_generate:
    run: ../tools/star_fusion_1.10.1_gen_reference.cwl
    in:
      gencode_version: gencode_version
      date: date_today
      resource_manifest: download_gencode/manifest
      ctat_source: download_gencode/ctat_resource
      genome_fa: download_gencode/gencode_genome
      reference_gtf: download_gencode/gencode_annotation
      fusion_annot_lib: download_gencode/ctat_fusion
      annot_filter_rule: download_gencode/annot_filter_rule
      pfam_db: download_gencode/pfam
      dfam_db: download_gencode/dfam
      output_dir:
        valueFrom: |
          $("GRCh38_v" + inputs.gencode_version + "_CTAT_lib_" + inputs.date + ".CUSTOM")
    out: [star_fusion_reference]
  generate_fusion_annot:
    run: ../tools/generate_fusion_annot.cwl
    in:
      fusion_tar: star_fusion_genome_generate/star_fusion_reference
      output_basename:
        valueFrom: |
          $(inputs.fusion_tar.basename.replace(/CUSTOM.tar.gz/, "fusion_annot"))
    out: [fusion_annot_ref]
hints:
- class: "sbg:maxNumberOfParallelInstances"
  value: 3
sbg:license: Apache License 2.0
sbg:publisher: KFDRC
$namespaces:
  sbg: https://sevenbridges.com
sbg:categories:
- FUSION
- GENCODE
- GTEX
- KALLISTO
- RSEM
- STAR
- STARFUSION
- T1K
