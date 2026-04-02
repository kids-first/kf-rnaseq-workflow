cwlVersion: v1.2
class: Workflow
id: prepare_gencode_refs
requirements:
- class: ScatterFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
inputs:
  gencode_version: { type: 'string', doc: "Version of GENCODE to download from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/" }
  ctat_resource_version: { type: 'string', doc: "Version of CTAT resource SOURCE file to download from https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/" }
  ctat_fusion_version: { type: 'string', doc: "Version of CTAT fusion dat.gz to download from https://github.com/FusionAnnotator/CTAT_HumanFusionLib/releases" }
  hla_version: { type: 'string', doc: "Version of HLA to download from https://github.com/ANHIG/IMGTHLA/releases" }
  dfam_version: { type: 'string', doc: "Version of DFAM to download from https://dfam.org/releases/" }
  pfam_version: { type: 'string', doc: "Version of PFAM to download from https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/" }
  date_today: { type: 'string', doc: "Todays date in the following format: Mar012021 (MonthDayYear)" }
outputs:
  resource_manifest: { type: 'File', outputSource: download_gencode/manifest }
  gencode_genome: { type: 'File', outputSource: download_gencode/gencode_genome }
  gencode_annotation: { type: 'File', outputSource: download_gencode/gencode_annotation }
  gtex_collapsed_annotation: { type: 'File', outputSource: gtex_collapse_annotation/collapsed_gtf }
  rsem_genome: { type: 'File', outputSource: rsem_generate_genome/genome_tar }
  kallisto_idx: { type: 'File', outputSource: kallisto_index/index }
  star_genome: { type: 'File', outputSource: star_genome_generate/star_ref }
  star_fusion_genome: { type: 'File', outputSource: star_fusion_genome_generate/star_fusion_reference }
  star_fusion_annot: { type: 'File', outputSource: generate_fusion_annot/fusion_annot_ref }
  hla_rna_ref_seqs: {type: 'File', outputSource: t1k_build/rna_ref_seqs }
  hla_rna_gene_coords: {type: 'File', outputSource: t1k_build/rna_gene_coords }
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
    out: [ manifest, annot_filter_rule, gencode_genome, gencode_transcripts, gencode_annotation, ctat_resource, ctat_fusion, hla, pfam, dfam ]
  rsem_generate_genome:
    run: ../tools/rsem_prepare_reference.cwl
    in:
      gencode_version: gencode_version
      gtf: download_gencode/gencode_annotation
      reference: download_gencode/gencode_genome
      output_prefix:
        valueFrom: |
          $("RSEM_GENCODE" + inputs.gencode_version)
    out: [ genome_tar ]
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
    out: [ collapsed_gtf ]
  star_genome_generate:
    run: ../tools/star_2.7.10a_genome_generate.cwl
    in:
      gencode_version: gencode_version
      genome_fa: download_gencode/gencode_genome
      gtf: download_gencode/gencode_annotation
      genomeDir:
        valueFrom: |
          $("STAR_2.7.10a_GENCODE" + inputs.gencode_version)
    out: [ star_ref ]
  kallisto_index:
    run: ../tools/kallisto_index.cwl
    in:
      gencode_version: gencode_version
      transcripts_fasta: download_gencode/gencode_transcripts
      output_filename:
        valueFrom: |
          $("RSEM_GENCODE" + inputs.gencode_version + ".transcripts.kallisto.idx")
    out: [ index ]
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
    out: [ star_fusion_reference ]
  generate_fusion_annot:
    run: ../tools/generate_fusion_annot.cwl
    in:
      fusion_tar: star_fusion_genome_generate/star_fusion_reference
      output_basename:
        valueFrom: |
          $(inputs.fusion_tar.basename.replace(/CUSTOM.tar.gz/, "fusion_annot"))
    out: [ fusion_annot_ref ]
hints:
- class: "sbg:maxNumberOfParallelInstances"
  value: 3
sbg:license: Apache License 2.0
sbg:publisher: KFDRC
$namespaces:
  sbg: https://sevenbridges.com
