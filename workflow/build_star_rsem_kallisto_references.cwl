cwlVersion: v1.2
class: Workflow
id: build_star_rsem_kallisto_references
label: "build_star_rsem_kallisto_references"
doc: |
  Build references for STAR, RSEM, and Kallisto for GRCm39

requirements:
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs:
  reference_fasta:
    type: File
    doc: "Reference fasta file"

  reference_name:
    type: string
    doc: "Output file prefix. Recommend format: RSEM_<SOURCE><Version>/"

  reference_gtf:
    type: File?
    doc: "Gene model definitions. This OR GFF required"

  reference_gff:
    type: File?
    doc: "Gene model definitions. This OR GTF required"

  transcript_idx:
    type: File
    doc: "Kallisto index file name"

outputs:
  rsem_reference_file:
    type: File
    outputSource: rsem_prepare_reference/rsem_reference

  kallisto_index_output:
    type: File
    outputSource: kallisto_index/index_out

  star_index_output:
    type: Directory
    outputSource: star_index/star_index

steps:
  rsem_prepare_reference:
    run: ../tools/rsem_prepare_reference.cwl
    in:
      reference_fasta: reference_fasta
      reference_name: reference_name
      reference_gff: reference_gff
      reference_gtf: reference_gtf
    out: [rsem_reference, rsem_fasta]

  kallisto_index:
    run: ../tools/kallisto_index.cwl
    in:
      transcript_idx: transcript_idx
      transcript_fasta: rsem_prepare_reference/rsem_fasta
    out: [index_out]

  star_index:
    run: ../tools/star_2.7.10a_genome_generate.cwl
    in:
      genomeFasta: reference_fasta
      annotationGTF: reference_gtf
      sjdbOverhang: { default: 100 }
      outDir: { default: "star_index_GRCm39" }
      runThreadN: { default: 8 }
    out: [star_index]

