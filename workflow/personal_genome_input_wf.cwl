cwlVersion: v1.2
class: Workflow
id: personal-genome-input-wf
label: Personal Genome Input Workflow
doc: "Filter an input VCF and remove annotations"

inputs:
  # Strip, subset and PASS vars
  input_vcf: {type: File, secondaryFiles: ['.tbi']}
  strip_info: {type: 'string?', doc: "If given, remove previous annotation information based on INFO file, i.e. to strip VEP info, use INFO/ANN"}
  output_basename: string
  include_expression: { type: 'string?', doc: "See bcftools docs for valid expression. Can't be used at the same time as exclude_expression"}
  exclude_expression: {type: 'string?', doc: "See bcftools docs for valid expression. Can't be used at the same time as include_expression"}
  sample_name: { type: 'string?', doc: "csv string of samples if user wishes to apply filtering to and output specific samples"}
  filter_type: { type: 'string?', doc: "Apply a FILTER value expression", default: "PASS"}
  subtract_bed: {type: 'File?', doc: "Supply if you want to remove regions for any reason, like low complexity or repeat mask, etc" }
  # Genome gen vars
  genome_fa: { type: File, doc: "Fasta file to index. Recommend from GENCODE, PRI assembly. Must unzip first if compressed" }
  genomeTransformType: { type: [ 'null', {type: enum, name: genomeTransformType, symbols: [
      "None",
      "Haploid",
      "Diploid"
      ]}],
  doc: "type of genome transformation - None: no transformation. Haploid: replace reference alleles with alternative alleles from VCF file (e.g. consensus allele) \
  Diploid: create two haplotypes for each chromosome listed in VCF file, for genotypes 1—2, assumes perfect phasing (e.g. personal genome)" }
  gtf: { type: File, doc: "Matched GTF file to index. Recommend from GENCODE, PRI assembly" }
  runThreadN: { type: 'int?', default: 16 }
  memory: { type: 'int?', doc: "Mem in GB required. With no VCF, 60DB is fine, need more with VCF", default: 60 }
  sjdbOverhang: { type: 'int?', default: 100, doc: "Ideal value is read len minus 1, but default 100 ok for most cases" }
outputs:
  star_ref: { type: File, outputSource: star_personal_genome_generate/star_ref }
  debug_log: { type: File, outputSource: star_personal_genome_generate/debug_log }
steps:
  # Really just here to limit downstream uncompressed size
  bcftools_strip_info:
    run: ../tools/bcftools_strip_ann.cwl
    when: $(inputs.strip_info != null)
    in:
      input_vcf: input_vcf
      output_basename: output_basename
      tool_name:
        valueFrom: "dna_in"
      strip_info: strip_info
    out: [stripped_vcf]
  bedtools_subtract:
    run: ../tools/bedtools_subtract.cwl
    when: $(inputs.subtract_bed != null)
    in:
      input_vcf: 
        source: [bcftools_strip_info/stripped_vcf, input_vcf]
        pickValue: first_non_null
      subtract_bed: subtract_bed
      output_basename: output_basename
    out: [subtracted_vcf]
  bcftools_subset_vcf:
    run: ../tools/bcftools_filter_vcf.cwl
    when: $(inputs.include_expression != null || inputs.exclude_expression != null || inputs.filter_type != null)
    in:
      input_vcf: 
        source: [bedtools_subtract/subtracted_vcf, bcftools_strip_info/stripped_vcf, input_vcf]
        pickValue: first_non_null
      include_expression: include_expression
      exclude_expression: exclude_expression
      sample_name: sample_name
      filter_type: filter_type
      output_type:
        valueFrom: "v"
      output_basename: output_basename
    out: [filtered_vcf]
  star_personal_genome_generate:
    run: ../tools/star_2.7.11b_personal_genome_generate.cwl
    in:
      genomeDir:
        source: output_basename
        valueFrom: $(self + ".STAR_2.7.11b_diploid_genome")
      genome_fa: genome_fa
      genomeTransformVCF: bcftools_subset_vcf/filtered_vcf
      genomeTransformType: genomeTransformType
      gtf: gtf
      runThreadN: runThreadN
      memory: memory
      sjdbOverhang: sjdbOverhang
    out: [star_ref, debug_log]
