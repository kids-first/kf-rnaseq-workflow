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
  tool_name: string
  include_expression: 'string?'
  exclude_expression: 'string?'
  # Genome gen vars
  genomeDir: { type: string, doc: "Output dirname. Recommend STAR_{version}_GENCODE{version num}_{Patient/sample id}" }
  genome_fa: { type: File, doc: "Fasta file to index. Recommend from GENCODE, PRI assembly. Must unzip first if compressed" }
  genomeTransformType: { type: [ 'null', {type: enum, name: genomeTransformType, symbols: [
      "None",
      "Haploid",
      "Diploid"
      ]}],
  doc: "type of genome transformation - None: no transformation. Haploid: eplace reference alleles with alternative alleles from VCF file (e.g. consensus allele) \
  Diploid: create two haplotypes for each chromosome listed in VCF file, for genotypes 1—2, assumes perfect phasing (e.g. personal genome)" }
  gtf: { type: File, doc: "Matched GTF file to index. Recommend from GENCODE, PRI assembly" }
  runThreadN: { type: 'int?', default: 16 }
  memory: { type: 'int?', doc: "Mem in GB required. With no VCF, 60DB is fine, need more with VCF", default: 60}
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
      tool_name: tool_name
      strip_info: strip_info
    out: [stripped_vcf]
  bcftools_subset_vcf:
    run: ../tools/bcftools_filter_vcf.cwl
    when: $(inputs.include_expression != null || inputs.exclude_expression != null)
    in:
      input_vcf: 
        source: [bcftools_strip_info/stripped_vcf, input_vcf]
        pickValue: first_non_null
      include_expression: include_expression
      exclude_expression: exclude_expression
      output_basename: output_basename
    out: [filtered_vcf]
  gatk4_selectvariants:
    run: ../tools/gatk_selectvariants.cwl
    in:
      input_vcf:
        source: [bcftools_subset_vcf/filtered_vcf, bcftools_strip_info/stripped_vcf, input_vcf]
        pickValue: first_non_null
      output_basename: output_basename
      tool_name: tool_name
    out: [pass_vcf]
  star_personal_genome_generate:
    run: ../tools/star_2.7.11b_personal_genome_generate.cwl
    in:
      genomeDir: genomeDir
      genome_fa: genome_fa
      genomeTransformVCF: gatk4_selectvariants/pass_vcf
      genomeTransformType: genomeTransformType
      gtf: gtf
      runThreadN: runThreadN
      memory: memory
      sjdbOverhang: sjdbOverhang
    out: [star_ref, debug_log]
