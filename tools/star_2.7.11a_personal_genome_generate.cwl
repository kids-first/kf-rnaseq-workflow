cwlVersion: v1.2
class: CommandLineTool
id: star_2-7-11a_build_personal_ref
label: "STAR Personal Reference Genome Index Generator"
doc: " Allow user to create a custom genome to align aganist"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/star:2.7.11a'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.runThreadN)
    ramMin: 60000

baseCommand: [mkdir]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.genomeDir)
  - position: 2
    shellQuote: false
    valueFrom: >-
      && STAR --runMode genomeGenerate
  - position: 4
    shellQuote: false
    valueFrom: >-
      && tar -czf $(inputs.genomeDir).tar.gz ./$(inputs.genomeDir)

inputs:
  genomeDir: { type: string, doc: "Output dirname. Recommend STAR_{version}_GENCODE{version num}", inputBinding: { position: 3, prefix: "--genomeDir" } }
  genome_fa: { type: File, doc: "Fasta file to index. Recommend from GENCODE, PRI assembly. Must unzip first if compressed", inputBinding: { position: 3, prefix: "--genomeFastaFiles" } }
  genomeTransformVCF: { type: 'File?', doc: "path to UNCOMPRESSED VCF file for genome transformation", inputBinding: { position: 3, prefix: "--genomeTransformVCF"} }
  genomeTransformType: { type: [ 'null', {type: enum, name: genomeTransformType, symbols: [
      "None",
      "Haploid",
      "Diploid"
      ]}],
  doc: "type of genome transformation - None: no transformation. Haploid: eplace reference alleles with alternative alleles from VCF file (e.g. consensus allele) \
  Diploid: create two haplotypes for each chromosome listed in VCF file, for genotypes 1—2, assumes perfect phasing (e.g. personal genome)",
  inputBinding: { position: 3, prefix: "--genomeTransformType", shellQuote: false } }
  gtf: { type: File, doc: "Matched GTF file to index. Recommend from GENCODE, PRI assembly", inputBinding: { position: 3, prefix: "--sjdbGTFfile" } }
  runThreadN: { type: 'int?', default: 16, inputBinding: { position: 3, prefix: "--runThreadN" } }
  sjdbOverhang: { type: 'int?', default: 100, doc: "Ideal value is read len minus 1, but default 100 ok for most cases", inputBinding: { position: 3, prefix: "--sjdbOverhang" } }

outputs:
  star_ref:
    type: File
    outputBinding:
      glob: $(inputs.genomeDir).tar.gz
  debug_log:
    type: File
    outputBinding:
      glob: $(inputs.genomeDir)/Log.out

