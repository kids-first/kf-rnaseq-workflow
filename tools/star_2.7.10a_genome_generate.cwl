cwlVersion: v1.2
class: CommandLineTool
id: star_2-7-10a_build_ref
label: "STAR Reference Genome Index Generator"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/star:2.7.10a'
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
  genomeDir: { type: string, doc: "Output dirname. Recommend STAR_{version}_GENCODE{version num}", inputBinding: { position: 3, prefix: '--genomeDir' } }
  genome_fa: { type: File, doc: "Fasta file to index. Recommend from GENCODE, PRI assembly. Must unzip first if compressed", inputBinding: { position: 3, prefix: '--genomeFastaFiles' } }
  gtf: { type: File, doc: "Matched GTF file to index. Recommend from GENCODE, PRI assembly", inputBinding: { position: 3, prefix: '--sjdbGTFfile' } }
  runThreadN: { type: 'int?', default: 16, inputBinding: { position: 3, prefix: '--runThreadN' } }
  sjdbOverhang: { type: 'int?', default: 100, doc: "Ideal value is read len minus 1, but default 100 ok for most cases", inputBinding: { position: 3, prefix: '--sjdbOverhang' } }

outputs:
  star_ref:
    type: File
    outputBinding:
      glob: $(inputs.genomeDir).tar.gz
