cwlVersion: v1.2
class: CommandLineTool
id: star_2-7-10a_build_ref
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/star:2.7.10a'
  - class: ResourceRequirement
    coresMin: $(inputs.cpu)
    ramMin: $(inputs.ram * 1000)

baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      mkdir $(inputs.genomeDir)
  - position: 10
    shellQuote: false
    valueFrom: >-
      && STAR --runMode genomeGenerate
  - position: 20
    shellQuote: false
    valueFrom: >-
      && tar -czf $(inputs.genomeDir).tar.gz ./$(inputs.genomeDir)

inputs:
  genomeDir: { type: string, inputBinding: { position: 12, prefix: "--genomeDir" }, doc: "Output dirname. Recommend STAR_{version}_GENCODE{version num}" }
  genome_fa: { type: File, inputBinding: { position: 12, prefix: "--genomeFastaFiles" }, doc: "Fasta file to index. Recommend from GENCODE, PRI assembly. MUST BE UNZIPPED" }
  gtf: { type: File, inputBinding: { position: 12, prefix: "--sjdbGTFfile" }, doc: "GTF file (matched to genome_fa) to index. Recommend from GENCODE, PRI assembly. MUST BE UNZIPPED" }
  runThreadN: { type: 'int?', default: 16, inputBinding: { position: 12, prefix: '--runThreadN' } }
  sjdbOverhang: { type: 'int?', default: 100, inputBinding: { position: 12, prefix: '--sjdbOverhang' }, doc: "Ideal value is read len minus 1, but default 100 ok for most cases" }
  cpu: { type: 'int?', default: 16, inputBinding: { position: 12, prefix: '--runThreadN' }, doc: "CPUs to allocate to this task" }
  ram: { type: 'int?', default: 60, doc: "GB of RAM to allocate to this task" }

outputs:
  star_ref:
    type: File
    outputBinding:
      glob: $(inputs.genomeDir).tar.gz
