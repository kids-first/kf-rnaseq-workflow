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

baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      mkdir $(inputs.genomeDir)
      && gunzip -c $(inputs.genome_fa.path) > genome.fa
      && gunzip -c $(inputs.gtf.path) > annotations.gtf
  - position: 2
    shellQuote: false
    valueFrom: >-
      && STAR --runMode genomeGenerate --genomeFastaFiles genome.fa --sjdbGTFfile annotations.gtf
  - position: 4
    shellQuote: false
    valueFrom: >-
      && tar -czf $(inputs.genomeDir).tar.gz ./$(inputs.genomeDir)

inputs:
  genomeDir: { type: string, doc: "Output dirname. Recommend STAR_{version}_GENCODE{version num}", inputBinding: { position: 3, prefix: '--genomeDir' } }
  genome_fa: { type: File, doc: "Fasta file to index. Recommend from GENCODE, PRI assembly." }
  gtf: { type: File, doc: "Matched GTF file to index. Recommend from GENCODE, PRI assembly" }
  runThreadN: { type: 'int?', default: 16, inputBinding: { position: 3, prefix: '--runThreadN' } }
  sjdbOverhang: { type: 'int?', default: 100, doc: "Ideal value is read len minus 1, but default 100 ok for most cases", inputBinding: { position: 3, prefix: '--sjdbOverhang' } }

outputs:
  star_ref:
    type: File
    outputBinding:
      glob: $(inputs.genomeDir).tar.gz
