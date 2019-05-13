cwlVersion: v1.0
class: CommandLineTool
id: build_star_align_ref
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/star:latest'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: "$(inputs.runThreadN ? inputs.runThreadN : 16)"
    ramMin: 60000

baseCommand: [mkdir]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.genomeDir) &&
      STAR --runMode genomeGenerate --runThreadN $(inputs.runThreadN) --genomeDir $(inputs.genomeDir) --genomeFastaFiles $(inputs.genome_fa.path) --sjdbGTFfile $(inputs.gtf.path) --sjdbOverhang $(inputs.sjdbOverhang)
      && tar -czf $(inputs.genomeDir).tar.gz ./$(inputs.genomeDir)

inputs:
  genomeDir: string
  genome_fa: File
  gtf: File
  runThreadN: int
  sjdbOverhang: int

outputs:
  star_ref:
    type: File
    outputBinding:
      glob: $(inputs.genomeDir).tar.gz