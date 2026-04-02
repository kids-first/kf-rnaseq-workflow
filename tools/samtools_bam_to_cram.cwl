cwlVersion: v1.2
class: CommandLineTool
id: samtools_bam_to_cram
doc: |-
  This tool converts the input BAM into a CRAM.
  The following programs are run in this tool:
    - samtools view
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cores)
    ramMin: 16000
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/samtools:1.20-multi-arch'
baseCommand: [samtools, view]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -C -o $(inputs.input_bam.nameroot).cram##idx##$(inputs.input_bam.nameroot).cram.crai
inputs:
  reference: {type: File, secondaryFiles: [.fai], doc: "Reference fasta with associated fai index",
    inputBinding: {prefix: "-T", position: 2 }}
  input_bam: {type: File, doc: "Input bam file",
    inputBinding: { position: 3 }}
  cores: {type: 'int?', default: 8, inputBinding: {prefix: "-@", position: 2}}
outputs:
  output: { type: File, outputBinding: { glob: '*.cram' }, secondaryFiles: [.crai] }
