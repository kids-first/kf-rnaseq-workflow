cwlVersion: v1.2
class: CommandLineTool
id: samtools_bam_to_cram
doc: |-
  This tool converts the input BAM into a CRAM.
  The following programs are run in this tool:
    - samtools view
    - samtools index
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cores)
    ramMin: $(inputs.ram * 1000)
    https://platform.illumina.com/rdf/ica/resources:tier: economy
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9'
baseCommand: [samtools, view]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -C -o $(inputs.input_bam.nameroot).cram
  - position: 4
    shellQuote: false
    valueFrom: >-
      && samtools index $(inputs.input_bam.nameroot).cram
inputs:
  reference: {type: File, secondaryFiles: [.fai], doc: "Reference fasta with associated fai index",
    inputBinding: {prefix: "-T", position: 2 }}
  input_bam: {type: File, doc: "Input bam file",
    inputBinding: { position: 3 }}
  cores: {type: 'int?', default: 16, inputBinding: {prefix: "-@", position: 2}}
  ram: {type: 'int?', default: 32, doc: "GB of RAM to allocate to this task"}
outputs:
  output: { type: File, outputBinding: { glob: '*.cram' }, secondaryFiles: [.crai] }
