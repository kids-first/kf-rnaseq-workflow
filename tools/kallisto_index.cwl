cwlVersion: v1.2
class: CommandLineTool
id: kallisto_index
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cpu)
    ramMin: $(inputs.ram * 1000)
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/uros_sipetic/kallisto:0.43.1'

baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      kallisto index

inputs:
  transcripts_fasta: { type: 'File', inputBinding: { position: 9 }, doc: "Transcripts Fasta" }
  output_filename: { type: 'string?', default: "kallisto.idx", inputBinding: { position: 2, prefix: "-i" }, doc: "Name for output" }
  cpu: { type: 'int?', default: 8, doc: "CPUs to use for this task." }
  ram: { type: 'int?', default: 16, doc: "GB of RAM to use for this task." }

outputs:
  index:
    type: File
    outputBinding:
      glob: '$(inputs.output_filename)'
