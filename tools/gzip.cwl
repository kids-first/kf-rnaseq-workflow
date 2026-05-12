cwlVersion: v1.2
class: CommandLineTool
id: gzip
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.cpu)
  - class: DockerRequirement
    dockerPull: ubuntu:24.04
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      gzip -c
  - position: 10
    shellQuote: false
    valueFrom: >-
      $(inputs.infile.path) > $(inputs.decompress ? inputs.infile.basename.replace(/.gz$/, "") : inputs.infile.basename + ".gz") 
inputs:
  infile: { type: 'File', doc: "reference sequence file" }
  decompress: { type: 'boolean?', inputBinding: { position: 2, prefix: "--decompress" }, doc: "decompress" }
  name: { type: 'boolean?', inputBinding: { position: 2, prefix: "--name" }, doc: "save or restore the original name and timestamp" }
  no_name: { type: 'boolean?', inputBinding: { position: 2, prefix: "--no-name" }, doc: "do not save or restore the original name and timestamp" }
  fast: { type: 'boolean?', inputBinding: { position: 2, prefix: "--fast" }, doc: "compress faster" }
  best: { type: 'boolean?', inputBinding: { position: 2, prefix: "--best" }, doc: "compress better" }
  ram: { type: 'int?', default: 2, doc: "GB memory to allocate to this task" }
  cpu: { type: 'int?', default: 1, doc: "CPUs to allocate to this task." }
outputs:
  outfile:
    type: File
    outputBinding:
      glob: |
        $(inputs.decompress ? inputs.infile.basename.replace(/.gz$/, "") : inputs.infile.basename + ".gz")
