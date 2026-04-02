cwlVersion: v1.2
class: CommandLineTool
id: gtex_collapse_annotation
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.cpu)
  - class: DockerRequirement
    dockerPull: pgc-images.sbgenomics.com/danmiller/gtex_collapse:0.1.0
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      collapse_annotation.py

inputs:
  gtf: { type: 'File', inputBinding: { position: 8 }, doc: "Input GTF file. Can be GZIP." }
  output_filename: { type: 'string?', default: "output.gtf", inputBinding: { position: 9 }, doc: "Name for output GTF file. Program will GZIP if this ends in .gz" }
  transcript_blacklist: { type: 'File?', inputBinding: { position: 2, prefix: '--transcript_blacklist' }, doc: "List of transcripts to exclude (e.g., unannotated readthroughs)" }
  collapse_only: { type: 'boolean?', inputBinding: { position: 2, prefix: '--collapse_only' }, doc: "Only collapse transcripts of each gene, do not remove overlaps." }
  stranded: { type: 'boolean?', inputBinding: { position: 2, prefix: '--stranded' }, doc: "Only consider genes on the same strand when removing overlaps." }
  cpu: { type: 'int?', doc: "Num processing threads to use", default: 2 }
  ram: { type: 'int?', doc: "Num GB memory to make available", default: 4 }
outputs:
  collapsed_gtf: { type: 'File', outputBinding: { glob: '$(inputs.output_filename)' } }
