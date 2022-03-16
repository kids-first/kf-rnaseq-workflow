cwlVersion: v1.0
class: CommandLineTool
id: rsem-prepare-reference
label: "RSEM v1.3.1 Prepare Reference"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/uros_sipetic/rsem:1.3.1'
  - class: InlineJavascriptRequirement

baseCommand: [mkdir]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.reference_name)
      && rsem-prepare-reference
  - position: 4
    shellQuote: false
    valueFrom: >-
      && mv $(inputs.reference_name).* $(inputs.reference_name)/
      && tar -czf $(inputs.reference_name).tar.gz $(inputs.reference_name)

inputs:
  reference_fasta: { type: File, doc: "Reference fasta file", inputBinding: { position: 3} }
  reference_name: { type: string, doc: "Output file prefix. Recommend format: RSEM_<SOURCE><Version>/", inputBinding: { position: 3}}
  reference_gtf: { type: 'File?', doc: "gene model definitions. This OR gff required" , inputBinding: { position: 2, prefix: '--gtf' } }
  reference_gff: { type: 'File?', doc: "gene model definitions. This OR gtf required" , inputBinding: { position: 2, prefix: '--gff' } }

outputs:
  rsem_reference:
    type: File
    outputBinding: 
      glob: '*tar.gz'

