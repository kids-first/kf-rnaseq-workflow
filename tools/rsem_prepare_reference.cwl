cwlVersion: v1.2
class: CommandLineTool
id: rsem_prepare_reference
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cpu)
    ramMin: $(inputs.ram * 1000)
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/uros_sipetic/rsem:1.3.1'

baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      mkdir $(inputs.output_prefix)
      && gunzip -c $(inputs.gtf.path) > annotations.gtf
      && gunzip -c $(inputs.reference.path) > genome.fa
      && rsem-prepare-reference --gtf annotations.gtf
  - position: 9
    shellQuote: false
    valueFrom: >-
      genome.fa $(inputs.output_prefix)/$(inputs.output_prefix)
  - position: 10
    shellQuote: false
    valueFrom: >-
      && tar -czf $(inputs.output_prefix).tar.gz $(inputs.output_prefix)

inputs:
  gtf: { type: 'File', doc: "GTF file. MUST BE GZIPPED." }
  reference: { type: 'File', doc: "Reference fasta. MUST BE GZIPPED." }
  output_prefix: { type: 'string?', default: "rsem1.3.1", doc: "Name for output" }
  cpu: { type: 'int?', default: 8, inputBinding: { position: 2, prefix: "--num-threads" }, doc: "CPUs to allocate to this task." }
  ram: { type: 'int?', default: 16, doc: "GB of RAM to allocate to this task." }

outputs:
  genome_tar:
    type: File
    outputBinding: 
      glob: '*tar.gz'
  transcripts_fasta:
    type: File
    outputBinding:
      glob: '$(inputs.output_prefix)/*transcripts.fa'
