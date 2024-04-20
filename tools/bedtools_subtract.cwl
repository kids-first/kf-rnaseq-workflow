cwlVersion: v1.2
class: CommandLineTool
id: bedtools_subtract
doc: "Remove entries from input vcf file based on bed file with coordinates to subtract"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 4
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/vcfutils:latest'

baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      'bedtools subtract -wa
  - position: 2
    shellQuote: false
    valueFrom: >-
      | bgzip -@ 4 -c > $(inputs.output_basename).bed_subtract.vcf.gz && tabix $(inputs.output_basename).bed_subtract.vcf.gz'

inputs:
    input_vcf: { type: File, secondaryFiles: ['.tbi'], doc: "Input VCF file.",
      inputBinding: { position: 1, prefix: "-a"} }
    subtract_bed: { type: File, doc: "Bed file to subtract entries from VCF",
      inputBinding: { position: 1, prefix : "-b"} }
    output_basename: string
    output_header_flag: { type: 'boolean?', doc: "Output header from input file", default: true,
      inputBinding: { position: 1, prefix: "-header"} }

outputs:
  subtracted_vcf:
    type: File
    outputBinding:
      glob: '*.bed_subtract.vcf.gz'
    secondaryFiles: ['.tbi']