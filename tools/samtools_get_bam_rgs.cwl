cwlVersion: v1.0
class: CommandLineTool
id: samtools_get_bam_rgs
doc: |-
  This tool takes the input bam and returns a file containing the @RG lines from the header.
requirements:
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
  - class: ShellCommandRequirement
baseCommand: ["/bin/bash","-c"]
arguments:
  - position: 0
    shellQuote: true
    valueFrom: |-
      set -eo pipefail

      samtools view -H $(inputs.input_bam.path) | grep ^@RG > rg.txt
inputs:
  input_bam: { type: File, doc: "Input BAM file" }
outputs:
  output:
    type: File
    outputBinding:
      glob: "rg.txt"
