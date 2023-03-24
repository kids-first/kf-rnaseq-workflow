cwlVersion: v1.2
class: CommandLineTool
id: bam_readlength
doc: |-
  samtools view file.bam | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cores)
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9'
baseCommand: [samtools,view]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.input_bam.path) | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c > $(inputs.input_bam.nameroot).bam_readlength
  
inputs:
  input_bam: { type: File, doc: "Input bam file"}

outputs:
  output:
    type: File
    outputBinding:
      glob: "*.bam_readlength"