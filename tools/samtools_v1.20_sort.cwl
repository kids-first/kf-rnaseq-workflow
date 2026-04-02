cwlVersion: v1.2
class: CommandLineTool
id: samtools_v1-20_sort
label: "Samtools v1.20 coordinate sort"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/samtools:1.20-multi-arch'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cores)
    ramMin: ${ return inputs.cores * 1500 }

baseCommand: [samtools]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      sort $(inputs.unsorted_bam.path)
      -@ $(inputs.cores)
      -m 1500M
      -O bam
      --write-index
      -o $(inputs.unsorted_bam.nameroot).sorted.bam##idx##$(inputs.unsorted_bam.nameroot).sorted.bai
inputs:
  cores: { type: 'int?', doc: "Num cores to use for sorting", default: 8}
  unsorted_bam: { type: File, doc: "Bam to sort, likely from STAR" }

outputs:
  sorted_bam:
    type: File
    outputBinding:
      glob: '*.sorted.bam'
  sorted_bai:
    type: File
    outputBinding:
      glob: '*.bai'
