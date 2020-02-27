cwlVersion: v1.0
class: CommandLineTool
id: annoFuse
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'gaonkark/annofuse:latest'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 4
    ramMin: 8000

baseCommand: [Rscript]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      /rocker-build/annoFusePerSample.R
      --fusionfileArriba $(inputs.arriba_formatted_fusions.path)
      --fusionfileStarFusion $(inputs.starfusion_formatted_fusions.path)
      --outputfile $(inputs.output_basename).annoFuse_filter.tsv

inputs:
  arriba_formatted_fusions: {type: File, doc: "arriba fusion file formatted by format_fusion_file.cwl, and annotated by fusion_annotator.cwl"}
  starfusion_formatted_fusions: {type: File, doc: "STARFusion file formatted by format_fusion_file.cwl"}
  output_basename: string

outputs:
  filtered_fusions_tsv:
    type: File
    outputBinding:
      glob: '*.tsv'
