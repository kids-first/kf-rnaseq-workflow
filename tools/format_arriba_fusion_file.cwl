cwlVersion: v1.0
class: CommandLineTool
id: format_arriba_fusion_file
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 4000
    coresMin: 2
    https://platform.illumina.com/rdf/ica/resources:tier: economy
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/annofuse:0.92.0'
  - class: InitialWorkDirRequirement
    listing:
      - entryname: formatArribaFusionCalls.R
        entry:
          $include: ../scripts/formatArribaFusionCalls.R

baseCommand: [Rscript, formatArribaFusionCalls.R]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      --outputfile $(inputs.sample_name).arriba_formatted.tsv

inputs:
  arriba_fusion_file: { type: File, doc: "Arriba fusion file", inputBinding: { prefix: "--fusionfile", position: 1 } }
  sample_name: { type: string, doc: "Sample name", inputBinding: { prefix: "--tumorid", position: 1 } }

outputs:
  formatted_fusion_tsv:
    type: File
    outputBinding:
      glob: '*.tsv'
