cwlVersion: v1.0
class: CommandLineTool
id: format_fusion_file
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/annofuse:0.90.0'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 4
    ramMin: 8000

baseCommand: [Rscript]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      /rocker-build/formatFusionCalls.R
      --fusionfile $(inputs.input_caller_fusion_file.path)
      --tumorid $(inputs.sample_name)
      --caller $(inputs.caller)
      --outputfile $(inputs.sample_name).$(inputs.caller)_formatted.tsv

inputs:
  input_caller_fusion_file: File
  sample_name: string
  caller: {type: [{type: enum, name: caller, symbols: ["arriba", "starfusion"]}], doc: "Source of calls, currently support algorithms: arriba and STARFusion"}

outputs:
  formatted_fusion_tsv:
    type: File
    outputBinding:
      glob: '*.tsv'
