cwlVersion: v1.0
class: CommandLineTool
id: annoFuse
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/annofuse:0.90.0'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 4
    ramMin: 8000

baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      A_CT=`wc -l $(inputs.arriba_formatted_fusions.path) | cut -f 1 -d " "`

      S_CT=`wc -l $(inputs.starfusion_formatted_fusions.path) | cut -f 1 -d " "`

      if [ $A_CT -eq 1 ] && [ $S_CT -eq 1 ]; then
        echo "Both inputs are empty, will skip processing as there no fusions." >&2;
        exit 0;
      fi

      Rscript /rocker-build/annoFusePerSample.R
      --fusionfileArriba $(inputs.arriba_formatted_fusions.path)
      --fusionfileStarFusion $(inputs.starfusion_formatted_fusions.path)
      --expressionFile $(inputs.rsem_expr_file.path)
      --tumorID $(inputs.sample_name)
      --outputfile $(inputs.output_basename).annoFuse_filter.tsv

inputs:
  arriba_formatted_fusions: {type: File, doc: "arriba fusion file formatted by format_fusion_file.cwl, and annotated by fusion_annotator.cwl"}
  starfusion_formatted_fusions: {type: File, doc: "STARFusion file formatted by format_fusion_file.cwl"}
  rsem_expr_file: {type: File, doc: "gzipped rsem gene expression file"}
  sample_name: string
  output_basename: string

outputs:
  filtered_fusions_tsv:
    type: ['null', File]
    outputBinding:
      glob: '*.tsv'
