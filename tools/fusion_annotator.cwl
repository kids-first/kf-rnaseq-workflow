cwlVersion: v1.0
class: CommandLineTool
id: fusion_annotator
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/fusionannotator:0.1.1'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 4
    ramMin: 8000
    https://platform.illumina.com/rdf/ica/resources:tier: economy

baseCommand: [tar]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -I pigz -xf $(inputs.genome_tar.path) &&
      /opt/FusionAnnotator/FusionAnnotator
      --genome_lib_dir ./$(inputs.genome_untar_path)
      --annotate $(inputs.input_fusion_file.path)
      --fusion_name_col $(inputs.col_num)
      > $(inputs.output_basename).annotated.tsv

inputs:
  input_fusion_file: {type: File, doc: "Fusion file formatted by format_fusion_file.cwl, likely from arriba input"}
  genome_tar: File
  genome_untar_path: {type: ['null', string], doc: "This is what the path will be when genome_tar is unpackaged", default: "GRCh38_v39_CTAT_lib_Mar242022.CUSTOM"}
  col_num: {type: ['null', int], doc: "column number in file of fusion name, use 24 for arriba v1.1, 30 for v2", default: 30}
  output_basename: string

outputs:
  annotated_tsv:
    type: File
    outputBinding:
      glob: '*.tsv'
