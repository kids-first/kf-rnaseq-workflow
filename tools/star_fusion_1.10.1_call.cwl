cwlVersion: v1.2
class: CommandLineTool
id: star_fusion_1-10-1_call
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'migbro/star_fusion:1.10.1'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cores)
    ramMin: 64000

baseCommand: [tar, -I pigz, -xf]
arguments:
  - position: 2
    shellQuote: false
    valueFrom: >-
      && /usr/local/STAR-Fusion/STAR-Fusion
      --genome_lib_dir ./$(inputs.genome_untar_path)
      --output_dir STAR-Fusion_outdir
  - position: 3
    shellQuote: false
    valueFrom: >-
      && mv STAR-Fusion_outdir/star-fusion.fusion_predictions.abridged.coding_effect.tsv $(inputs.output_basename).STAR-1.10.1.fusion_predictions.abridged.coding_effect.tsv &&
      pigz -c $(inputs.Chimeric_junction.path) > $(inputs.Chimeric_junction.basename).gz
      

inputs:
  genome_tar: { type: File, doc: "CTAT library, either downloaded or generated through an arduous process", inputBinding: { position: 1 } }
  Chimeric_junction: { type: File, doc: "Output junction file from STAR", inputBinding: { prefix: '-J', position: 2 } }
  genome_untar_path: {type: ['null', string], doc: "This is what the path will be when genome_tar is unpackaged", default: "GRCh38_v38_CTAT_lib_Feb262022.CUSTOM"}
  examine_coding_effect: { type: 'string?', doc: "Flag to activate examining coding effect. Use null to skip",
  default: "--examine_coding_effect", inputBinding: { position: 2 } }
  cores: { type: 'int?', doc: "Num cpus to use, >16 has diminishing returns", default: 16, inputBinding: { prefix: '--CPU', position: 2 } }
  output_basename: string

outputs:
  abridged_coding:
    type: File
    outputBinding:
      glob: '*.fusion_predictions.abridged.coding_effect.tsv'
  chimeric_junction_compressed:
    type: File
    outputBinding:
      glob: "$(inputs.Chimeric_junction.basename).gz"
