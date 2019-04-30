cwlVersion: v1.0
class: CommandLineTool
id: star_fusion
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/star-fusion:1.5.0'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 16
    ramMin: 64000

baseCommand: [tar]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -zxf $(inputs.genomeDir.path) &&
      /usr/local/STAR-Fusion/STAR-Fusion
      --genome_lib_dir ./GRCh38_v27_CTAT_lib_Feb092018/ctat_genome_lib_build_dir
      -J $(inputs.Chimeric_junction.path)
      --output_dir STAR-Fusion_outdir
      --examine_coding_effect
      --CPU 16 &&
      mv STAR-Fusion_outdir/star-fusion.fusion_predictions.abridged.coding_effect.tsv $(inputs.SampleID).STAR.fusion_predictions.abridged.coding_effect.tsv &&
      gzip -c $(inputs.Chimeric_junction.path) > $(inputs.Chimeric_junction.basename).gz
      

inputs:
  Chimeric_junction: File
  genomeDir: File
  SampleID: string

outputs:
  abridged_coding:
    type: File
    outputBinding:
      glob: '*.fusion_predictions.abridged.coding_effect.tsv'
  chimeric_junction_compressed:
    type: File
    outputBinding:
      glob: "$(inputs.Chimeric_junction.basename).gz"
