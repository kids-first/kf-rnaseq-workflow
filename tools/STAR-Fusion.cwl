cwlVersion: v1.0
class: CommandLineTool
id: star_fusion
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'trinityctat/ctatfusion:1.5.0'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 64000

baseCommand: [tar]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -zxf $(inputs.genomeDir.path) &&
      /usr/local/src/STAR-Fusion/STAR-Fusion
      --left_fq $(inputs.readFilesIn1.path)
      --right_fq $(inputs.readFilesIn2.path)
      --genome_lib_dir ./GRCh38_v27_CTAT_lib_Feb092018/ctat_genome_lib_build_dir
      -J $(inputs.Chimeric_junction.path)
      --output_dir STAR-Fusion_outdir
      --examine_coding_effect --denovo_reconstruct --FusionInspector inspect
      --CPU 8 &&
      mv STAR-Fusion_outdir/star-fusion.fusion_predictions.abridged.coding_effect.tsv $(inputs.SampleID).STAR.fusion_predictions.abridged.coding_effect.tsv
      

inputs:
  readFilesIn1: File
  readFilesIn2: File
  Chimeric_junction: File
  genomeDir: File
  SampleID: string

outputs:
  abridged_coding:
    type: File
    outputBinding:
      glob: '*.fusion_predictions.abridged.coding_effect.tsv'
