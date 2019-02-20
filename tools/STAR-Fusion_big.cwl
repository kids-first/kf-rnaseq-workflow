cwlVersion: v1.0
class: CommandLineTool
id: star_fusion
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'trinityctat/ctatfusion:1.5.0'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 16
    ramMin: 100000

baseCommand: [tar]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -zxf $(inputs.genomeDir.path) &&
      /usr/local/src/STAR-Fusion/STAR-Fusion
      --genome_lib_dir ./GRCh38_v27_CTAT_lib_Feb092018/ctat_genome_lib_build_dir
      -J $(inputs.Chimeric.path)
      --output_dir STAR-Fusion_outdir
      --CPU 16 &&
      mv STAR-Fusion_outdir/star-fusion.fusion_predictions.tsv $(inputs.SampleID).fusion_predictions.tsv


inputs:
  Chimeric: File
  genomeDir: File
  SampleID: string

outputs:
  fusion_out:
    type: File
    outputBinding:
      glob: '*.fusion_predictions.tsv'


