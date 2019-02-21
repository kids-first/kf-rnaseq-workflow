cwlVersion: v1.0
class: CommandLineTool
id: fusion_catcher
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'migbro/fusioncatcher:latest'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 36
    ramMin: 60000

baseCommand: [tar, -xzf]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.ensembl_genome.path) &&
      mkdir INPUTS &&
      ln -s $(inputs.readFilesIn1.path) INPUTS/ &&
      ln -s $(inputs.readFilesIn2.path) INPUTS/ &&
      /opt/fusioncatcher/v1.00/bin/fusioncatcher
      -p 36
      -o OUTPUT
      -i INPUTS
      -d ./$(inputs.ensembl_genome.nameroot)/
      --skip-blat

inputs:
  ensembl_genome: File
  readFilesIn1: File
  readFilesIn2: File

outputs:
  fusion_caught:
    type: File
    outputBinding:
      glob: "./OUTPUT/*"