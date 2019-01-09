cwlVersion: v1.0
class: CommandLineTool
id: cutadapter
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'migbro/cutadapt:latest'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 16
    ramMin: 24000

baseCommand: [cutadapt]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -j 16
      -m 20
      --quality-base=33
      -q 20
      -a $(inputs.r1_adapter)
      -A $(inputs.r2_adapter)
      -o TRIMMED.$(inputs.readFilesIn1.basename)
      -p TRIMMED.$(inputs.readFilesIn2.basename)
      $(inputs.readFilesIn1.path)
      $(inputs.readFilesIn2.path)
      > $(inputs.sample_name).cutadapt_results.txt

inputs:
  r1_adapter: string
  r2_adapter: string
  readFilesIn1: File
  readFilesIn2: File
  sample_name: string

outputs:
  trimmedReadsR1:
    type: File
    outputBinding:
      glob: TRIMMED.$(inputs.readFilesIn1.basename)
  trimmedReadsR2:
    type: File
    outputBinding:
      glob: TRIMMED.$(inputs.readFilesIn2.basename)
  cutadapt_stats:
    type: File
    outputBinding:
      glob: $(inputs.sample_name).cutadapt_results.txt
