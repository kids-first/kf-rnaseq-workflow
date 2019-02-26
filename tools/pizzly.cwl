cwlVersion: v1.0
class: CommandLineTool
id: pizzly
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'migbro/pizzly:latest'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 4
    ramMin: 24000

baseCommand: [pizzly]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -k 31 --gtf $(inputs.GTF.path)
      --align-score 2
      --insert-size 400
      --fasta $(inputs.transcript_fa.path)
      --output $(inputs.SampleID)
      $(inputs.kallisto_fusion.path) &&
      python3 /pizzly/scripts/get_fragment_length.py $(inputs.SampleID).json > $(inputs.SampleID).pizzly.flattened.txt

inputs:
  transcript_fa: File
  GTF: File
  kallisto_fusion: File
  SampleID: string

outputs:
  fusions_flattened:
    type: File
    outputBinding:
      glob: "$(inputs.SampleID).pizzly.flattened.txt"

