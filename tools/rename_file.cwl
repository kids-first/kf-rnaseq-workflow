cwlVersion: v1
class: CommandLineTool
baseCommand: [bash, -c]
inputs:
  input_file:
    type: File
    inputBinding:
      position: 1
  new_name:
    type: string
    inputBinding:
      position: 2
outputs:
  renamed_file:
    type: File
    outputBinding:
      glob: $(inputs.new_name)
arguments:
  - valueFrom: |
      cp $(inputs.input_file.path) $(inputs.new_name)

