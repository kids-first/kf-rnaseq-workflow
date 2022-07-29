cwlVersion: v1.0
class: ExpressionTool
id: basename-picker
requirements:
  - class: InlineJavascriptRequirement

inputs:
  read1_filename: string
  output_basename: 'string?'

outputs:
  output:
    type: string

expression:
  "${
    return {'output' : inputs.output_basename ? inputs.output_basename : inputs.read1_filename};
    }"
