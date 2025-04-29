cwlVersion: v1.2
class: CommandLineTool
id: basename_picker
requirements:
  - class: InlineJavascriptRequirement

baseCommand: [echo, done]

inputs:
  root_name: { type: 'string' }
  output_basename: { type: 'string?'} 
  sample_name: { type: 'string?'}
  star_rg_line: { type: 'string?'}

outputs:
  outname:
    type: string
    outputBinding:
      outputEval: |
        $(inputs.output_basename ? inputs.output_basename : inputs.root_name)
  outsample:
    type: string
    outputBinding:
      outputEval: |
        $(inputs.sample_name ? inputs.sample_name : inputs.root_name)
  outrg:
    type: string
    outputBinding:
      outputEval: |
        ${
          var sample = inputs.sample_name ? inputs.sample_name : inputs.root_name;
          var rgid = "ID:" + sample + "_1"
          var rgid = "ID:" + sample + "_1"
          var rglb = "LB:" + sample
          var rgsm = "SM:" + sample
          var rgpl = "PL:Illumina"
          var rgds = "DS:\"Values for this read group were auto-generated and may not reflect the true read group information.\""
          var rg = inputs.star_rg_line ? inputs.star_rg_line : [rgid, rgpl, rglb, rgsm, rgds].join("\t")
          return rg
        }
