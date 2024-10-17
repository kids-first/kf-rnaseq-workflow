cwlVersion: v1.2
class: ExpressionTool
id: basename_picker
requirements:
  - class: InlineJavascriptRequirement

inputs:
  root_name: string
  output_basename: 'string?'
  sample_name: 'string?'
  star_rg_line: 'string?'

outputs:
  outname:
    type: string
  outsample:
    type: string
  outrg:
    type: string

expression: |
  ${
    var name = inputs.output_basename ? inputs.output_basename : inputs.root_name;
    var sample = inputs.sample_name ? inputs.sample_name : inputs.root_name;
    var rgid = "ID:" + sample + "_1"
    var rglb = "LB:" + sample
    var rgsm = "SM:" + sample
    var rgpl = "PL:Illumina"
    var rgds = "DS:\"Values for this read group were auto-generated and may not reflect the true read group information.\""
    var rg = inputs.star_rg_line ? inputs.star_rg_line : [rgid, rgpl, rglb, rgsm, rgds].join("\t")
    return {
      'outname' : name,
      'outsample' : sample,
      'outrg' : rg
    }
  }
