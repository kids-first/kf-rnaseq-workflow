cwlVersion: v1.2
class: CommandLineTool
id: bcftools_filter_vcf
doc: "More generic tool to take in an include expression and optionally an exclude expresssion to filter a vcf"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/bcftools:1.20'
  - class: ResourceRequirement
    ramMin: 16000
    coresMin: 8
  - class: InlineJavascriptRequirement
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      ${
        var cmd  = "bcftools view"
        if ( inputs.sample_name != null ){
          cmd += " --threads " + inputs.threads + " -s " + inputs.sample_name + " " + inputs.input_vcf.path + " | bcftools view ";
        }
        return cmd;
      }
  - position: 2
    shellQuote: false
    valueFrom: >-
      ${
        var arg = " -o " + inputs.output_basename + ".bcf_filtered"
        if (inputs.output_type == "v"){
          arg += ".vcf"
        } else if (inputs.output_type == "z"){
          arg += ".vcf.gz && tabix " + inputs.output_basename + ".bcf_filtered.vcf.gz"
        } else if (inputs.output_type == "b"){
          arg += ".bcf.gz"
        } else{
          arg += ".bcf"
        }
        if (inputs.sample_name == null){
          arg = inputs.input_vcf.path + arg;
        }
        return arg;
      }

inputs:
  input_vcf: File
  include_expression: { type: 'string?', doc: "See bcftools docs for valid expression. Can't be used at the same time as exclude_expression",
    inputBinding: { position: 1, prefix: "--include", shellQuote: true} }
  threads: { type: 'int?', default: 4, inputBinding: {position: 1, prefix: "--threads"} }
  exclude_expression: { type: 'string?', doc: "See bcftools docs for valid expression. Can't be used at the same time as include_expression",
    inputBinding: { position: 1, prefix: "--exclude", shellQuote: true} }
  filter_expression: { type: 'string?', doc: "Add values from FILTER field to subset on",
    inputBinding: { position: 1, prefix: "-f"}}
  output_type: { type: [ 'null', {type: enum, name: output_type, symbols: [ "u", "b", "v", "z"]}],
    inputBinding: { position: 1, prefix: "-O"}, default: "z" }
  sample_name: { type: 'string?', doc: "csv string of samples if user wishes to apply filtering to and output specific samples"}
  output_basename: string
outputs:
  filtered_vcf:
    type: File
    outputBinding:
      glob: "*.{v,b}cf{,.gz}"
    secondaryFiles: ['.tbi?']
