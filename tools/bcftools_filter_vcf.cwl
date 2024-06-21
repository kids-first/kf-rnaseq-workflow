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
  - class: InitialWorkDirRequirement
    listing:
      - entryname: "run_filter.sh"
        entry: |
          #!/usr/bin/env bash
          set -xeo pipefail

          command=cat
          if [[ $(inputs.sample_name) != null ]]
          then
            command="bcftools view --threads $(inputs.threads) -s $(inputs.sample_name)"
          fi
          $command $(inputs.input_vcf.path) \
            | bcftools view \
              --threads $(inputs.threads) \
              -O $(inputs.output_type) \
              -o $(inputs.output_basename).bcf_filtered.$(inputs.output_type == "v" ? "vcf" : inputs.output_type == "z" ? "vcf.gz" : inputs.output_type == "b" ? "bcf.gz" : "bcf") \
              $(inputs.exclude_expression == null ? "" : "--exclude " + "'" + inputs.exclude_expression + "'") \
              $(inputs.include_expression == null ? "" : "--include " + "'" + inputs.include_expression + "'") \
              $(inputs.filter_expression == null ? "" : "-f " + inputs.filter_expression)
          if [[ $(inputs.output_type) == z ]]
          then
            tabix $(inputs.output_basename).bcf_filtered.vcf.gz
          fi

baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      /bin/bash run_filter.sh

inputs:
  input_vcf: File
  include_expression: { type: 'string?', doc: "See bcftools docs for valid expression. Can't be used at the same time as exclude_expression. Use double quotes when a string needs to be quoted"}
  threads: { type: 'int?', default: 4 }
  exclude_expression: { type: 'string?', doc: "See bcftools docs for valid expression. Can't be used at the same time as include_expression.  Use double quotes when a string needs to be quoted"}
  filter_expression: { type: 'string?', doc: "Add values from FILTER field to subset on"}
  output_type: { type: [ 'null', {type: enum, name: output_type, symbols: [ "u", "b", "v", "z"]}], default: "z"}
  sample_name: { type: 'string?', doc: "csv string of samples if user wishes to apply filtering to and output specific samples"}
  output_basename: string
outputs:
  filtered_vcf:
    type: File
    outputBinding:
      glob: "*.{v,b}cf{,.gz}"
    secondaryFiles: ['.tbi?']
  debug_run_filter_sh:
    type: File
    outputBinding:
      glob: "run_filter.sh"
