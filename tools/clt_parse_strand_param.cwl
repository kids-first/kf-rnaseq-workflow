cwlVersion: v1.2
class: CommandLineTool
id: clt_strand_params
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement

baseCommand: [echo, done]

inputs:
  wf_strand_param: {type: ['null', {type: enum, name: wf_strand_param, symbols: ["default", "rf-stranded", "fr-stranded"]}], default: "default", doc: "use 'default' for unstranded/auto, rf_stranded if read1 in the fastq read pairs is reverse complement to the transcript, fr-stranded if read1 same sense as transcript"}

outputs:
  rsem_std:
    type:
    - type: enum
      name: rsem_std
      symbols: ["none", "forward", "reverse"]
    outputBinding:
      outputEval: |
        ${
          var parse_dict = {
            'default': 'none',
            'rf-stranded': 'reverse',
            'fr-stranded': 'forward'
          }
          return parse_dict[inputs.wf_strand_param]
        }
  kallisto_std:
    type: string
    outputBinding:
      outputEval: |
        ${
          var parse_dict = {
            'default': 'default',
            'rf-stranded': 'rf-stranded',
            'fr-stranded': 'fr-stranded'
          }
          return parse_dict[inputs.wf_strand_param]
        }
  rnaseqc_std:
    type:
    - 'null'
    - type: enum
      name: rnaseqc_std
      symbols: ["rf", "fr"]
    outputBinding:
      outputEval: |
        ${
          var parse_dict = {
            'default': null,
            'rf-stranded': 'rf',
            'fr-stranded': 'fr'
          }
          return parse_dict[inputs.wf_strand_param]
        }
  arriba_std:
    type:
    - type: enum
      name: arriba_std
      symbols: ["auto", "reverse", "yes"]
    outputBinding:
      outputEval: |
        ${
          var parse_dict = {
            'default': 'auto',
            'rf-stranded': 'reverse',
            'fr-stranded': 'yes'
          }
          return parse_dict[inputs.wf_strand_param]
        }
  rmats_std:
    type:
    - 'null'
    - type: enum
      name: rmats_std
      symbols: ["fr-firststrand", "fr-secondstrand", "fr-unstranded"]
    outputBinding:
      outputEval: |
        ${
          var parse_dict = {
            "rf-stranded": "fr-firststrand",
            "fr-stranded": "fr-secondstrand",
            "default": "fr-unstranded"
          }
          return parse_dict[inputs.wf_strand_param]
        }
