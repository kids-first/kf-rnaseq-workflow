cwlVersion: v1.0
class: ExpressionTool
id: expression_strand_params
requirements:
  - class: InlineJavascriptRequirement

inputs:
  wf_strand_param: {type: [{type: enum, name: wf_strand_param, symbols: ["default", "rf-stranded", "fr-stranded"]}], doc: "use 'default' for unstranded/auto, rf_stranded if read1 in the fastq read pairs is reverse complement to the transcript, fr-stranded if read1 same sense as transcript"}

outputs:
  rsem_std:
    type: string
  kallisto_std:
    type: string
  rnaseqc_std:
    type: string
  arriba_std:
    type: string

expression:
  "${
      var parse_dict = {
          'default': {'rsem_std': 'none', 'kallisto_std': 'default', 'rnaseqc_std': 'default', 'arriba_std': 'auto'},
          'rf-stranded': {'rsem_std': 'reverse', 'kallisto_std': 'rf-stranded', 'rnaseqc_std': 'rf', 'arriba_std': 'reverse'},
          'fr-stranded': {'rsem_std': 'forward', 'kallisto_std': 'fr-stranded', 'rnaseqc_std': 'fr', 'arriba_std': 'yes'}
          };
      return parse_dict[inputs.wf_strand_param];
  }"
