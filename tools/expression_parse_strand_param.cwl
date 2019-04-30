cwlVersion: v1.0
class: ExpressionTool
id: expression_strand_params
requirements:
  - class: InlineJavascriptRequirement

inputs:
  strand:
    type: ['null', string]
    doc: "use 'default' for unstranded/auto, rf_stranded if read1 in the fastq read pairs is reverse complement to the transcript, fr-stranded if read1 same sense as transcript"

outputs:
  rsem_std: string
  kallisto_std: string
  rnaseqc_std: string
  arriba_std: string

expression:
  "${
      var strand = 'default';
      if (inputs.strand != null){
        strand = inputs.strand;
      }
      var parse_dict = {
          'default': {'rsem_std': 'none', 'kallisto_std': 'default', 'rnaseqc_std': 'default', 'arriba_std': 'auto'},
          'rf-stranded': {'rsem_std': 'reverse', 'kallisto_std': 'rf-stranded', 'rnaseqc_std': 'rf', 'arriba_std': 'reverse'},
          'fr-stranded': {'rsem_std': 'forward', 'kallisto_std': 'fr-stranded', 'rnaseqc_std': 'fr', 'arriba_std': 'yes'}
          };
      if (strand in parse_dict){
        return parse_dict[strand];
        
      }
      else{
        throw new Error(strand + ' is a not a valid strand param. Use one of default, rf-stranded, fr-stranded');
      }
  }"
