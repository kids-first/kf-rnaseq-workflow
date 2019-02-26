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
  rsem_std: double
  kallisto_std: string
  rnaseqc_std: string
  arriba_std: string

expression:
  "${
      var strand = inputs.strand;
      if (inputs.strand == null){
        strand = 'default';
      }
      var parse_dict = {
          'default': {'rsem_std': null, 'kallisto_std': null, 'rnaseqc_std': null, 'arriba_std': null},
          'rf-stranded': {'rsem_std': 0, 'kallisto_std': 'rf-stranded', 'rnaseqc_std': 'rf', 'arriba_std': 'reverse'},
          'fr-stranded': {'rsem_std': 1, 'kallisto_std': 'fr-stranded', 'rnaseqc_std': 'fr', 'arriba_std': 'yes'}
          };

        return parse_dict[inputs.strand];

  }"
