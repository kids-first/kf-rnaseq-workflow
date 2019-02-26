cwlVersion: v1.0
class: CommandLineTool
id: arriba_fusion
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/arriba:latest'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 32000

baseCommand: [/arriba_v1.0.1/arriba]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -x $(inputs.genome_aligned_bam.path)
      -a $(inputs.reference_fasta.path)
      -g $(inputs.gtf_anno.path)
      -o $(inputs.outFileNamePrefix).arriba.fusions.tsv
      -O $(inputs.outFileNamePrefix).arriba.discarded_fusions.tsv
      -b /arriba_v1.0.1/database/blacklist_hg38_GRCh38_2018-04-04.tsv.gz
      -T
      -P
      ${
        if(inputs.arriba_strand_flag == null){
          return "-s auto";
        }
        else{
          return "-s " + inputs.arriba_strand_flag;
        }
      }

inputs:
  genome_aligned_bam: File
  reference_fasta: File
  gtf_anno: File
  outFileNamePrefix: string
  arriba_strand_flag: ['null', string]

outputs:
  arriba_fusions:
    type: File
    outputBinding:
      glob: "$(inputs.outFileNamePrefix).arriba.fusions.tsv"
