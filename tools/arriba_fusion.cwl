cwlVersion: v1.2
class: CommandLineTool
id: arriba_fusion
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/arriba:1.1.0'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 64000

baseCommand: [/arriba_v1.1.0/arriba]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      -x $(inputs.genome_aligned_bam.path)
      -a $(inputs.reference_fasta.path)
      -g $(inputs.gtf_anno.path)
      -o $(inputs.outFileNamePrefix).arriba.fusions.tsv
      -O $(inputs.outFileNamePrefix).arriba.discarded_fusions.tsv
      -b /arriba_v1.1.0/database/blacklist_hg38_GRCh38_2018-11-04.tsv.gz
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
  - position: 2
    shellQuote: false
    valueFrom: >-
      &&
      /arriba_v1.1.0/draw_fusions.R
      --annotation=$(inputs.gtf_anno.path)
      --fusions=$(inputs.outFileNamePrefix).arriba.fusions.tsv
      --alignments=$(inputs.genome_aligned_bam.path)
      --cytobands=/arriba_v1.1.0/database/cytobands_hg38_GRCh38_2018-02-23.tsv
      --proteinDomains=/arriba_v1.1.0/database/protein_domains_hg38_GRCh38_2018-03-06.gff3
      --output=$(inputs.outFileNamePrefix).arriba.fusions.pdf

inputs:
  genome_aligned_bam: { type: File, secondaryFiles: [ { pattern: ".bai", required: false },  { pattern: "^.bai", required: false } ] }
  genome_aligned_bai: 'File?'
  chimeric_sam_out: { type: 'File?', doc: "If older version of STAR used, separate sam/bam file might exist", inputBinding: { prefix: '-c', position: 1 } }
  reference_fasta: File
  gtf_anno: File
  outFileNamePrefix: string
  arriba_strand_flag: ['null', string]

outputs:
  arriba_fusions:
    type: File
    outputBinding:
      glob: "$(inputs.outFileNamePrefix).arriba.fusions.tsv"
  arriba_pdf:
    type: File
    outputBinding:
      glob: "$(inputs.outFileNamePrefix).arriba.fusions.pdf"
