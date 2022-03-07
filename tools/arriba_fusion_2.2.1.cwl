cwlVersion: v1.2
class: CommandLineTool
id: arriba_fusion_2-2-1
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'migbro/arriba:2.2.1'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: ${ return inputs.memory * 1000 }

baseCommand: [/arriba_v2.2.1/arriba]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      -o $(inputs.outFileNamePrefix).arriba_2.2.1.fusions.tsv
      -O $(inputs.outFileNamePrefix).arriba_2.2.1.discarded_fusions.tsv
      -b /arriba_v2.2.1/database/blacklist_hg38_GRCh38_v2.2.1.tsv.gz
      -k /arriba_v2.2.1/database/known_fusions_hg38_GRCh38_v2.2.1.tsv.gz
      -t /arriba_v2.2.1/database/known_fusions_hg38_GRCh38_v2.2.1.tsv.gz
      -p /arriba_v2.2.1/database/protein_domains_hg38_GRCh38_v2.2.1.gff3
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
      /arriba_v2.2.1/draw_fusions.R
      --annotation=$(inputs.gtf_anno.path)
      --fusions=$(inputs.outFileNamePrefix).arriba_2.2.1.fusions.tsv
      --alignments=$(inputs.genome_aligned_bam.path)
      --cytobands=/arriba_v2.2.1/database/cytobands_hg38_GRCh38_v2.2.1.tsv
      --proteinDomains=/arriba_v2.2.1/database/protein_domains_hg38_GRCh38_v2.2.1.gff3
      --output=$(inputs.outFileNamePrefix).arriba_2.2.1.fusions.pdf

inputs:
  genome_aligned_bam: { type: File, doc: "STAR-aligned, coordinate sorted bam file",
  secondaryFiles: [ { pattern: ".bai", required: false },  { pattern: "^.bai", required: false } ],
  inputBinding: { prefix: '-x', position: 1 } }
  memory: { type: 'int?', doc: "Mem intensive tool. Set in GB", default: 64 }
  reference_fasta: { type: File, doc: "Fasta reference file used for alignment", inputBinding: { prefix: '-a', position: 1 } }
  gtf_anno: { type: File, doc: "GTF file used for alignment indexing",  inputBinding: { prefix: '-g', position: 1 } }
  outFileNamePrefix: string
  arriba_strand_flag: ['null', string]

outputs:
  arriba_fusions:
    type: File
    outputBinding:
      glob: "$(inputs.outFileNamePrefix).arriba_2.2.1.fusions.tsv"
  arriba_pdf:
    type: File
    outputBinding:
      glob: "$(inputs.outFileNamePrefix).arriba_2.2.1.fusions.pdf"
