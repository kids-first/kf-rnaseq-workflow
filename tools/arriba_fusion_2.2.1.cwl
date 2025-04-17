cwlVersion: v1.2
class: CommandLineTool
id: arriba_fusion_2-2-1
label: "Arriba Fusion Caller v2.2.1"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/arriba:2.2.1'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: ${ return inputs.memory * 1000 }
    https://platform.illumina.com/rdf/ica/resources:tier: economy

baseCommand: [/arriba_v2.2.1/arriba]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      -o $(inputs.outFileNamePrefix).arriba_2.2.1.fusions.tsv
      -O $(inputs.outFileNamePrefix).arriba_2.2.1.discarded_fusions.tsv

inputs:
  genome_aligned_bam: { type: File, doc: "STAR-aligned, coordinate sorted bam file",
  secondaryFiles: [ { pattern: ".bai", required: false },  { pattern: "^.bai", required: false } ],
  inputBinding: { prefix: '-x', position: 1 } }
  memory: { type: 'int?', doc: "Mem intensive tool. Set in GB", default: 64 }
  reference_fasta: { type: File, doc: "Fasta reference file used for alignment", inputBinding: { prefix: '-a', position: 1 } }
  gtf_anno: { type: File, doc: "GTF file used for alignment indexing",  inputBinding: { prefix: '-g', position: 1 } }
  outFileNamePrefix: string
  arriba_strand_flag: {  type: ['null', {type: enum, name: arriba_strand_flag, symbols: ["auto", "yes", "no", "reverse"]}],
  default: 'auto', doc: "input strandedness flag",
  inputBinding: { prefix: '-s', position: 1 } }
  blacklist:
    type:
      - 'null'
      - type: enum
        name: blacklist
        symbols:
          - 'hg38_GRCh38'
          - 'hg19_hs37d5_GRCh37'
          - 'mm10_GRCm38'
          - 'mm39_GRCm39'
    default: 'hg38_GRCh38'
    doc: "Path to built-in blacklist to use"
    inputBinding:
      position: 1
      prefix: '-b'
      valueFrom: >-
        /arriba_v2.2.1/database/blacklist_$(self)_v2.2.1.tsv.gz
      shellQuote: false
  known:
    type:
      - 'null'
      - type: enum
        name: known
        symbols:
          - 'hg38_GRCh38'
          - 'hg19_hs37d5_GRCh37'
          - 'mm10_GRCm38'
          - 'mm39_GRCm39'
    default: 'hg38_GRCh38'
    doc: "Path to built-in known/recurrent fusions"
    inputBinding:
      position: 1
      prefix: '-k'
      valueFrom: >-
        /arriba_v2.2.1/database/known_fusions_$(self)_v2.2.1.tsv.gz
      shellQuote: false
  tags:
    type:
      - 'null'
      - type: enum
        name: tags
        symbols:
          - 'hg38_GRCh38'
          - 'hg19_hs37d5_GRCh37'
          - 'mm10_GRCm38'
          - 'mm39_GRCm39'
    default: 'hg38_GRCh38'
    doc: "Path to built-in tag files. Typically same as known input"
    inputBinding:
      position: 1
      prefix: '-t'
      valueFrom: >-
        /arriba_v2.2.1/database/known_fusions_$(self)_v2.2.1.tsv.gz
      shellQuote: false
  protein_domains:
    type:
      - 'null'
      - type: enum
        name: protein_domains
        symbols:
          - 'hg38_GRCh38'
          - 'hg19_hs37d5_GRCh37'
          - 'mm10_GRCm38'
          - 'mm39_GRCm39'
    default: 'hg38_GRCh38'
    doc: "Path to built-in protein domain annotation"
    inputBinding:
      position: 1
      prefix: '-p'
      valueFrom: >-
        /arriba_v2.2.1/database/protein_domains_$(self)_v2.2.1.gff3
      shellQuote: false

outputs:
  arriba_fusions:
    type: File
    outputBinding:
      glob: "$(inputs.outFileNamePrefix).arriba_2.2.1.fusions.tsv"
