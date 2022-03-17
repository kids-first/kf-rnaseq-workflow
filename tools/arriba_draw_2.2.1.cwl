cwlVersion: v1.2
class: CommandLineTool
id: arriba_draw_2-2-1
label: "Arriba Draw v2.2.1"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/arriba:2.2.1'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: ${ return inputs.memory * 1000 }

baseCommand: [/arriba_v2.2.1/draw_fusions.R]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      --output=$(inputs.fusions.nameroot).pdf

inputs:
  fusions: { type: File, doc: "Fusion calls from Arriba OR STAR-Fusion", inputBinding: { prefix: '--fusions=', separate: false, position: 1 } }
  genome_aligned_bam: { type: File, doc: "STAR-aligned, coordinate sorted bam file",
  secondaryFiles: [ { pattern: ".bai", required: false },  { pattern: "^.bai", required: false } ],
  inputBinding: { prefix: '--alignments=', separate: false, position: 1 } }
  memory: { type: 'int?', doc: "Mem intensive tool. Set in GB", default: 16 }
  gtf_anno: { type: File, doc: "GTF file used for alignment indexing",  inputBinding: { prefix: '--annotation=', separate: false, position: 1 } }
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
      prefix: '--proteinDomains='
      separate: false
      valueFrom: >-
        /arriba_v2.2.1/database/protein_domains_$(self)_v2.2.1.gff3
      shellQuote: false
  cytobands:
    type:
      - 'null'
      - type: enum
        name: cytobands
        symbols:
          - 'hg38_GRCh38'
          - 'hg19_hs37d5_GRCh37'
          - 'mm10_GRCm38'
          - 'mm39_GRCm39'
    default: 'hg38_GRCh38'
    doc: "Path to built-in coordinates of the Giemsa staining bands to draw ideograms"
    inputBinding:
      position: 1
      prefix: '--cytobands='
      separate: false
      valueFrom: >-
        /arriba_v2.2.1/database/cytobands_$(self)_v2.2.1.tsv
      shellQuote: false

outputs:
  arriba_pdf:
    type: File
    outputBinding:
      glob: "$(inputs.fusions.nameroot).pdf"
