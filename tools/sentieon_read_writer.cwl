cwlVersion: v1.2
class: CommandLineTool
id: sentieon-bam-to-cram
label: "Sentieon bam to cram"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/hdchen/sentieon:202503.01'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.threads)
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      export SENTIEON_LICENSE=10.5.64.221:8990
  - position: 2
    shellQuote: false
    valueFrom: >-
      && sentieon driver
  - position: 10
    shellQuote: false
    valueFrom: >-
      --algo ReadWriter
      $(inputs.input_bam.basename).cram


inputs:
  input_bam: { type: File, doc: "Input bam to convert", 
    inputBinding: { position: 3, prefix: "-i"}, secondaryFiles: [{"pattern": "^.bai", required: false}, {"pattern": ".bai", required: false}] }
  threads: { type: 'int?', default: 8,
    inputBinding: { position: 3, prefix: "-t"} }
  reference: { type: File, doc: "Reference FASTA used",
    inputBinding: { position: 3, prefix: "-r" }, secondaryFiles: [.fai] }
 

outputs:
  cram: { type: File, doc: "Converted CRAM", outputBinding: {glob: '*.cram'}, secondaryFiles: [.crai]  }

$namespaces:
  sbg: https://sevenbridges.com
