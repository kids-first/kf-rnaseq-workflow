cwlVersion: v1.0
class: CommandLineTool
id: rsem-calculate-expression
label: "RSEM v1.3.1 Calulate Expression"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/uros_sipetic/rsem:1.3.1'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.num_threads)
    ramMin: 24000

baseCommand: [tar]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      -zxf $(inputs.genomeDir.path)
  - position: 1
    shellQuote: false
    valueFrom: >-
      && rsem-calculate-expression
      --no-bam-output
      --alignments
  - position: 3
    shellQuote: false
    valueFrom: >-
      ./$(inputs.genomeDir.nameroot.replace(".tar", ""))/$(inputs.genomeDir.nameroot.replace(".tar", ""))
      $(inputs.outFileNamePrefix).rsem
  - position: 4
    shellQuote: false
    valueFrom: >-
      && gzip *results

inputs:
  paired-end: {type: 'boolean?', doc: "If input is paired-end, add this flag", default: true, inputBinding: {position: 1, prefix: '--paired-end'} }
  num_threads: { type: 'int?', doc: "Num threads to use", default: 16, inputBinding: { position: 1, prefix: '--num-threads'} }
  strandedness: { type: [ 'null', {type: enum, name: strandedness, symbols: ["none", "forward", "reverse"]}], default: "none",
  doc: "'none' refers to non-strand-specific protocols. 'forward' means all (upstream) reads are derived from the forward strand. 'reverse' means all (upstream) reads
  are derived from the reverse strand",
  inputBinding: { position: 1, prefix: '--strandedness', shellQuote: false } }
  append_names: { type: 'boolean?', doc: "If available, append gene/tx name to gene/tx id", default: true, inputBinding: {position: 1, prefix: '--append-names'} }
  estimate_rspd: { type: 'boolean?', doc: "Set this option if you want to estimate the read start position distribution (RSPD) from data", default: false, inputBinding: {position: 1, prefix: '--estimate-rspd'} }
  bam: { type: File, doc: "Aligned transcriptome bam", inputBinding: { position: 2 } }
  genomeDir: { type: File, doc: "RSEM reference tar ball" }
  outFileNamePrefix: {type: string, doc: "String to prepend output file names with" }

outputs:
  gene_out:
    type: File
    outputBinding: 
      glob: '*genes.results.gz'
  isoform_out:
    type: File
    outputBinding:
      glob: '*isoforms.results.gz'
