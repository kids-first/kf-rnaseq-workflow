cwlVersion: v1.2
class: CommandLineTool
id: seqkit_readlength_fastq
doc: |-
  seqkit fx2tab -nl NA18152.fastq
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cores)
  - class: DockerRequirement
    dockerPull: zqs891011/seqkit:2.3.1
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      seqkit fx2tab -nl $(inputs.input_fastq.path) | head -n 1000000 | cut -f 2 | sort | uniq -c > $(inputs.input_fastq.nameroot).fastq_readlength

inputs:
  input_fastq: { type: File, doc: "Input fastq file"}

outputs:
  output:
    type: File
    outputBinding:
      glob: "*.fastq_readlength"
  top_readlength:
    type: int
    outputBinding:
      glob: "*.fastq_readlength"
      loadContents: true
      outputEval: |
        ${
          var rows = self[0].contents.split(/\r?\n/).slice(0,-1);
          return rows[0].split(/\s/).pop();
        }
  variable_readlength:
    type: boolean
    outputBinding:
      glob: "*.fastq_readlength"
      loadContents: true
      outputEval: |
        ${
          var rows = self[0].contents.split(/\r?\n/).slice(0,-1);
          return rows.length > 1;
        }
