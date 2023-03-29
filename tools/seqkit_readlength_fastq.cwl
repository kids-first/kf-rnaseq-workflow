cwlVersion: v1.2
class: CommandLineTool
id: seqkit_readlength_fastq
doc: |-
  Given a fastq input:
  - Use seqkit + head to grab the first 1000000 records and convert that to a TSV
  - Grab the column that contains the read length from the TSV (column 2)
  - Get the unique read lengths as well as their counts
  - Sort those read lengths by count in descending order
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cores)
  - class: DockerRequirement
    dockerPull: pgc-images.sbgenomics.com/d3b-bixu/seqkit:2.3.1
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      seqkit fx2tab -nl $(inputs.input_fastq.path) | head -n 1000000 | cut -f 2 | sort | uniq -c | sort -nr -k1,1 > $(inputs.input_fastq.nameroot).fastq_readlength

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
