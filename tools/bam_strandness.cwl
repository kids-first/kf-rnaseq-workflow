cwlVersion: v1.2
class: CommandLineTool
id: bam_strandness
doc: |-
  Converts the input BAM into a SAM file of n_reads * 2 length. Then it runs check_strandedness on that tiny SAM file.
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
      - |
        /**
        * Given the contents of a file and a string keyword. Return the post-colon value on the line containing the keyWord
        * For instance if the line contents were KeyWord:Value, the function will return Value.
        *
        * @param {String} fileContents - The contents of the file
        * @param {String} keyWord - The exact keyWord for which we will be checking
        * @return {String, null} the string after the keyWord and colon; null if no keyWord found
        */
        function returnKeyValue (fileContents, keyWord) {
          var rows = fileContents.split(/\r?\n/).slice(0,-1);
          for (var row in rows) {
            if (rows[row].search(keyWord) == 0) {
              return rows[row].split(':')[1];
            }
          }
          return null;
        }
  - class: ResourceRequirement
    coresMin: $(inputs.cores)
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/danmiller/stranded:1.1.2'
baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      samtools view -h $(inputs.input_bam.path) | head -n $(inputs.n_reads * 2) |
      $(inputs.paired_end ? "samtools sort -n -@ " + (inputs.cpu - 3) + " |" : "")
      samtools fastq -@ 2 -n $(inputs.paired_end ? "-1 1.fastq -2 2.fastq -0 /dev/null" : "-0 1.fastq") -s /dev/null
  - position: 10
    shellQuote: false
    prefix: "&&"
    valueFrom: >-
      check_strandedness --gtf $(inputs.annotation_gtf.path) --kallisto_index $(inputs.kallisto_idx.path) --reads_1 1.fastq $(inputs.paired_end ? "--reads_2 2.fastq" : "") --nreads $(inputs.n_reads) > $(inputs.input_bam.nameroot).bam.strandness

inputs:
  input_bam: {type: File, doc: "Input bam file"}
  annotation_gtf: {type: 'File', doc: "gtf file from `gtf_anno` is the primary gtf from gencode"}
  kallisto_idx: {type: 'File', doc: "Specialized index of a **transcriptome** fasta file for kallisto"}
  n_reads: {type: 'int?', doc: "number of reads to sample", default: 200000}
  paired_end: { type: 'boolean?', doc: "Set to true if reads are paired end." }
  cpu: { type: 'int?', default: 16, doc: "CPUs to allocate to this task" }
outputs:
  output:
    type: File
    outputBinding:
      glob: "*.bam.strandness"
  read_length_median:
    type: int?
    outputBinding:
      glob: "*.bam.strandness"
      loadContents: true
      outputEval: |
        $(parseInt(returnKeyValue(self[0].contents.trim(), "MedianReadLength")))
  read_length_stddev:
    type: float?
    outputBinding:
      glob: "*.bam.strandness"
      loadContents: true
      outputEval: |
        $(parseFloat(returnKeyValue(self[0].contents.trim(), "StddevReadLength")))
  strandedness:
    type:
      - 'null'
      - type: enum
        name: strandedness
        symbols:
          - "default"
          - "rf-stranded"
          - "fr-stranded"
    outputBinding:
      glob: "*.bam.strandness"
      loadContents: true
      outputEval: |
        ${
          var rows = self[0].contents.split(/\r?\n/).slice(0,-1);
          var lastword = rows.pop().split(/\s/).pop();
          switch(lastword) {
            case "FR/fr-stranded":
            case "FR/fr-secondstrand":
              return "fr-stranded";
            case "RF/rf-stranded":
            case "RF/fr-firststrand":
              return "rf-stranded";
            case "unstranded":
              return "default";
            default:
              return null;
          }
        }
