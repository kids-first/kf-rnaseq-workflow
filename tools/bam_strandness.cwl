cwlVersion: v1.2
class: CommandLineTool
id: bam_strandness
doc: |-
  Converts the input BAM into a SAM file of n_reads * 2 length. Then it runs check_strandedness on that tiny SAM file.
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cores)
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/stranded:1.0.0'
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
