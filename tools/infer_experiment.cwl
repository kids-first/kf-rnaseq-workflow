cwlVersion: v1.2
class: CommandLineTool
id: infer_experiment 
doc: |-
  Determines whether a given BAM/SAM is paired or single end.
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cores)
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/danmiller/stranded:1.1.2'
baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      gtf2bed --gtf $(inputs.input_gtf.path) --bed $(inputs.input_gtf.basename.replace(/.gtf$/, '.bed'))
  - position: 11
    shellQuote: false
    prefix: "&&"
    valueFrom: >-
      infer_experiment.py -i $(inputs.input_reads.path) -r $(inputs.input_gtf.basename.replace(/.gtf$/, '.bed')) -s $(inputs.n_reads) > $(inputs.input_reads.basename.replace(/.(b|cr|s)am$/, '.infered_experiment.txt'))

inputs:
  input_reads: {type: File, doc: "Input BAM/SAM file"}
  input_gtf: {type: 'File', doc: "gtf file from `gtf_anno` is the primary gtf from gencode"}
  n_reads: {type: 'int?', doc: "number of reads to sample", default: 200000}
  cpu: { type: 'int?', default: 16, doc: "CPUs to allocate to this task" }
outputs:
  infered_experiment:
    type: File
    outputBinding:
      loadContents: true
      glob: "*.infered_experiment.txt"
  is_paired_end:
    type: boolean?
    outputBinding:
      loadContents: true
      glob: "*.infered_experiment.txt"
      outputEval: |
        ${
          var rows = self[0].contents.trim().split(/\r?\n/);
          for (var rowNum in rows) {
            if (rows[rowNum].search("PairEnd") != -1) { return true; }
            if (rows[rowNum].search("SingleEnd") != -1) { return false; }
          }
          return null;
        }
