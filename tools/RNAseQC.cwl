cwlVersion: v1.0
class: CommandLineTool
id: rnaseqc
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'gcr.io/broad-cga-aarong-gtex/rnaseqc:latest'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 10000

baseCommand: [rnaseqc]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.collapsed_gtf.path)
      $(inputs.Aligned_sorted_bam.path)
      output/
      ${
        var cmd = "--legacy";
        if (inputs.strand != null or inputs.strand != "default"){
          cmd += " --stranded=" + inputs.strand;
        }
        return cmd;
      }

inputs:
  Aligned_sorted_bam: File
  collapsed_gtf: File
  strand: {type: ['null', string]}

outputs:
  Metrics:
    type: File
    outputBinding:
      glob: 'output/*.metrics.tsv'
  Gene_TPM:
    type: File
    outputBinding:
      glob: 'output/*.gene_tpm.gct'
  Gene_count:
    type: File
    outputBinding:
      glob: 'output/*.gene_reads.gct'
  Exon_count:
    type: File
    outputBinding:
      glob: 'output/*.exon_reads.gct'
