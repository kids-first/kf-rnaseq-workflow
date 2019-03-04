cwlVersion: v1.0
class: CommandLineTool
id: kallisto
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/uros_sipetic/kallisto:0.43.1'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 10000

baseCommand: [kallisto]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      ${
        var cmd = "quant -i " + inputs.transcript_idx.path + " -o output --fusion -b 10 -t 8";
        if (inputs.strand != null || inputs.strand != "default"){
          cmd += " --" + inputs.strand;
        }
        cmd += " " + inputs.reads1.path + " " + inputs.reads2.path;
        return cmd;
      }

      mv output/abundance.tsv $(inputs.SampleID).kallisto.abundance.tsv &&
      gzip $(inputs.SampleID).kallisto.abundance.tsv &&
      mv output/fusion.txt $(inputs.SampleID).kallisto.fusion.txt

inputs:
  transcript_idx: File
  strand: {type: ['null', string], doc: "input none if unstranded, otherwise rf-stranded or fr-stranded"}
  reads1: File
  reads2: File
  SampleID: string

outputs:
  abundance_out:
    type: File
    outputBinding:
      glob: '*.abundance.tsv.gz'

  fusion_out:
    type: File
    outputBinding:
      glob: '*.fusion.txt'
