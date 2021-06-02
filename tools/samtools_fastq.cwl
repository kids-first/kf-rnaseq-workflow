cwlVersion: v1.0
class: CommandLineTool
id: bam2fastq
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.runThreadN)

baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail

      ${
          if(inputs.input_type == "PEBAM"){
              var command = "samtools sort -m 1G -n -O SAM -@ " + inputs.runThreadN + " " + inputs.input_reads_1.path + " | samtools fastq -c 2 -1 " + inputs.SampleID + ".converted_1.fastq.gz -2 " + inputs.SampleID + ".converted_2.fastq.gz -@ " + inputs.runThreadN+ " -"
              return command
          }
          else if(inputs.input_type == "SEBAM"){
              var command = "samtools sort -m 1G -n -O SAM -@ " + inputs.runThreadN + " " + inputs.input_reads_1.path + " | samtools fastq -c 2 -@ " + inputs.runThreadN + " - > " + inputs.SampleID + ".converted_1.fastq && bgzip " + inputs.SampleID + ".converted_1.fastq"
              return command
          }
          else if(inputs.input_type == "FASTQ" && inputs.input_reads_2 != null){
              var extr1 = (inputs.input_reads_1.nameext == '.gz' ? '.fastq.gz' : inputs.input_reads_1.nameext)
              var extr2 = (inputs.input_reads_2.nameext == '.gz' ? '.fastq.gz' : inputs.input_reads_2.nameext)
              var command =  "cp " + inputs.input_reads_1.path + " " + inputs.input_reads_1.nameroot + ".converted_1" + extr1 + " && cp " + inputs.input_reads_2.path + " " + inputs.input_reads_2.nameroot + ".converted_2" + extr2
              return command
          }
          else if(inputs.input_type == "FASTQ" && inputs.input_reads_2 == null){
              var extr1 = (inputs.input_reads_1.nameext == '.gz' ? '.fastq.gz' : inputs.input_reads_1.nameext)
              var command =  "cp " + inputs.input_reads_1.path + " " + inputs.input_reads_1.nameroot + ".converted_1" + extr1
              return command
          }
      }


inputs:
  input_reads_1: {type: File, doc: "For FASTQ input, please enter reads 1 here. For BAM input, please enter reads here."}
  input_reads_2: {type: 'File?', doc: "For FASTQ input, please enter reads 2 here. For BAM input, leave empty."}
  SampleID: string
  runThreadN: int
  input_type: {type: [{type: enum, name: input_type, symbols: ["PEBAM", "SEBAM", "FASTQ"]}], doc: "Please select one option for input file type, PEBAM (paired-end BAM), SEBAM (single-end BAM), or FASTQ."}


outputs:
  fq1:
    type: File
    outputBinding:
      glob: '*.converted_1.*'

  fq2:
    type: File?
    outputBinding:
      glob: '*.converted_2.*'
