cwlVersion: v1.2
class: CommandLineTool
id: bam2fastq
label: "Samtools bam-to-fastq"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cores)
    ramMin: $(inputs.cores * 1000)

baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail

      ${
          if(inputs.input_type == "PEBAM"){
              var command = "samtools sort -m 1G -n -O SAM -@ " + inputs.cores + " " + inputs.input_reads_1.path + " | samtools fastq -c 2 -1 " + inputs.SampleID + ".converted_1.fastq.gz -2 " + inputs.SampleID + ".converted_2.fastq.gz -@ " + inputs.cores+ " -"
              return command
          }
          else if(inputs.input_type == "SEBAM"){
              var command = "samtools sort -m 1G -n -O SAM -@ " + inputs.cores + " " + inputs.input_reads_1.path + " | samtools fastq -c 2 -@ " + inputs.cores + " - > " + inputs.SampleID + ".converted_1.fastq && bgzip " + inputs.SampleID + ".converted_1.fastq"
              return command
          }
      }

inputs:
  input_reads_1: {type: File, doc: "For FASTQ input, please enter reads 1 here. For BAM input, please enter reads here."}
  SampleID: string
  cores: { type: 'int?', default: 36 } 
  input_type: {type: [{type: enum, name: input_type, symbols: ["PEBAM", "SEBAM"]}], doc: "Please select one option for input file type, PEBAM (paired-end BAM), SEBAM (single-end BAM), or FASTQ."}


outputs:
  fq1:
    type: File
    outputBinding:
      glob: '*.converted_1.*'

  fq2:
    type: 'File?'
    outputBinding:
      glob: '*.converted_2.*'
