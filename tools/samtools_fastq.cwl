cwlVersion: v1.2
class: CommandLineTool
id: align2fastq
label: "Samtools bam/cram-to-fastq"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cores)
    ramMin: $(inputs.ram * 1000)
    https://platform.illumina.com/rdf/ica/resources:tier: economy

baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      set -eo pipefail

      ${
          if(inputs.is_paired_end){
              var command = "samtools sort -m 1G -n -O SAM -@ " + inputs.cores
              if (inputs.cram_reference){
                command += " --reference " + inputs.cram_reference.path
              }
              command += " " + inputs.input_reads_1.path + " | samtools fastq -c 2 -1 " + inputs.SampleID + ".converted_1.fastq.gz -2 " + inputs.SampleID + ".converted_2.fastq.gz -@ " + inputs.cores+ " -"
              return command
          }
          else {
              var command = "samtools sort -m 1G -n -O SAM -@ " + inputs.cores
              if (inputs.cram_reference){
                command += " --reference " + inputs.cram_reference.path
              }
              command += " " + inputs.input_reads_1.path + " | samtools fastq -c 2 -@ " + inputs.cores + " - > " + inputs.SampleID + ".converted_1.fastq && bgzip " + inputs.SampleID + ".converted_1.fastq"
              return command
          }
      }

inputs:
  input_reads_1: {type: File, doc: "Input alignment file"}
  SampleID: string
  cores: { type: 'int?', default: 16 } 
  ram: { type: 'int?', default: 32, doc: "GB of RAM to allocate to this task" } 
  is_paired_end: { type: boolean, doc: "Is the input_reads_1 file paired end?" }
  cram_reference: { type: 'File?', secondaryFiles: [.fai], doc: "If input align is cram and you are uncertain all contigs are registered at http://www.ebi.ac.uk/ena/cram/md5/, provide here" }

outputs:
  fq1:
    type: File
    outputBinding:
      glob: '*.converted_1.*'
  fq2:
    type: 'File?'
    outputBinding:
      glob: '*.converted_2.*'
