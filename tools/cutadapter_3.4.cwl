cwlVersion: v1.2
class: CommandLineTool
id: cutadapter_3.4
label: "Cutadapt v3.4 Trim"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/cutadapt:3.4'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 16000

baseCommand: [cutadapt, -j 8
, -m 20, --quality-base=33, -q 20]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -o TRIMMED.$(inputs.readFilesIn1.basename)
  - position: 3
    shellQuote: false
    valueFrom: >-
      ${
        var arg = "";
        if (inputs.r2_adapter && inputs.readFilesIn2){
          arg = " -p TRIMMED." + inputs.readFilesIn2.basename;
        }
        return arg;
      }
  - position: 5
    shellQuote: false
    valueFrom: >-
      > $(inputs.sample_name).cutadapt_results.txt

inputs:
  r1_adapter: { type: string, doc: "read1 adapter sequence", inputBinding: {prefix: "-a", position: 1} }
  r2_adapter: { type: 'string?', doc: "read2 adapter sequence, if paired", inputBinding: {prefix: "-A", position: 2} }
  readFilesIn1: { type: File, doc: "read1 fastq file", inputBinding: {position: 4} }
  readFilesIn2: { type: 'File?', doc: "read2 fastq file, if paired", inputBinding: {position: 4} }
  sample_name: string

outputs:
  trimmedReadsR1:
    type: File
    outputBinding:
      glob: $("*TRIMMED." + inputs.readFilesIn1.basename)
  trimmedReadsR2:
    type: 'File?'
    outputBinding:
      glob: ${
          if (inputs.readFilesIn2){
              return "*TRIMMED." + inputs.readFilesIn2.basename;
          }
          else{
              return null;
          }
        }
  cutadapt_stats:
    type: File
    outputBinding:
      glob: $(inputs.sample_name).cutadapt_results.txt
