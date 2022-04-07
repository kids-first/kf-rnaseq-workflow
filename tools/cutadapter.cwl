cwlVersion: v1.2
class: CommandLineTool
id: cutadapter
label: "Cutadapt Trim"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/cutadapt:latest'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 16
    ramMin: 24000

baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      ${
        if (inputs.r1_adapter == null){
          var cmd = "cp " + inputs.readFilesIn1.path + " ./UNTRIMMED." + inputs.readFilesIn1.basename
          if (inputs.readFilesIn2 != null){
            cmd += ";cp " + inputs.readFilesIn2.path + " ./UNTRIMMED." + inputs.readFilesIn2.basename;
        }
        return cmd;
        }
        else{
          var cmd = "cutadapt -j 16 -m 20 --quality-base=33 -q 20 -a " + inputs.r1_adapter;
          if (inputs.r2_adapter && inputs.readFilesIn2){
            cmd += " -A " + inputs.r2_adapter + " -p TRIMMED." + inputs.readFilesIn2.basename;
          }
          cmd += " -o TRIMMED." + inputs.readFilesIn1.basename + " " + inputs.readFilesIn1.path + " ";
          if (inputs.r2_adapter && inputs.readFilesIn2){
            cmd += inputs.readFilesIn2.path
          }
          cmd += " > " + inputs.sample_name + ".cutadapt_results.txt";
          return cmd;
        }
      }

inputs:
  r1_adapter: ['null', string]
  r2_adapter: ['null', string]
  readFilesIn1: File
  readFilesIn2: ['null', File]
  sample_name: string

outputs:
  trimmedReadsR1:
    type: File
    outputBinding:
      glob: $("*TRIMMED." + inputs.readFilesIn1.basename)
  trimmedReadsR2:
    type: ['null', File]
    outputBinding:
      glob: ${
          if (inputs.readFilesIn2){
              return "*TRIMMED." + inputs.readFilesIn2.basename
          }
          else{
              return "placeholder"
          }
        }
  cutadapt_stats:
    type: ['null', File]
    outputBinding:
      glob: $(inputs.sample_name).cutadapt_results.txt
