class: CommandLineTool
cwlVersion: v1.0
id: awk-junction-filtering
doc: "Takes JC.txt file output from rmats and removes calls that have junction counts (sum of fields: IJC_SAMPLE_1,SJC_SAMPLE_1,IJC_SAMPLE_2,SJC_SAMPLE_2) less than 10"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'ubuntu:20.04'
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.threads)
    https://platform.illumina.com/rdf/ica/resources:tier: economy
baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      awk -F'\t' 'NR==1 { for (i=1; i<=NF; i++) { f[$i] = i } print $line } $f["IJC_SAMPLE_1"]+$f["IJC_SAMPLE_2"]+$f["SJC_SAMPLE_1"]+$f["SJC_SAMPLE_2"] >= 10 {print $line}' $(inputs.input_jc_file.path) > ${ var arr = inputs.input_jc_file.basename.split('.'); arr.splice(-4,0,'filtered'); return arr.join('.') }
inputs:
  input_jc_file: { type: 'File', doc: "JC.txt file output form rmats" }
  threads: { type: 'int?', default: 1, doc: "The number of threads. The optimal number of threads should be equal to the number of CPU cores." }
  ram: { type: 'int?', default: 2, doc: "GB of RAM to allocate to this task." }
outputs:
  output: { type: 'File', outputBinding: { glob: '*filtered.*MATS.JC.txt' } }
