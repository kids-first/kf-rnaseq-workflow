cwlVersion: v1.2
class: CommandLineTool
id: corr_heatmap
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cpu)
    ramMin: $(inputs.ram * 1000)
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/danmiller/plot_gencode_rnaseq_r:4.5.0'
  - class: InitialWorkDirRequirement
    listing:
      - entryname: corr_heatmap.R
        entry:
          $include: ../scripts/corr_heatmap.R
baseCommand: []
arguments:
- position: 1
  shellQuote: false
  valueFrom: >-
    Rscript corr_heatmap.R
inputs:
  gene_list: {type: 'File', inputBinding: {position: 2, prefix: "--gene_list_path" }, doc: "TSV containing ENSG names on which the TSNE will be performed" }
  a_files: {type: 'File[]', inputBinding: {position: 2, prefix: "--files_setA"}, doc: "One or more TSV files belonging to set A" }
  b_files: {type: 'File[]', inputBinding: {position: 2, prefix: "--files_setB"}, doc: "One or more TSV files belonging to set B" }
  a_label: {type: 'string?', inputBinding: {position: 2, prefix: "--setA_name"}, doc: "Label to use for set A in the legend and filenames" }
  b_label: {type: 'string?', inputBinding: {position: 2, prefix: "--setB_name"}, doc: "Label to use for set B in the legend and filenames" }
  output_basename: {type: 'string?', inputBinding: {position: 2, prefix: "--output_basename"}, doc: "Basename for output files (no extension). If omitted, defaults to <setA_name>v<setB_name>." }
  cpu: { type: 'int?', default: 1, doc: "CPUs to allocate to this task" }
  ram: { type: 'int?', default: 8, doc: "GB of RAM to allocate to this task" } 
outputs:
  plot:
    type: File
    outputBinding:
      glob: '*.sample_correlation_heatmap.pdf'
  csv:
    type: File
    outputBinding:
      glob: '*.sample_correlation_heatmap.csv'

