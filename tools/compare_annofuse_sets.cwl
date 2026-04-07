cwlVersion: v1.2
class: CommandLineTool
id: compare_annofuse_sets
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cpu)
    ramMin: $(inputs.ram * 1000)
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/danmiller/plot_gencode_rnaseq_python:3.13.0'
  - class: InitialWorkDirRequirement
    listing:
      - entryname: compare_sets.py
        entry:
          $include: ../scripts/compare_sets.py
baseCommand: []
arguments:
- position: 1
  shellQuote: false
  valueFrom: >-
    python compare_sets.py --dump-unique-dir unique
inputs:
  a_files: {type: 'File[]', inputBinding: {position: 2, prefix: "--file-a"}, doc: "One or more TSV files belonging to set A" }
  b_files: {type: 'File[]', inputBinding: {position: 2, prefix: "--file-b"}, doc: "One or more TSV files belonging to set B" }
  output_filename: {type: 'string', inputBinding: {position: 2, prefix: "--out"}, doc: "Output figure path (e.g., overlap.png, overlap.pdf, overlap.svg)" }
  a_label: {type: 'string?', inputBinding: {position: 2, prefix: "--label-a"}, doc: "Label to use for set A in the legend and filenames" }
  b_label: {type: 'string?', inputBinding: {position: 2, prefix: "--label-b"}, doc: "Label to use for set B in the legend and filenames" }
  title: {type: 'string?', inputBinding: {position: 2, prefix: "--title"}, doc: "Plot title" }
  tool: {type: 'string?', default: "annoFuse", inputBinding: {position: 2, prefix: "--tool"}, doc: "Tool in comparison" }
  percent: {type: 'boolean?', inputBinding: {position: 2, prefix: "--percent"}, doc: "Stack the bars by percentage rather than raw count" }
  cpu: { type: 'int?', default: 1, doc: "CPUs to allocate to this task" }
  ram: { type: 'int?', default: 8, doc: "GB of RAM to allocate to this task" } 
outputs:
  plot:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
  unique_fusions:
    type: File[]
    outputBinding:
      glob: 'unique/*'
