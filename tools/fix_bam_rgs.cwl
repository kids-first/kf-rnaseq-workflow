cwlVersion: v1.2
class: CommandLineTool
id: fix_bam_rgs
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.cpu)
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/pysam:0.23.3--py313hd07c5dd_2
  - class: InitialWorkDirRequirement
    listing:
      - entryname: fix_bam_rgs.py
        entry:
          $include: ../scripts/fix_bam_rgs.py
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      python fix_bam_rgs.py -o $(inputs.output_basename).$(inputs.cram_reference != null ? "cram" : "bam") --stats-out $(inputs.output_basename).rg_stats.json

inputs:
  inbam: { type: 'File', inputBinding: {position: 2, prefix: "--input"}, doc: "BAM/CRAM to fix" }
  cram_reference: { type: 'File?', inputBinding: {position: 2, prefix: "--reference"}, doc: "Reference to decode CRAMs" }
  output_basename: { type: 'string' }
  rg_sm: { type: 'string?', inputBinding: {position: 2, prefix: "--sm"}, doc: "Value for @RG SM (default: inherit from input header, else SAMPLE)" }
  rg_pl: { type: 'string?', inputBinding: {position: 2, prefix: "--pl"}, doc: "Value for @RG PL (default: inherit from input header, else ILLUMINA)" }
  rg_lb: { type: 'string?', inputBinding: {position: 2, prefix: "--lb"}, doc: "Value for @RG LB (default: inherit from input header, else LIB1)" }
  single_pass: { type: 'boolean?', inputBinding: {position: 2, prefix: "--single-pass"}, doc: "Use single-pass RG discovery with bounded buffering. Disables stats output" } 
  extra_args: { type: 'string?', inputBinding: {position: 3}, doc: "Any extra args for this task." }
  cpu: { type: 'int?', default: 8, inputBinding: {position: 2, prefix: "-p" }, doc: "Num processing threads to use" }
  ram: { type: 'int?', doc: "Num GB memory to make available", default: 16 }
outputs:
  fixed_am: { type: 'File', outputBinding: { glob: "*am" }}
  rg_stats: { type: 'File?', outputBining: { glob: "*.rg_stats.json" }}
