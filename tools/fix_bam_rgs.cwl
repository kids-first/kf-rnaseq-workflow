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
    dockerPull: pgc-images.sbgenomics.com/danmiller/pysam:0.1.0
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
      python fix_bam_rgs.py -o $(inputs.output_basename).bam

inputs:
  inbam: { type: 'File', inputBinding: {position: 2, prefix: "--input"}, doc: "BAM/CRAM to fix" }
  cram_reference: { type: 'File', inputBinding: {position: 2, prefix: "--reference"}, doc: "Reference to decode CRAMs" }
  output_basename: { type: 'string' }
  rg_sm: { type: 'string?', inputBinding: {position: 2, prefix: "--sm"}, doc: "Value for @RG SM (default: inherit from input header, else SAMPLE)" }
  rg_pl: { type: 'string?', inputBinding: {position: 2, prefix: "--pl"}, doc: "Value for @RG PL (default: inherit from input header, else SAMPLE)" }
  rg_lb: { type: 'string?', inputBinding: {position: 2, prefix: "--lb"}, doc: "Value for @RG LB (default: inherit from input header, else SAMPLE)" }
  cpu: { type: 'int?', default: 8, inputBinding: {position: 2, prefix: "--hts_threads" }, doc: "Num processing threads to use" }
  ram: { type: 'int?', doc: "Num GB memory to make available", default: 16 }
outputs:
  fixed_bam: { type: 'File', outputBinding: { glob: "*bam" }}
