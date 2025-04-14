class: CommandLineTool
cwlVersion: v1.2
id: samtools_split
doc: |-
  This tool splits the input bam input read group bams if it has more than one readgroup.
  Programs run in this tool:
    - samtools head | grep
    - samtools split
  Using samtools view and grep count the header lines starting with @RG. If that number is
  not one, split the bam file into read group bams using samtools.
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'staphb/samtools:1.20'
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.cpu)
    https://platform.illumina.com/rdf/ica/resources:tier: economy
    https://platform.illumina.com/rdf/ica/resources:type: hicpu
    https://platform.illumina.com/rdf/ica/resources:size: small
  - class: InitialWorkDirRequirement
    listing:
      - entryname: split_bam.sh
        entry: |-
          set -xeo pipefail
          RG_NUM=`samtools head $(inputs.input_reads.path) | grep -c ^@RG`
          if [ $RG_NUM != 1 ]; then
            samtools split -f '%*_%#.bam' -@ $(inputs.cpu) $(inputs.reference ? '--reference ' + inputs.reference.path : '') $(inputs.input_reads.path)
          fi
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      /bin/bash split_bam.sh
inputs:
  input_reads: { type: File, doc: "Input bam file" }
  reference: { type: 'File?', doc: "Reference fasta file" }
  ram: { type: 'int?', default: 16, doc: "GB of RAM to allocate to the task." }
  cpu: { type: 'int?', default: 32, doc: "Minimum reserved number of CPU cores for the task." }
outputs:
  bam_files:
    type: File[]
    outputBinding:
      glob: '*.bam'
      outputEval: |-
        ${
          if (self.length == 0) return [inputs.input_reads]
          else return self
        }
