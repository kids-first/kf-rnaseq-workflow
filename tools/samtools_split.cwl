class: CommandLineTool
cwlVersion: v1.2
id: samtools_split
doc: |-
  This tool splits the input bam input read group bams if it has more than one readgroup.
  Programs run in this tool:
    - samtools addreplacerg
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
  - class: InitialWorkDirRequirement
    listing:
      - entryname: split_bam.sh
        entry: |-
          set -xeo pipefail
          $(inputs.orphan_readgroup_text ? ["samtools addreplacerg -m orphan_only -o rg_rescue.bam", "--threads", inputs.cpu, "-r", inputs.orphan_readgroup_text.replace(/\\t/g,"\\\\t"), (inputs.reference ? "--reference " + inputs.reference.path : ""), inputs.input_reads.path].join(' ') : "echo 'Orphan ReadGroup Text Not Provided'")
          RG_NUM=`samtools head $(inputs.orphan_readgroup_text ? "rg_rescue.bam" : inputs.input_reads.path) | grep -c ^@RG` || :
          if [ $RG_NUM = 0 ]; then
            echo "INPUT HAS NO READGROUPS! READGROUPS ARE REQUIRED FOR BAM INPUTS!"
            exit 1
          elif [ $RG_NUM -gt 1 ]; then
            samtools split -f '%*_%#.bam' -@ $(inputs.cpu) $(inputs.reference ? '--reference ' + inputs.reference.path : '') $(inputs.orphan_readgroup_text ? "rg_rescue.bam && rm rg_rescue.bam" : inputs.input_reads.path)
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
  orphan_readgroup_text: { type: 'string?', doc: "@RG line text to use for orphan reads."}
  ram: { type: 'int?', default: 16, doc: "GB of RAM to allocate to the task." }
  cpu: { type: 'int?', default: 8, doc: "Minimum reserved number of CPU cores for the task." }
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
