cwlVersion: v1.2
class: CommandLineTool
id: samtools_readlength_bam
doc: |-
  samtools view file.bam | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c | sort -nr -k1,1
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9'
baseCommand: [samtools,view]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.input_bam.path) | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c > $(inputs.input_bam.nameroot).bam_readlength | sort -nr -k1,1

inputs:
  input_bam: { type: File, doc: "Input bam file"}

outputs:
  output:
    type: File
    outputBinding:
      glob: "*.bam_readlength"
  top_readlength:
    type: int
    outputBinding:
      glob: "*.bam_readlength"
      loadContents: true
      outputEval: |
        ${
          var rows = self[0].contents.split(/\r?\n/).slice(0,-1);
          return rows[0].split(/\s/).pop();
        }
  variable_readlength:
    type: boolean
    outputBinding:
      glob: "*.bam_readlength"
      loadContents: true
      outputEval: |
        ${
          var rows = self[0].contents.split(/\r?\n/).slice(0,-1);
          return rows.length > 1;
        }
