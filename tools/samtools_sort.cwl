cwlVersion: v1.0
class: CommandLineTool
id: samtools_sort
label: "Samtools coordinate sort bam"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 16
    ramMin: 24000

baseCommand: [samtools]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      sort $(inputs.unsorted_bam.path)
      -@ 16
      -m 1G
      -O bam
      > $(inputs.unsorted_bam.nameroot).sorted.bam &&
      samtools
      index
      -@ 16
      $(inputs.unsorted_bam.nameroot).sorted.bam
      $(inputs.unsorted_bam.nameroot).sorted.bai &&
      ${
        var cmd = "echo skip sorting chimeric bam";
        if (inputs.chimeric_sam_out !== null){
          var cmd = "samtools view -bh -@ 16 " + inputs.chimeric_sam_out.path + " -o " + inputs.chimeric_sam_out.nameroot + ".bam";
        }

        return cmd;
      }

inputs:
  unsorted_bam: { type: File, doc: "Bam to sort, likely from STAR" }
  chimeric_sam_out: { type: 'File?', doc: "chimeric bam file - created using v2.6 STAR, probably not in 2.7 STAR" }

outputs:
  sorted_bam:
    type: File
    outputBinding:
      glob: '*.sorted.bam'
  sorted_bai:
    type: File
    outputBinding:
      glob: '*.sorted.bai'
  chimeric_bam_out:
    type: 'File?'
    outputBinding:
      glob: "${var out = ((inputs.chimeric_sam_out === null) ? null : inputs.chimeric_sam_out.nameroot + '.bam'); return out}"
