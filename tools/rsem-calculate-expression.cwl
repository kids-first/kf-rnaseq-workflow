cwlVersion: v1.0
class: CommandLineTool
id: rsem-calculate-expression
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/uros_sipetic/rsem:1.3.1'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 16
    ramMin: 24000

baseCommand: [tar]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -zxf $(inputs.genomeDir.path) &&
      ${
        var cmd = "rsem-calculate-expression --paired-end --alignments --append-names --no-bam-output -p 16";
        if (inputs.forward_prob != null){
          cmd += " --forward-prob " + inputs.forward_prob;
        }
        cmd += " " + inputs.bam.path + " ./" + inputs.genomeDir.nameroot.split('.')[0] + "/" + inputs.genomeDir.nameroot.split('.')[0] + " " +  inputs.outFileNamePrefix + ".rsem";
        return cmd
      }

inputs:
  bam: File
  genomeDir: File
  outFileNamePrefix: string
  forward_prob: {type: ['null', double], doc: "Leave blank if unstranded, 1 if an upstream would be read forward, 0 if read reverse"}

outputs:
  gene_out:
    type: File
    outputBinding: 
      glob: '*genes.results'

  isoform_out:
    type: File
    outputBinding:
      glob: '*isoforms.results'
