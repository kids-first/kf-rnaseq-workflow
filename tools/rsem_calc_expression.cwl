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
        var cmd = "rsem-calculate-expression --alignments --append-names --no-bam-output -p 16";
        if (inputs.strandedness != null && inputs.strandedness != "default"){
          cmd += " --strandedness " + inputs.strandedness;
        }
        if (inputs.input_reads_2 != null){
            cmd += " --paired-end";
        }
        cmd += " " + inputs.bam.path + " ./" + inputs.genomeDir.nameroot.replace(".tar", "") + "/" + inputs.genomeDir.nameroot.replace(".tar", "") + " " +  inputs.outFileNamePrefix + ".rsem";
        return cmd
      } &&
      gzip *results

inputs:
  input_reads_2: {type: ['null', File], doc: "Enter paired reads file if applicable, if not leave as null"}
  bam: File
  genomeDir: File
  outFileNamePrefix: string
  strandedness: {type: ['null', string], doc: "Options relative to upstream reads - none, forward, reverse"}

outputs:
  gene_out:
    type: File
    outputBinding: 
      glob: '*genes.results.gz'

  isoform_out:
    type: File
    outputBinding:
      glob: '*isoforms.results.gz'
