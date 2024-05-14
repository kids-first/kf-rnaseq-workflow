cwlVersion: v1.2
class: CommandLineTool
id: rsem-compatibility
doc: "Removes insertions, deletions, and soft-clipped reads for RSEM compatibility"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'staphb/samtools:1.19'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.threads)
    ramMin: 16000
  - class: InitialWorkDirRequirement
    listing:
      - entryname: run_remove_IDS_reads.sh
        entry:
          $include: ../scripts/run_remove_IDS_reads.sh
baseCommand: [bash run_remove_IDS_reads.sh]

inputs:
  input_bam: {type: File, doc: "Input bam file to clean up", inputBinding: { position: 1 }}
  output_basename: {type: string, inputBinding: { position: 2 }}
  threads: {type: 'int?', default: 8, inputBinding: { position: 3 }}


outputs:
  rsem_compatible_bam:
    type: File
    outputBinding:
      glob: '*.Aligned.toTranscriptome_noIDS.out.bam'
