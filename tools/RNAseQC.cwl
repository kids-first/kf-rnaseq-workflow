cwlVersion: v1.0
class: CommandLineTool
id: rnaseqc
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'broadinstitute/gtex_rnaseq:V8'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 10000

baseCommand: [bash, -c]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      /src/run_rnaseqc.py
      $(inputs.Aligned_bam.path)
      $(inputs.GTF.path)
      $(input.GenomeRef.path)
      $(inputs.SampleID)
      --output_dir output/

inputs:
  Aligned_bam: File
  GTF: File
  GenomeRef:
    type: File
    secondaryFiles: [.fai, ^.dict]
  SampleID: string

outputs:
  Metrics:
    type: File
    outputBinding:
      glob: '*.metrics.tsv'

  Gene_TPM:
    type: File
    outputBinding:
      glob: '*.gene_tpm.gct'

  Gene_count:
    type: File
    outputBinding:
      glob: '*.gene_reads.gct'

  Exon_count:
    type: File
    outputBinding:
      glob: '*.exon_reads.gct'

