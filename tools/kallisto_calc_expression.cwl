cwlVersion: v1.2
class: CommandLineTool
id: kallisto
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/uros_sipetic/kallisto:0.43.1'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 10000
    https://platform.illumina.com/rdf/ica/resources:tier: economy
  - class: SchemaDefRequirement
    types:
    - $import: ../schema/reads_record_type.yml

baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      kallisto quant
      -i $(inputs.transcript_idx.path)
      -o output
      -b 10
      -t 8
      $(inputs.strand ? inputs.strand == "default" ? "" : "--"+inputs.strand : "")
      $(inputs.reads_records.every(function(e) { return !e.is_paired_end }) ? ["--single", "-l", inputs.avg_frag_len, "-s", inputs.std_dev].join(' ') : "")
      $(inputs.reads_records.map(function(e) { return (e.reads2 != null ? [e.reads1.path, e.reads2.path].join(' ') : e.reads1.path) }).join(' '))
      && mv output/abundance.tsv $(inputs.SampleID).kallisto.abundance.tsv
      && gzip $(inputs.SampleID).kallisto.abundance.tsv

inputs:
  transcript_idx: File
  strand: {type: ['null', string], doc: "input none if unstranded, otherwise rf-stranded or fr-stranded"}
  reads_records:
    type:
      type: array
      items: ../schema/reads_record_type.yml#reads_record
  SampleID: string
  std_dev: long?
  avg_frag_len: int?

outputs:
  abundance_out:
    type: File
    outputBinding:
      glob: '*.abundance.tsv.gz'
