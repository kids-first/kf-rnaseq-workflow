cwlVersion: v1.2
class: Workflow
id: preprocess_reads
doc: |
  Preprocess RNAseq reads
  Pick basename
  Check for pairedness
  Convert to FASTQ
  Cutadapt
requirements:
- class: ScatterFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
- class: MultipleInputFeatureRequirement
- class: SchemaDefRequirement
  types:
  - $import: ../schema/reads_record_type.yml
inputs:
  reads_record: {type: ../schema/reads_record_type.yml#reads_record}
  sample_name: {type: 'string?', doc: "Sample ID of the input reads. If not provided, will use reads1 file basename."}
  output_basename: {type: 'string?', doc: "String to use as basename for outputs. Will use read1 file basename if null"}
  samtools_fastq_cores: {type: 'int?', doc: "Num cores for align2fastq conversion, if input is an alignment file", default: 16}
outputs:
  processed_reads_record:
    type: ../schema/reads_record_type.yml#reads_record
    outputSource: build_out_record/out_rr
  cutadapt_stats: {type: 'File?', outputSource: cutadapt_3-4/cutadapt_stats, doc: "Cutadapt stats output, only if adapter is supplied."}
steps:
  alignmentfile_pairedness:
    run: ../tools/alignmentfile_pairedness.cwl
    in:
      input_reads:
        source: reads_record
        valueFrom: $(self.reads1)
      input_reference:
        source: reads_record
        valueFrom: $(self.cram_reference)
    out: [is_paired_end]
  align2fastq:
    run: ../tools/samtools_fastq.cwl
    in:
      input_reads_1:
        source: reads_record
        valueFrom: $(self.reads1)
      SampleID: basename_picker/outname
      cores: samtools_fastq_cores
      is_paired_end:
        source: [reads_record, alignmentfile_pairedness/is_paired_end]
        valueFrom: |
          $(self[0].is_paired_end != null ? self[0].is_paired_end : self[1])
      cram_reference:
        source: reads_record
        valueFrom: $(self.cram_reference)
    out: [fq1, fq2]
  build_out_record:
    run: ../tools/build_reads_record.cwl
    in:
      reads1:
        source: [cutadapt_3-4/trimmedReadsR1, align2fastq/fq1, reads_record]
        valueFrom: |
          $(self[0] != null ? self[0] : self[1] != null ? self[1] : self[2].reads1)
      reads2:
        source: [cutadapt_3-4/trimmedReadsR2, align2fastq/fq2, reads_record]
        valueFrom: |
          $(self[0] != null ? self[0] : self[1] != null ? self[1] : self[2].reads2)
      is_paired_end:
        source: [reads_record, alignmentfile_pairedness/is_paired_end]
        valueFrom: |
          $(self[0].reads2 != null ? true : self[0].is_paired_end != null ? self[0].is_paired_end : self[1] != null ? self[1] : false)
      outSAMattrRGline: basename_picker/outrg
    out: [out_rr]
$namespaces:
  sbg: https://sevenbridges.com
hints:
- class: "sbg:maxNumberOfParallelInstances"
  value: 2
