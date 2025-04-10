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
- class: SubworkflowFeatureRequirement
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
  basename_picker:
    run: ../tools/basename_picker.cwl
    in:
      root_name:
        source: reads_record
        valueFrom: $(self.reads1.basename.split('.')[0])
      sample_name: sample_name
      star_rg_line:
        source: reads_record
        valueFrom: $(self.outSAMattrRGline)
    out: [outname, outsample, outrg]
  prepare_aligned_reads:
    run: ../subworkflows/prepare_aligned_reads.cwl
    when: $(inputs.reads_record.reads1.basename.search(/.(b|cr|s)am$/) != -1)
    in:
      reads_record: reads_record
      sample_name: basename_picker/outsample
      output_basename: basename_picker/outname
      samtools_fastq_cores: samtools_fastq_cores
    out: [reads1, reads2, is_paired_end, rg_string]
  cutadapt_3-4:
    # Skip if no adapter given, get fastq from prev step if not null or wf input
    run: ../tools/cutadapter_3.4.cwl
    when: $(inputs.reads_record.r1_adapter != null)
    in:
      reads_record: reads_record
      readFilesIn1:
        source: prepare_aligned_reads/reads1
        valueFrom: |
          $(self != null ? self : inputs.reads_record.reads1)
      readFilesIn2:
        source: prepare_aligned_reads/reads2
        valueFrom: |
          $(self != null ? self : inputs.reads_record.reads2)
      r1_adapter:
        valueFrom: $(inputs.reads_record.r1_adapter)
      r2_adapter:
        valueFrom: $(inputs.reads_record.r2_adapter)
      min_len:
        valueFrom: $(inputs.reads_record.min_len)
      quality_base:
        valueFrom: $(inputs.reads_record.quality_base)
      quality_cutoff:
        valueFrom: $(inputs.reads_record.quality_cutoff)
      sample_name:
        source: [output_basename, basename_picker/outname]
        valueFrom: |
          $(self.filter(function(e) { return e != null }).join('.'))
    out: [trimmedReadsR1, trimmedReadsR2, cutadapt_stats]
  build_out_record:
    run: ../tools/build_reads_record.cwl
    in:
      reads1:
        source: [cutadapt_3-4/trimmedReadsR1, prepare_aligned_reads/reads1, reads_record]
        valueFrom: |
          $(self[0] != null ? self[0] : self[1] != null ? self[1] : self[2].reads1)
      reads2:
        source: [cutadapt_3-4/trimmedReadsR2, prepare_aligned_reads/reads2, reads_record]
        valueFrom: |
          $(self[0] != null ? self[0] : self[1] != null ? self[1] : self[2].reads2)
      is_paired_end:
        source: [reads_record, prepare_aligned_reads/is_paired_end]
        valueFrom: |
          $(self[0].reads2 != null ? true : self[0].is_paired_end != null ? self[0].is_paired_end : self[1] != null ? self[1] : false)
      outSAMattrRGline:
        source: [prepare_aligned_reads/rg_string, basename_picker/outrg]
        pickValue: first_non_null
    out: [out_rr]
$namespaces:
  sbg: https://sevenbridges.com
hints:
- class: "sbg:maxNumberOfParallelInstances"
  value: 2
