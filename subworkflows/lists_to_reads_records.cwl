cwlVersion: v1.2
class: Workflow
id: lists_to_reads_records
doc: |
  Convert three reads list inputs into reads records.

  The three sets of reads lists are:
  - input_aligned_reads
  - input_se_reads and input_se_rg_strs
  - input_pe_reads, input_pe_mates, and input_se_rg_strs

  This workflow will bind those elements together into reads records.
  The reads records contain the following information:
  - reads1
  - reads2
  - outSAMattrRGline
  - cram_reference
  - r1_adapter
  - r2_adapter
  - min_len
  - quality_base
  - quality_cutoff
requirements:
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: SubworkflowFeatureRequirement
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement
- class: SchemaDefRequirement
  types:
  - $import: ../schema/reads_record_type.yml
inputs:
  input_alignment_files: { type: 'File[]?', default: [], secondaryFiles: [{"pattern": "^.bai", required: false }, {"pattern": ".bai", required: false }, {"pattern": "^.crai", required: false }, {"pattern": ".crai", required: false}] }
  input_pe_reads: { type: 'File[]?', default: [] }
  input_pe_mates: { type: 'File[]?', default: [] }
  input_se_reads: { type: 'File[]?', default: [] }
  input_pe_rg_strs: { type: 'string[]?', default: [] }
  input_se_rg_strs: { type: 'string[]?', default: [] }
  is_paired_end: { type: 'boolean?', doc: "Are the alignment files provided paried end?" }
  cram_reference: { type: 'File?', secondaryFiles: [.fai], doc: "If any input alignment files are CRAM, provide the reference used to create them" }
  r1_adapter: { type: 'string?', doc: "!Warning this will be applied to all R1 reads (PE, SE, and reads from alignment files)! If you have multiple adapters, manually trim your reads before input. If they share the same adapter, supply adapter here" }
  r2_adapter: { type: 'string?', doc: "!Warning this will be applied to all R2 reads (PE and reads from alignment files)! If you have multiple adapters, manually trim your reads before input. If they share the same adapter, supply adapter here" }
  min_len: { type: 'int?', doc: "If trimming adapters, what is the minimum length reads should have post trimming" }
  quality_base: { type: 'int?', doc: "Phred scale used for quality scores of the reads" }
  quality_cutoff: { type: 'int[]?', doc: "Quality trim cutoff, see https://cutadapt.readthedocs.io/en/v3.4/guide.html#quality-trimming for how 5' 3' is handled" }
outputs:
  am_reads_records:
    type:
    - 'null'
    - type: array
      items: ../schema/reads_record_type.yml#reads_record
    outputSource: create_reads_records_am/out_rr
    doc: "Reads records made from Aligned Reads input lists"
  pe_fq_reads_records:
    type:
    - 'null'
    - type: array
      items: ../schema/reads_record_type.yml#reads_record
    outputSource: create_reads_records_pe_fq/out_rr
    doc: "Reads records made from Aligned Reads input lists"
  se_fq_reads_records:
    type:
    - 'null'
    - type: array
      items: ../schema/reads_record_type.yml#reads_record
    outputSource: create_reads_records_se_fq/out_rr
    doc: "Reads records made from Aligned Reads input lists"
steps:
  create_reads_records_am:
    run: ../tools/build_reads_record.cwl
    scatter: [reads1]
    in:
      reads1: input_alignment_files
      cram_reference: cram_reference
      is_paired_end: is_paired_end
      r1_adapter: r1_adapter
      r2_adapter: r2_adapter
      min_len: min_len
      quality_base: quality_base
      quality_cutoff: quality_cutoff
    out: [out_rr]
  create_se_reads_null_array:
    run:
      class: CommandLineTool
      cwlVersion: v1.2
      baseCommand: [echo, done]
      requirements: [{class: InlineJavascriptRequirement}]
      inputs:
        in_filelist: { type: { type: array, items: ['null', string] } }
      outputs:
        out_filelist: { type: { type: array, items: ['null', string] }, outputBinding: { outputEval: $(inputs.in_filelist) } }
    in:
      in_filelist:
        source: [input_se_rg_strs, input_se_reads]
        valueFrom: |
          $(self[0].length > 0 ? self[0] : self[1].map(function(e) { return null }))
    out: [out_filelist]
  create_reads_records_se_fq:
    run: ../tools/build_reads_record.cwl
    scatter: [reads1, outSAMattrRGline]
    scatterMethod: dotproduct
    in:
      reads1: input_se_reads
      outSAMattrRGline: create_se_reads_null_array/out_filelist
      r1_adapter: r1_adapter
      min_len: min_len
      quality_base: quality_base
      quality_cutoff: quality_cutoff
    out: [out_rr]
  create_pe_reads_null_array:
    run:
      class: CommandLineTool
      cwlVersion: v1.2
      baseCommand: [echo, done]
      requirements: [{class: InlineJavascriptRequirement}]
      inputs:
        in_filelist: { type: { type: array, items: ['null', string] } }
      outputs:
        out_filelist: { type: { type: array, items: ['null', string] }, outputBinding: { outputEval: $(inputs.in_filelist) } }
    in:
      in_filelist:
        source: [input_pe_rg_strs, input_pe_reads]
        valueFrom: |
          $(self[0].length > 0 ? self[0] : self[1].map(function(e) { return null }))
    out: [out_filelist]
  create_reads_records_pe_fq:
    run: ../tools/build_reads_record.cwl
    scatter: [reads1, reads2, outSAMattrRGline]
    scatterMethod: dotproduct
    in:
      reads1: input_pe_reads
      reads2: input_pe_mates
      outSAMattrRGline: create_pe_reads_null_array/out_filelist
      r1_adapter: r1_adapter
      r2_adapter: r2_adapter
      min_len: min_len
      quality_base: quality_base
      quality_cutoff: quality_cutoff
    out: [out_rr]
$namespaces:
  sbg: https://sevenbridges.com
hints:
- class: "sbg:maxNumberOfParallelInstances"
  value: 2
