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
inputs:
  reads_record:
    type:
      type: record
      fields:
        reads1:
          type: File 
        reads2:
          type: File?
        is_paired_end:
          type: boolean?
        cram_reference:
          type: File?
        outSAMattrRGline:
          type: string?
        r1_adapter:
          type: string?
        r2_adapter:
          type: string?
        min_len:
          type: int?
        quality_base:
          type: int?
        quality_cutoff:
          type: int[]?
  sample_name: {type: 'string?', doc: "Sample ID of the input reads. If not provided, will use reads1 file basename."}
  output_basename: {type: 'string?', doc: "String to use as basename for outputs. Will use read1 file basename if null"}  
  samtools_fastq_cores: {type: 'int?', doc: "Num cores for align2fastq conversion, if input is an alignment file", default: 16}

outputs:
  processed_reads_record:
    type:
      type: record
      fields:
        reads1:
          type: File
        reads2:
          type: File?
        is_paired_end:
          type: boolean
        outSAMattrRGline:
          type: string
    outputSource: construct_out_record/out_record
  cutadapt_stats: {type: 'File?', outputSource: cutadapt_3-4/cutadapt_stats, doc: "Cutadapt stats output, only if adapter is supplied."}
steps:
  basename_picker:
    run: ../tools/basename_picker.cwl
    in:
      root_name:
        source: reads_record
        valueFrom: $(self.reads1.basename.split('.')[0])
      output_basename: output_basename
      sample_name: sample_name
      star_rg_line:
        source: reads_record
        valueFrom: $(self.outSAMattrRGline)
    out: [outname, outsample, outrg]
  alignmentfile_pairedness:
    run: ../tools/alignmentfile_pairedness.cwl
    when: $(inputs.input_reads.reads1.basename.search(/.(b|cr|s)am$/) != -1)
    in:
      input_reads:
        source: reads_record
        valueFrom: $(self.reads1)
      input_reference:
        source: reads_record
        valueFrom: $(self.cram_reference)
    out: [is_paired_end]
  align2fastq:
    # Skip if input is FASTQ already
    run: ../tools/samtools_fastq.cwl
    when: $(inputs.input_reads_1.reads1.basename.search(/.(b|cr|s)am$/) != -1)
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
  cutadapt_3-4:
    # Skip if no adapter given, get fastq from prev step if not null or wf input
    run: ../tools/cutadapter_3.4.cwl
    when: $(inputs.r1_adapter.r1_adapter != null)
    in:
      readFilesIn1:
        source: [align2fastq/fq1, reads_record]
        valueFrom: |
          $(self[0] != null ? self[0] : self[1].reads1)
      readFilesIn2:
        source: [align2fastq/fq2, reads_record]
        valueFrom: |
          $(self[0] != null ? self[0] : self[1].reads2)
      r1_adapter:
        source: reads_record
        valueFrom: $(self.r1_adapter)
      r2_adapter:
        source: reads_record
        valueFrom: $(self.r2_adapter)
      min_len:
        source: reads_record
        valueFrom: $(self.min_len)
      quality_base:
        source: reads_record
        valueFrom: $(self.quality_base)
      quality_cutoff:
        source: reads_record
        valueFrom: $(self.quality_cutoff)
      sample_name: basename_picker/outname
    out: [trimmedReadsR1, trimmedReadsR2, cutadapt_stats]
  construct_out_record:
    run:
      class: CommandLineTool 
      cwlVersion: v1.2
      baseCommand: [echo, done]
      inputs:
        in_r1: File
        in_r2: File?
        in_pe: boolean
        in_rg: string
      outputs:
        out_record:
          type:
            type: record
            fields:
              reads1:
                type: File
              reads2:
                type: File?
              is_paired_end:
                type: boolean
              outSAMattrRGline:
                type: string
          outputBinding:
            outputEval: |
              $({"reads1": inputs.in_r1, "reads2": inputs.in_r2, "is_paired_end": inputs.in_pe, "outSAMattrRGline": inputs.in_rg})
    in:
      in_r1:
        source: [cutadapt_3-4/trimmedReadsR1, align2fastq/fq1, reads_record]
        valueFrom: |
          $(self[0] != null ? self[0] : self[1] != null ? self[1] : self[2].reads1)
      in_r2:
        source: [cutadapt_3-4/trimmedReadsR2, align2fastq/fq2, reads_record]
        valueFrom: |
          $(self[0] != null ? self[0] : self[1] != null ? self[1] : self[2].reads2)
      in_pe:
        source: [reads_record, alignmentfile_pairedness/is_paired_end]
        valueFrom: |
          $(self[0].reads2 != null ? true : self[0].is_paired_end != null ? self[0].is_paired_end : self[1] != null ? self[1] : false)
      in_rg: basename_picker/outrg 
    out: [out_record]

sbg:license: Apache License 2.0
sbg:publisher: KFDRC
$namespaces:
  sbg: https://sevenbridges.com
