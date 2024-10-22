cwlVersion: v1.2
class: Workflow
id: prepare_aligned_reads
requirements:
- class: ScatterFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
- class: SchemaDefRequirement
  types:
  - $import: ../schema/reads_record_type.yml
inputs:
  reads_record: {type: ../schema/reads_record_type.yml#reads_record}
  sample_name: string
  output_basename: string
  samtools_fastq_cores: int?
outputs:
  reads1:
    type: File
    outputSource: align2fastq/fq1
  reads2:
    type: File?
    outputSource: align2fastq/fq2
  is_paired_end:
    type: boolean
    outputSource: alignmentfile_pairedness/is_paired_end
  rg_string:
    type: string
    outputSource: create_star_rg_line/rg_str
steps:
  samtools_head_rg:
    run: ../tools/samtools_head.cwl
    in:
      input_bam:
        source: reads_record
        valueFrom: $(self.reads1)
      line_filter:
        valueFrom: "^@RG"
    out: [header_file]
  create_star_rg_line:
    run:
      cwlVersion: v1.2
      class: CommandLineTool
      requirements: [{class: InlineJavascriptRequirement}]
      baseCommand: [echo, done]
      inputs:
        rg: {type: 'File', inputBinding: {loadContents: true}}
        sample: {type: string}
      outputs:
        rg_str:
          type: string
          outputBinding:
            outputEval: |-
              ${
                var rgline = inputs.rg.contents.trim().split('\n')[0];
                var fix_sample_rgline = rgline.replace(/\tSM:.+?\t/, "\tSM:" + inputs.sample + "\t");
                var star_rgline = fix_sample_rgline.replace(/^@RG\t/, "");
                return star_rgline;
              }
    in:
      rg: samtools_head_rg/header_file
      sample: sample_name
    out: [rg_str]
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
      SampleID: output_basename
      cores: samtools_fastq_cores
      is_paired_end:
        source: [reads_record, alignmentfile_pairedness/is_paired_end]
        valueFrom: |
          $(self[0].is_paired_end != null ? self[0].is_paired_end : self[1])
      cram_reference:
        source: reads_record
        valueFrom: $(self.cram_reference)
    out: [fq1, fq2]

$namespaces:
  sbg: https://sevenbridges.com
