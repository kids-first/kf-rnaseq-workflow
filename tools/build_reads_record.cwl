cwlVersion: v1.2
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: SchemaDefRequirement
  types:
  - $import: ../schema/reads_record_type.yml 
- class: ResourceRequirement
  coresMin: 4
  ramMin: 4000
  https://platform.illumina.com/rdf/ica/resources:tier: economy
  https://platform.illumina.com/rdf/ica/resources:type: standard
  https://platform.illumina.com/rdf/ica/resources:size: medium
baseCommand: [echo, done]
inputs:
  reads1: { type: 'File', secondaryFiles: [{"pattern": "^.bai", required: false }, {"pattern": ".bai", required: false }, {"pattern": "^.crai", required: false }, {"pattern": ".crai", required: false}] } 
  reads2: File?
  outSAMattrRGline: string?
  cram_reference: File?
  is_paired_end: boolean?
  r1_adapter: string?
  r2_adapter: string?
  min_len: int?
  quality_base: int?
  quality_cutoff: int[]?
outputs:
  out_rr:
    type: ../schema/reads_record_type.yml#reads_record
    outputBinding:
      outputEval: |
        $({"reads1": inputs.reads1, "reads2": inputs.reads2, "outSAMattrRGline": inputs.outSAMattrRGline, "cram_reference": inputs.cram_reference, "r1_adapter": inputs.r1_adapter, "r2_adapter": inputs.r2_adapter, "min_len": inputs.min_len, "quality_base": inputs.quality_base, "quality_cutoff": inputs.quality_cutoff, "is_paired_end": inputs.is_paired_end})
