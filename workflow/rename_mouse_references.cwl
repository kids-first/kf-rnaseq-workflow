cwlVersion: v1
class: Workflow
id: rename_mouse_references
label: "Rename STAR reference outputs with _GRCm39"
doc: |
  Rename STAR reference files (e.g., chrLength.txt → chrLength_GRCm39.txt)

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

inputs:
  chrLength: File
  geneInfo: File
  exonInfo: File
  transcriptInfo: File

outputs:
  chrLength_renamed:
    type: File
    outputSource: rename_chrLength/renamed_file
  geneInfo_renamed:
    type: File
    outputSource: rename_geneInfo/renamed_file
  exonInfo_renamed:
    type: File
    outputSource: rename_exonInfo/renamed_file
  transcriptInfo_renamed:
    type: File
    outputSource: rename_transcriptInfo/renamed_file

steps:
  rename_chrLength:
    run: ../tools/rename_file.cwl
    in:
      input_file: chrLength
      new_name: { valueFrom: $(inputs.input_file.basename.replace(/\.txt$/, '_GRCm39.txt')) }
    out: [renamed_file]

  rename_geneInfo:
    run: ../tools/rename_file.cwl
    in:
      input_file: geneInfo
      new_name: { valueFrom: $(inputs.input_file.basename.replace(/\.tab$/, '_GRCm39.tab')) }
    out: [renamed_file]

  rename_exonInfo:
    run: ../tools/rename_file.cwl
    in:
      input_file: exonInfo
      new_name: { valueFrom: $(inputs.input_file.basename.replace(/\.tab$/, '_GRCm39.tab')) }
    out: [renamed_file]

  rename_transcriptInfo:
    run: ../tools/rename_file.cwl
    in:
      input_file: transcriptInfo
      new_name: { valueFrom: $(inputs.input_file.basename.replace(/\.tab$/, '_GRCm39.tab')) }
    out: [renamed_file]

