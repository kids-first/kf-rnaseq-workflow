cwlVersion: v1.2
class: CommandLineTool
id: alignmentfile_pairedness
doc: |-
  Determines whether a given BAM/SAM is paired or single end.
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
      - |
        /**
        * Given the contents of a file and a string keyword. Return the post-colon value on the line containing the keyWord
        * For instance if the line contents were KeyWord:Value, the function will return Value.
        *
        * @param {String} fileContents - The contents of the file
        * @param {String} keyWord - The exact keyWord for which we will be checking
        * @return {String, null} the string after the keyWord and colon; null if no keyWord found
        */
        function returnKeyValue (fileContents, keyWord) {
          var rows = fileContents.split(/\r?\n/);
          for (var row in rows) {
            if (rows[row].search(keyWord) == 0) {
              return rows[row].split(':')[1];
            }
          }
          return null;
        }
  - class: ResourceRequirement
    coresMin: $(inputs.cpu)
  - class: InitialWorkDirRequirement
    listing:
      - entryname: alignmentfile_pairedness.py
        entry:
          $include: ../scripts/alignmentfile_pairedness.py
  - class: DockerRequirement
    dockerPull: 'quay.io/biocontainers/pysam:0.22.0--py310h41dec4a_0'
baseCommand: [python, alignmentfile_pairedness.py]
stdout: $(inputs.output_filename)
inputs:
  input_reads: {type: File, inputBinding: {position: 2, prefix: "--input_reads"}, doc: "Input BAM/SAM file"}
  input_reference: {type: 'File?', secondaryFiles: [{pattern: '.fai', required: true}], inputBinding: {position: 2, prefix: "--input_reference"}, doc: "For CRAM only, provide the reference file used when making the input_reads"}
  max_reads: {type: 'int?', inputBinding: {position: 2, prefix: "--max_reads"}, doc: "The max number of reads to examine to make PAIRED/SINGLE determination"}
  output_filename: {type: 'string?', default: "pairedness.txt", inputBinding: {position: 9, shellQuote: false, prefix: ">"}, doc: "String to use for output filename"}
  cpu: { type: 'int?', default: 8, inputBinding: {position: 2, prefix: "--threads"}, doc: "CPUs to allocate to this task" }
outputs:
  pairedness_stdout:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
      loadContents: true
  is_paired_end:
    type: boolean
    outputBinding:
      glob: $(inputs.output_filename)
      loadContents: true
      outputEval: |
        $(returnKeyValue(self[0].contents.trim(), "ReadType") == "PAIRED")
