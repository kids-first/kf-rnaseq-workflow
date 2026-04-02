cwlVersion: v1.2
class: CommandLineTool
id: generate_fusion_annot_ref
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.cpu)
  - class: DockerRequirement
    dockerPull: ubuntu:24.04
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      tar xf $(inputs.fusion_tar.path) $(inputs.fusion_tar.basename.replace(/.tar.gz$/,""))/fusion_annot_lib.idx $(inputs.fusion_tar.basename.replace(/.tar.gz$/,""))/blast_pairs.idx
      && tar czf $(inputs.output_basename).tar.gz $(inputs.fusion_tar.basename.replace(/.tar.gz$/,""))/fusion_annot_lib.idx $(inputs.fusion_tar.basename.replace(/.tar.gz$/,""))/blast_pairs.idx

inputs:
  fusion_tar: { type: 'File', doc: "Fusion tar from fusion generate genome" }
  output_basename: { type: 'string?', default: "", doc: "Output name for pared tar" } 
  cpu: { type: 'int?', doc: "Num processing threads to use", default: 1 }
  ram: { type: 'int?', doc: "Num GB memory to make available", default: 2 }
outputs:
  fusion_annot_ref: { type: 'File', outputBinding: { glob: '$(inputs.output_basename).tar.gz' } }
