cwlVersion: v1.0
class: Workflow
id: kfdrc-annofuse-wf
requirements:
  - class: MultipleInputFeatureRequirement

inputs:
  sample_name: { type: 'string', doc: "Sample name used for file base name of all outputs" }
  FusionGenome: { type: 'File', doc: "GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz", sbg:suggestedValue: { class: 'File', path: '5d9c8d04e4b0950cce147f94', name: 'GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz' }}
  genome_untar_path: { type: 'string?', doc: "This is what the path will be when genome_tar is unpackaged", default: "GRCh38_v27_CTAT_lib_Feb092018/ctat_genome_lib_build_dir" }
  rsem_expr_file: { type: 'File', doc: "gzipped rsem gene expression file" }
  arriba_output_file: { type: 'File', doc: "Output from arriba, usually extension arriba.fusions.tsv" }
  col_num: { type: 'int?', doc: "column number in file of fusion name." }
  star_fusion_output_file: { type: 'File', doc: "Output from arriba, usually extension STAR.fusion_predictions.abridged.coding_effect.tsv" }
  output_basename: { type: 'string', doc: "String to use as basename for outputs" }

outputs:
  annofuse_filtered_fusions_tsv: { type: 'File?', outputSource: annoFuse_filter/filtered_fusions_tsv, doc: "Filtered output of formatted and annotated Star Fusion and arriba results" }

steps:
  format_arriba_output:
    run: ../tools/format_fusion_file.cwl
    in:
      input_caller_fusion_file: arriba_output_file
      sample_name: sample_name
      caller:
        valueFrom: ${ return "arriba" }
    out:
      [formatted_fusion_tsv]

  format_starfusion_output:
    run: ../tools/format_fusion_file.cwl
    in:
      input_caller_fusion_file: star_fusion_output_file
      sample_name: sample_name
      caller:
        valueFrom: ${ return "starfusion" }
    out:
      [formatted_fusion_tsv]

  annotate_arriba:
    run: ../tools/fusion_annotator.cwl
    in:
      input_fusion_file: format_arriba_output/formatted_fusion_tsv
      genome_tar: FusionGenome
      genome_untar_path: genome_untar_path
      col_num: col_num
      output_basename: output_basename
    out:
      [annotated_tsv]

  annoFuse_filter:
    run: ../tools/annoFuse.cwl
    in:
      arriba_formatted_fusions: annotate_arriba/annotated_tsv
      starfusion_formatted_fusions: format_starfusion_output/formatted_fusion_tsv
      rsem_expr_file: rsem_expr_file
      sample_name: sample_name
      output_basename: output_basename
    out:
      [filtered_fusions_tsv]

$namespaces:
  sbg: https://sevenbridges.com
