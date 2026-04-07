cwlVersion: v1.2
class: Workflow
id: gencode_benchmark
requirements:
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
inputs:
  baseline_rsems: {type: 'File[]', doc: "*.rsem.genes.results.gz files from baseline" }
  call_rsems: {type: 'File[]', doc: "*.rsem.genes.results.gz files from Set B" }
  baseline_annofuses: {type: 'File[]', doc: "*.annoFuse_filter.tsv files from baseline" }
  call_annofuses: {type: 'File[]', doc: "*.annoFuse_filter.tsv files from Set B" }
  baseline_label: {type: 'string', doc: "Description of baseline files" }
  call_label: {type: 'string', doc: "Description of call files" }
  target_genes: {type: 'File', doc: "TSV containing target ENSGs for heatmap and tsne plots. Recommendation: stable housekeeping gene list" }
  output_basename: {type: 'string', doc: "Basename for output files" }
outputs:
  heatmap: {type: 'File', outputSource: corr_heatmap/plot }
  corr: {type: 'File', outputSource: corr_heatmap/csv }
  tsne: {type: 'File', outputSource: tsne_plot/plot }
  fusion_comp_raw: {type: 'File', outputSource: compare_annofuse_raw/plot }
  fusion_comp_pct: {type: 'File', outputSource: compare_annofuse_pct/plot }
  fusion_comp_uniqs: {type: 'File[]', outputSource:  compare_annofuse_pct/unique_fusions }
steps:
  corr_heatmap:
    run: ../tools/corr_heatmap.cwl
    in:
      gene_list: target_genes
      a_files: baseline_rsems
      b_files: call_rsems
      a_label: baseline_label
      b_label: call_label
      output_basename: output_basename
    out: [ plot, csv ]
  tsne_plot:
    run: ../tools/tsne_plot.cwl
    in:
      gene_list: target_genes
      a_files: baseline_rsems
      b_files: call_rsems
      a_label: baseline_label
      b_label: call_label
      output_basename: output_basename
    out: [ plot ]
  compare_annofuse_raw:
    run: ../tools/compare_annofuse_sets.cwl
    in:
      a_files: baseline_annofuses
      b_files: call_annofuses
      a_label: baseline_label
      b_label: call_label
      output_filename: 
        source: output_basename
        valueFrom: '$(self).annofuse.raw_count.png'
      tool:
        valueFrom: 'annoFuse'
      percent:
        valueFrom: '$(false)'
    out: [ plot, unique_fusions ]
  compare_annofuse_pct:
    run: ../tools/compare_annofuse_sets.cwl
    in:
      a_files: baseline_annofuses
      b_files: call_annofuses
      a_label: baseline_label
      b_label: call_label
      output_filename:
        source: output_basename
        valueFrom: '$(self).annofuse.percent.png'
      tool:
        valueFrom: 'annoFuse'
      percent:
        valueFrom: '$(true)'
    out: [ plot, unique_fusions ]
sbg:license: Apache License 2.0
sbg:publisher: KFDRC
$namespaces:
  sbg: https://sevenbridges.com
