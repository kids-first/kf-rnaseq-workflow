cwlVersion: v1.2
class: Workflow
id: gencode_benchmark
label: GENCODE Bencharmking
doc: |
  # GENCODE Benchmark

  ## Introduction

  The GENCODE benchmark workflow generates a series of plots and other files used to compare two sets of data. While made specifically for GENCODE comparisons,
  the workflow will compare any two sets of KF RNAseq data. In total the workflow generates four plots; two are generated from the annoFuse inputs and two are
  generated from the RSEM gene results of the Kids First RNAseq pipeline. The first of these plots is a raw count comparison of the annoFuse fusions between 
  the two sets. For each sample, it shows the number of fusions that are shared or unique to each sample. The second of these plots is the first plot but
  rather than showing counts it shows percentage of total fusions. The third plot is a heatmap showing the correlation of RSEM TMP values among samples. The 
  fourth and final plot is an RSEM TPM tSNE plot for all the samples. In additon to the plot the workflow also returns the a CSV containing the raw correlation
  values used to generate the heatmap as well as a series of files that detail the fusions that are unique to each Sample and Set.

  These plots and other files are for evaluating the effect of going from one GENCODE version to another. The general idea is that greater differences between
  versions will result in:
  - more pronounced version-based tSNE clustering
  - lower same-sample heatmap correlation
  - greater fusion count differences

  ## Usage

  ### Inputs

  - `baseline_annofuses`: *.annoFuse_filter.tsv files from baseline samples
  - `baseline_rsems`: *.rsem.genes.results.gz files from baseline samples
  - `baseline_label`: Description of baseline samples
  - `call_annofuses`: *.annoFuse_filter.tsv files from call samples
  - `call_rsems`: *.rsem.genes.results.gz files from call samples
  - `call_label`: Description of samples
  - `target_genes`: TSV containing target ENSGs for heatmap and tsne plots. Recommendation: stable housekeeping gene list
  - `output_basename`: Basename for output files

  ### Outputs

  - `heatmap`: Heatmap plot comparing baseline and call samples by RSEM gene TPM values for ENSGs in the target_genes file
  - `corr`: Heatmap CSV comparing baseline and call samples by RSEM gene TPM values for ENSGs in the target_genes file
  - `tsne`: tSNE plot clustering baseline and call samples by RSEM gene TMP values for ENSGs in the target_genes files
  - `fusion_comp_raw`: Stacked bar chart showing shared and unique fusions for each sample
  - `fusion_comp_pct`: Stacked bar chart showing percente of shared and unique fusions for each sample
  - `fusion_comp_uniqs`: Files with unique fusions for each sample x label
requirements:
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
inputs:
  baseline_rsems: {type: 'File[]', doc: "*.rsem.genes.results.gz files from baseline"}
  call_rsems: {type: 'File[]', doc: "*.rsem.genes.results.gz files from Set B"}
  baseline_annofuses: {type: 'File[]', doc: "*.annoFuse_filter.tsv files from baseline"}
  call_annofuses: {type: 'File[]', doc: "*.annoFuse_filter.tsv files from Set B"}
  baseline_label: {type: 'string', doc: "Description of baseline files"}
  call_label: {type: 'string', doc: "Description of call files"}
  target_genes: {type: 'File', doc: "TSV containing target ENSGs for heatmap and tsne plots. Recommendation: stable housekeeping gene
      list"}
  output_basename: {type: 'string', doc: "Basename for output files"}
outputs:
  heatmap: {type: 'File', outputSource: corr_heatmap/plot}
  corr: {type: 'File', outputSource: corr_heatmap/csv}
  tsne: {type: 'File', outputSource: tsne_plot/plot}
  fusion_comp_raw: {type: 'File', outputSource: compare_annofuse_raw/plot}
  fusion_comp_pct: {type: 'File', outputSource: compare_annofuse_pct/plot}
  fusion_comp_uniqs: {type: 'File[]', outputSource: compare_annofuse_pct/unique_fusions}
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
    out: [plot, csv]
  tsne_plot:
    run: ../tools/tsne_plot.cwl
    in:
      gene_list: target_genes
      a_files: baseline_rsems
      b_files: call_rsems
      a_label: baseline_label
      b_label: call_label
      output_basename: output_basename
    out: [plot]
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
    out: [plot, unique_fusions]
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
    out: [plot, unique_fusions]
sbg:license: Apache License 2.0
sbg:publisher: KFDRC
$namespaces:
  sbg: https://sevenbridges.com
sbg:categories:
- ANNOFUSE
- CORR
- FUSION
- GENCODE
- GENE
- HEATMAP
- RNASEQ
- RSEM
- TSNE
