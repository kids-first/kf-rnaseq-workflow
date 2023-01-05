# D3b annoFuse Workflow

## Introduction

In this workflow, annoFuse performs standardization of StarFusion and arriba output files to retain information regarding fused genes, breakpoints, reading frame information as well as annotation from FusionAnnotator, output format description [here](https://github.com/d3b-center/annoFuse/wiki#1-standardize-calls-from-fusion-callers-to-retain-information-regarding-fused-genesbreakpoints-reading-frame-information-as-well-as-annotation-from-fusionannotator). Basic artifact filtering to remove fusions among gene paralogs, conjoined genes and fused genes found in normal samples is also performed by filtering fusions annotated by [FusionAnnotator](https://github.com/d3b-center/FusionAnnotator) with "GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap|ConjoinG". Each fusion call needs at least one junction reads support to be retained as true call. Additionally, if a fusion call has large number of spanning fragment reads compared to junction reads (spanning fragment minus junction read greater than ten), we remove these calls as potential false positives. An expression based filter is also applied, requiring a min FPKM value of 1 for the fusion genes in question.

Please refer to [annoFuse](https://github.com/d3b-center/annoFuse) R package for additional applications like putative oncogene annotations.

## Usage

### Inputs

```yaml
inputs:
  sample_name: { type: 'string', doc: "Sample name used for file base name of all outputs" }
  FusionGenome: { type: 'File', doc: "GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz", "sbg:suggestedValue": { class: 'File', path: '62853e7ad63f7c6d8d7ae5a8', name: 'GRCh38_v39_CTAT_lib_Mar242022.CUSTOM.tar.gz' }}
  genome_untar_path: { type: 'string?', doc: "This is what the path will be when genome_tar is unpackaged", default: "GRCh38_v39_CTAT_lib_Mar242022.CUSTOM" }
  rsem_expr_file: { type: 'File', doc: "gzipped rsem gene expression file" }
  arriba_output_file: { type: 'File', doc: "Output from arriba, usually extension arriba.fusions.tsv" }
  col_num: { type: 'int?', doc: "column number in file of fusion name." }
  star_fusion_output_file: { type: 'File', doc: "Output from arriba, usually extension STAR.fusion_predictions.abridged.coding_effect.tsv" }
  output_basename: { type: 'string', doc: "String to use as basename for outputs" }
```

### Run

1) Outputs from the arriba and STAR Fusion runs are required ahead of time (main RNAseq worflow output)
2) Gzipped rsem counts file, also generated in main RNAseq workflow
3) `FusionGenome` should match what was used to run STAR Fusion

### Outputs

```yaml
outputs:
  annofuse_filtered_fusions_tsv: {type: File, outputSource: annoFuse_filter/filtered_fusions_tsv, doc: "Filtered output of formatted and annotated Star Fusion and arriba results"}
```