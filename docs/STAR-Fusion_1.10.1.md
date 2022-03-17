# STAR-Fusion v1.10.1
RNA Fusion detection tool.
This tool was built with STAR aligner v2.7.10a.

## Required references
The main reference required is the `CTAT Genome Lib`.
It can be downloaded from [here](https://github.com/STAR-Fusion/STAR-Fusion/wiki/STAR-Fusion-release-and-CTAT-Genome-Lib-Compatibility-Matrix) or generated using the `tools/star_2.7.10a_genome_generate.cwl` tool.
We have already pre-built this library for GENCODE38 and it is strongly recommended to use a prebuilt one unless absolutely necessary, as the generation process can take 12+ hours.
See [Genome Generate](#genome-generate) section if you really must create one.

### STAR-Fusion
`tools/star_fusion_1.10.1_call.cwl` <br>
This tool uses chimeric junction output from STAR aligner to predict fusions from RNA data.
Estimated run time range: 7-25 minutes.
Estimated spot instance cost: $0.06 - $0.21

#### Inputs:
Along with the tar ball reference mentioned in [Required references](#required-references), the following are required inputs:
 - `Chimeric_junction`: Output junction file from STAR
 - `genome_untar_path`: String describing the path containing `CTAT Genome Lib` when uncompressed
 - `output_basename`: String used as an output file prefix

The following are optional and have default values:
 - `cores`: A resource parameter setting min number of processors available in run time environment. 16 is recommended as going beyond that has diminishing returns.
 - `examine_coding_effect`: Flag to activate examining coding effect. Use null to skip
 - `compress_chimeric_junction`: If part of a workflow, recommend to set this to `true` to take the input junction file and compress it.

#### Outputs:
This tool has just two outputs:
 - `abridged_coding`: tsv file with fusion call results, sans evidence fusion reads, as described [here](https://github.com/STAR-Fusion/STAR-Fusion/wiki#output-from-star-fusion)
 - `chimeric_junction_compressed`: compressed input junction file, if flag was given

### Genome Generate
`tools/star_fusion_1.10.1_gen_reference.cwl`

```yaml
inputs:
  ctat_source: { type: File, doc: "Resource file from https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/", inputBinding: { position: 1} }
  genome_fa: { type: File, doc: "Reference fasta file", inputBinding: { position: 3, prefix: "--genome_fa"} }
  fusion_annot_lib: { type: File, doc: "Can be extracted from ctat source, for from https://github.com/FusionAnnotator/CTAT_HumanFusionLib/releases",
  inputBinding: { position: 3, prefix: "--fusion_annot_lib"} }
  annot_filter_rule: {type: File, doc: "target AnnotFilterRule.pm, from https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB",
  inputBinding: { position: 3, prefix: "--annot_filter_rule"} }
  pfam_db: { type: 'string?', doc: "pfam database. If given will simply pull from current, and run hmmpress", default: 'current'
  inputBinding: { position: 3, prefix: '--pfam_db' } }
  dfam_db: { type: 'File[]',
  doc: "DNA transposable element database (Dfam.hmm), required for repeat masking. Obtain from http://dfam.org/releases/Dfam_3.1/infrastructure/dfamscan/" }
  cores: { type: 'int?', doc: "Num cores to use", default: 16, inputBinding: { position: 3, prefix: '--CPU' } }
  ram: { type: 'int?', doc: "Memory in GB to use", default: 64}
  reference_gtf: { type: File, doc: "gene model definitions." , inputBinding: { position: 2, prefix: '--gtf' } }
  human_gencode_filter: { type: 'string?', doc: "flag for customized prep operations for human/gencode genome and annotation data.",
  default: "--human_gencode_filter" }
  output_dir: { type: string, doc: "Output dirname", inputBinding: { position: 3, prefix: "--output_dir"} }
```

The `doc:` strings describe what is needed and where to obtain resources for generating a `CTAT Genome Lib`.
The output of this tool will be a tar ball of gzipped files, with the dirname being that of your chosen `output_dir`.

