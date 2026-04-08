# GENCODE Reference Builder

## Introduction

GENCODE is a project that aims to identify and classify all gene features in the human genome.
The tools of the Kids First RNAseq workflow makes extensive use of GENCODE-based references.
These references provide the necessary information to perform read assignment, quanfitifaciton,
and annoatation of RNAseq data.

GENCODE frequently updates their annotations. Over the past five years, GENCODE has released a
new version 2-3 times a year. While version-to-version changes tend to be minor, over the course
of multiple years, the differences can become significant.

For large projects, such as Kids First, it's important to have comparable data generated over large
spans of time. To achieve this we lock our version of GENCODE and generate stable references that
will be used for many years.

The pipeline boils down the creation of the Kids First GENCODE-based references to a series of
version declarations. The pipeline will download the appropriate reference file for each version
provided and perform the necessary processes with the necessary flags to build the KF references.

### Runtime and Instance Recommendation

If time is of the essence, we recommend using direct instances to limit retries of longer tasks.

Expected runtime for tools:
- Under 1 minute: `GTEX_COLLAPSE_ANNOTATION`, `T1K_BUILD`
- Under 10 minutes: `DOWNLOAD_GENCODE`, `KALLISO_INDEX`
- Under 2 hours: `RSEM_GENERATE_GENOME`, `STAR_GENERATE_GENOME`
- 24+ hours: `STAR_FUSION_GENOME_GENERATE`

`STAR_FUSION_GENOME_GENERATE` is an incredibly long running tool. Unfortunately, it seems like
there's little motivation to improve that.

### Software Changes

This pipeline is intended for changes of only the versions of the reference files. Changes in the
software used to build the GENCODE-based references is not possible in CWL. Currently, the pipeline
is using the versions of software used in the KF RNAseq pipeline. See the [dockers](./dockers_rnaseq.md).
If updating the Dockers of this pipeline, make sure the command line parameters have not changed.

## Usage

### Inputs

- `gencode_version`: Version of GENCODE to download from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/
- `ctat_resource_version`: Version of CTAT resource SOURCE file to download from https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/
- `ctat_fusion_version`: Version of CTAT fusion dat.gz to download from https://github.com/FusionAnnotator/CTAT_HumanFusionLib/releases/
- `hla_version`: Version of HLA to download from https://github.com/ANHIG/IMGTHLA/releases/
- `dfam_version`: Version of DFAM to download from https://dfam.org/releases/
- `pfam_version`: Version of PFAM to download from https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/
- `date_today`: Todays date in the following format: Mar012021 (MonthDayYear). Used to name output

### Outputs

- `resource_manifest`: This file contains the URLs of all files downloaded. Mostly for recordkeeping purposes, it also ends up in the `star_fusion_genome` TAR
- `gencode_genome`: GENCODE `*primary_assembly.fa` associated with the `gencode_version`
- `gencode_annotation`: GENCODE `*primary_assembly.annotation.gtf.gz` associated with the `gencode_version`
- `gtex_collapsed_annotation`: Collapsed (i.e., combining all isoforms of a gene into a single transcript) version of `gencode_annotation` created using https://github.com/broadinstitute/gtex-pipeline/tree/master/gene_model with the `--stranded` parameter set
- `rsem_genome`: TAR of transcript references for RSEM created using the `gencode_genome` and `gencode_annotation`
- `kallisto_idx`: IDX of the transcriptome FASTA within `rsem_genome` for Kallisto
- `star_genome`: TAR of binary genome sequence, suffix arrays, text chromosome names/lengths, splice junctions coordinates, and transcripts/genes information made using STAR in `--runMode genomeGenerate` from `gencode_genome` and `gencode_annotation`
- `star_fusion_genome`: STAR fusion genome TAR made from `gencode_genome`, `gencode_annotation`, `ctat_resource_version` SOURCE TAR + filtering rules, `ctat_fusion_version` TAR, `dfam_version` DB, and `pfam_version` DB
- `star_fusion_annot`: Pared version of `star_fusion_genome` containing only `fusion_annot_lib.idx` and `blast_pairs.idx`. Only used to reduce I/O for fusion annotation steps.
- `hla_rna_ref_seqs`: T1k HLA reference sequences made using `gencode_annotation` and `hla_version` release DAT file
- `hla_rna_gene_coords`: T1k HLA reference coordinates made using `gencode_annotation` and `hla_version` release DAT file
