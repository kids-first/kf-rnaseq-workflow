# kf-rnaseq-workflow
RNA-Seq workflow for Kids-First DRC

## Introduction
According to GDC's RNA-Seq workflow, this workflow is composed by several steps.

### STAR
RNA-Seq raw data alignment.
### RSEM
Calculation of gene expression.
### STAR-Fusion
Fusion detection for `STAR` chimeric reads.
### HTSeq
htseq-count tool for gene expression count using genomic bam file from `STAR`.
### RNA-SeQC
### Kallisto
Fusion detection using genomic bam file from `STAR`.
### Pizzly
Further fusion analysis tool and generate fusion fasta file from `Kallisto` result.


## Usage

### Inputs common:
```yaml
inputs:
  sample_name: string
  r1_adapter: {type: ['null', string]}
  r2_adapter: {type: ['null', string]}
  STARgenome: File
  RSEMgenome: File
  FusionGenome: File
  runThread: int
  RNAseQC_GTF: File
  GTF_Anno: File
  kallisto_idx: File
  pizzly_transcript_ref: File
```

### Bam input-specific:
```yaml
inputs:
    input_bam: File

```

### Fastq input-specific:
```yaml
inputs:
  reads1: File
  reads2: File
  STAR_outSAMattrRGline: string

```

### Run:

1) For fastq input, run `kfdrc-rnaseq-wf`, for bam input: `kfdrc-rnaseq-wf-bam-in`

2) r1_adapter and r2_adapter are OPTIONAL.  If the input reads have already been trimmed, leave these as null and 
cutadapt step will simple pass on the fastq files to STAR.  If they do need trimming, supply the adapters and the 
cutadapt step will trim, and pass trimmed fastqs along

3) Suggested inputs are:
```text
FusionGenome: GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz
GTF_Anno: gencode.v27.primary_assembly.annotation.gtf
RNAseQC_GTF: gencode.v27.primary_assembly.RNAseQC.gtf
RSEMgenome: RSEM_GENCODE27.tar.gz
STARgenome: STAR_GENCODE27.tar.gz
kallisto_idx: gencode.v27.kallisto.index
pizzly_transcript_ref: gencode.v27.transcripts.pizzly.fa.gz
```

### Outputs:
```yaml
outputs:
  cutadapt_stats: {type: File, outputSource: cutadapt/cutadapt_stats} # only if adapter supplied
  STAR_transcriptome_bam: {type: File, outputSource: star/transcriptome_bam_out}
  STAR_junctions: {type: File, outputSource: star/junctions_out}
  STAR_genomic_bam: {type: File, outputSource: star/genomic_bam_out}
  STAR_gene_counts: {type: File, outputSource: star/gene_counts}
  STAR_chimeric_junctions: {type: File, outputSource: star/chimeric_junctions}
  STAR_chimeric_sam: {type: File, outputSource: star/chimeric_sam_out}
  STAR_Fusion: {type: File, outputSource: star_fusion/fusion_out}
  RSEM_isoform: {type: File, outputSource: rsem/isoform_out}
  RSEM_gene: {type: File, outputSource: rsem/gene_out}
  RNASeQC_Metrics: {type: File, outputSource: rna_seqc/Metrics}
  RNASeQC_Gene_TPM: {type: File, outputSource: rna_seqc/Gene_TPM}
  RNASeQC_Gene_count: {type: File, outputSource: rna_seqc/Gene_count}
  RNASeQC_Exon_count: {type: File, outputSource: rna_seqc/Exon_count}
  kallisto_Abundance: {type: File, outputSource: kallisto/abundance_out}
  kallisto_Fusion: {type: File, outputSource: kallisto/fusion_out}
  Pizzly_Fasta: {type: File, outputSource: pizzly/fusions_fasta}
  Pizzly_unfiltered_Fasta: {type: File, outputSource: pizzly/unfiltered_fusion_fasta}
  ```