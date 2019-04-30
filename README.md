# kf-rnaseq-workflow
RNA-Seq workflow for Kids-First DRC

## Introduction
Follow mostly the workflow of the GDC.

### STAR
RNA-Seq raw data alignment.
### RSEM
Calculation of gene expression.
### STAR-Fusion
Fusion detection for `STAR` chimeric reads.
### RNA-SeQC
Generate metrics such as gene and transcript counts, sense/antisene mapping, mapping rates, etc
### Kallisto
Fusion detection using genomic bam file from `STAR`.
### Pizzly
Further fusion analysis tool and generate fusion txt (in tsv format) file from `Kallisto` result.
### arriba
Another fusion caller that uses star output.


## Usage

### Inputs common:
```yaml
inputs:
  sample_name: string
  r1_adapter: {type: ['null', string]}
  r2_adapter: {type: ['null', string]}
  STAR_outSAMattrRGline: string
  STARgenome: File
  RSEMgenome: File
  reference_fasta: File
  gtf_anno: File
  wf_strand_param: {type: ['null', string], doc: "use 'default' or leave blank for unstranded/auto, rf_stranded if read1 in the fastq read pairs is reverse complement to the transcript, fr-stranded if read1 same sense as transcript"}
  FusionGenome: File
  runThread: int
  RNAseQC_GTF: File
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

```

### Run:

1) For fastq input, run `kfdrc-rnaseq-wf`, for bam input: `kfdrc-rnaseq-wf-bam-in`

2) r1_adapter and r2_adapter are OPTIONAL.  If the input reads have already been trimmed, leave these as null and 
cutadapt step will simple pass on the fastq files to STAR.  If they do need trimming, supply the adapters and the 
cutadapt step will trim, and pass trimmed fastqs along

3) `wf_strand_param` is a workflow convenience param so that, if you input the following, the equivalent will propagate 
to the four tools that use that parameter:
    - `default`: 'rsem_std': null, 'kallisto_std': null, 'rnaseqc_std': null, 'arriba_std': null. This means unstranded or auto in the case of arriba.
    - `rf_stranded`: 'rsem_std': 0, 'kallisto_std': 'rf-stranded', 'rnaseqc_std': 'rf', 'arriba_std': 'reverse'.  This means if read1 in the input fastq/bam is reverse complement to the transcript that it maps to.
    - `fr-stranded`: 'rsem_std': 1, 'kallisto_std': 'fr-stranded', 'rnaseqc_std': 'fr', 'arriba_std': 'yes'. This means if read1 in the input fastq/bam is the same sense (maps 5' to 3') to the transcript that it maps to.

4) Suggested `STAR_outSAMattrRGline`, with **TABS SEPARATING THE TAGS**,  format is:
    `ID:sample_name LB:aliquot_id   PL:platform SM:BSID`
    for example `ID:7316-242   LB:750189 PL:ILLUMINA SM:BS_W72364MN`
5) Suggested inputs are:
```text
FusionGenome: GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz
gtf_anno: gencode.v27.primary_assembly.annotation.gtf
RNAseQC_GTF: gencode.v27.primary_assembly.RNAseQC.gtf
RSEMgenome: RSEM_GENCODE27.tar.gz
STARgenome: STAR_GENCODE27.tar.gz
reference_fasta: GRCh38.primary_assembly.genome.fa
kallisto_idx: gencode.v27.kallisto.index
pizzly_transcript_ref: gencode.v27.transcripts.pizzly.fa.gz
```

### Outputs:
```yaml
outputs:
  cutadapt_stats: {type: File, outputSource: cutadapt/cutadapt_stats} # only if adapter supplied
  STAR_transcriptome_bam: {type: File, outputSource: star/transcriptome_bam_out}
  STAR_sorted_genomic_bam: {type: File, outputSource: samtools_sort/sorted_bam}
  STAR_sorted_genomic_bai: {type: File, outputSource: samtools_sort/sorted_bai}
  STAR_chimeric_bam_out: {type: File, outputSource: samtools_sort/chimeric_bam_out}
  STAR_chimeric_junctions: {type: File, outputSource: star_fusion/chimeric_junction_compressed}
  STAR_gene_count: {type: File, outputSource: star/gene_counts}
  STAR_junctions_out: {type: File, outputSource: star/junctions_out}
  STAR_final_log: {type: File, outputSource: star/log_final_out}
  STAR-Fusion_results: {type: File, outputSource: star_fusion/abridged_coding}
  pizzly_fusion_results: {type: File, outputSource: pizzly/fusions_flattened}
  arriba_fusion_results: {type: File, outputSource: arriba_fusion/arriba_fusions}
  arriba_fusion_viz: {type: File, outputSource: arriba_fusion/arriba_pdf}
  RSEM_isoform: {type: File, outputSource: rsem/isoform_out}
  RSEM_gene: {type: File, outputSource: rsem/gene_out}
  RNASeQC_Metrics: {type: File, outputSource: rna_seqc/Metrics}
  RNASeQC_counts: {type: File, outputSource: supplemental/RNASeQC_counts} # contains gene tpm, gene read, and exon counts
  kallisto_Abundance: {type: File, outputSource: kallisto/abundance_out}
  ```
