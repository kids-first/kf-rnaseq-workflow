# Arriba v2.2.1
"[Arriba](https://arriba.readthedocs.io/en/v2.2.1/) is a command-line tool for the detection of gene fusions from RNA-Seq data.
It was developed for the use in a clinical research setting.
Therefore, short runtimes and high sensitivity were important design criteria"
It runs in about 30 minutes - 2 hours costing about $0.17 - $0.86 on a `r5.2xlarge` spot instance.

## Fusion detection
`tools/arriba_fusion_2.2.1.cwl`
### Required references:
 - `reference_fasta`: Fasta reference file used for alignment
 - `gtf_anno`: GTF file used for alignment indexing
### Required inputs:
 - `genome_aligned_bam`: STAR-aligned, **coordinate sorted** bam file
 - `outFileNamePrefix`: String to prepend output file names with
### Optional inputs:
 - `memory`: Amount of ram to dedicate to run the tool
 - `arriba_strand_flag`: auto means auto-detect whether the library is stranded and the type of strandedness.
yes means the library is stranded and the strand of the read designated as first-in-pair matches the transcribed strand.
no means the library is not stranded.
reverse means the library is stranded and the strand of the read designated as first-in-pair is the reverse of the transcribed strand
 - `blacklist`: Path to built-in blacklist to use. Choices are 'blacklist_hg38_GRCh38_v2.2.1.tsv.gz',
  'blacklist_hg19_hs37d5_GRCh37_v2.2.1.tsv.gz', 'blacklist_mm10_GRCm38_v2.2.1.tsv.gz', 'blacklist_mm39_GRCm39_v2.2.1.tsv.gz'.
 - `known`: Path to built-in known/recurrent fusions. Choices are 'known_fusions_hg38_GRCh38_v2.2.1.tsv.gz',
  'known_fusions_hg19_hs37d5_GRCh37_v2.2.1.tsv.gz', 'known_fusions_mm10_GRCm38_v2.2.1.tsv.gz', 'known_fusions_mm39_GRCm39_v2.2.1.tsv.gz'.
 - `tags`: Path to built-in tag files. Typically same as known input
 - `protein_domains`: Path to built-in protein domain annotation. Choices are 'protein_domains_hg38_GRCh38_v2.2.1.gff3',
  'protein_domains_hg19_hs37d5_GRCh37_v2.2.1.gff3', 'protein_domains_mm10_GRCm38_v2.2.1.gff3', 'known_fusions_mm39_GRCm39_v2.2.1.tsv.gz']}].
### Output files:
 - `arriba_fusions`: Fusion results that pass Arriba filters

## Fusion drawing
`tools/arriba_draw_2.2.1.cwl`<br>
This is a nifty tool that creates publication-quality figures of detected fusions.
It works with both Arriba *and* STAR-Fusion outputs
### Required references:
 - `gtf_anno`: GTF file used for alignment indexing
### Required inputs:
 - `genome_aligned_bam`: STAR-aligned, **coordinate sorted** bam file
 - `fusions`: Fusion calls from Arriba OR STAR-Fusion
### Optional inputs:
 - `memory`: Amount of ram to dedicate to run the tool
 - `protein_domains`: Path to built-in protein domain annotation. Choices are 'protein_domains_hg38_GRCh38_v2.2.1.gff3',
  'protein_domains_hg19_hs37d5_GRCh37_v2.2.1.gff3', 'protein_domains_mm10_GRCm38_v2.2.1.gff3', 'known_fusions_mm39_GRCm39_v2.2.1.tsv.gz']}].
 - `cytobands`: Path to built-in coordinates of the Giemsa staining bands to draw ideograms. Choices are 'cytobands_hg38_GRCh38_v2.2.1.tsv',
  'cytobands_hg19_hs37d5_GRCh37_v2.2.1.tsv', 'cytobands_mm10_GRCm38_v2.2.1.tsv', 'cytobands_mm39_GRCm39_v2.2.1.tsv'.
### Outputs:
 - `arriba_pdf`: pdf file with figures, one page per fusion
