# RSEM v1.3.1
[RNA-Seq by Expectation-Maximization](https://deweylab.github.io/RSEM/README.html).
This tool is used to estimate gene and isoform levels, outputting normalized read counts, fpkm, tpm.
It runs in about 1.5 - 4 hours, and spot instance cost of $0.60 - $1.50

## Calculate Expression
`tools/rsem_calc_expression.cwl`

### Required references:
 - `genomeDir`: RSEM reference tar ball. This can be created by running [prepare reference](#prepare-reference) tool
### Required inputs:
 - `bam`: Aligned transcriptome bam
 - `outFileNamePrefix`: String to prepend output file names with
 - `strandedness`: "'none' refers to non-strand-specific protocols.
 'forward' means all (upstream) reads are derived from the forward strand.
 'reverse' means all (upstream) reads are derived from the reverse strand"
 - `paired-end`: If input is paired-end, add this flag
### Optional inputs
 - `num_threads`: Num threads to use
 - `append_names`: If available, append gene/tx name to gene/tx id
 - `estimate_rspd`: Set this option if you want to estimate the read start position distribution (RSPD) from data
 - `fragment_length_max`: Maximum read/insert length allowed
### Outputs:
 - `gene_out`: Gene expression estimation
 - `isoform_out`: Transcript isoform level expression estimation

## Prepare reference
`tools/rsem_prepare_reference.cwl` <br>
This tool is used to create the necessary reference for RSEM.
You should only need to do this once for each gene model, fasta reference.

### Required references:
 - `reference_fasta`: Reference fasta file
 - `reference_gtf` **OR** `reference_gtf` : gene model definitions
### Required inputs:
 - `reference_name`: Output file prefix. Recommend format: RSEM_\<SOURCE\>\<Version\>/

### Outputs:
 - `rsem_reference`: RSEM reference tar ball