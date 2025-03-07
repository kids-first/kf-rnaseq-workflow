# Kids First RNA-Seq Workflow V4

This is the Kids First RNA-Seq pipeline, which calculates gene and transcript isoform expression, detects fusions and splice junctions.
We have transitioned to this current version which upgrades several software components.
Our legacy workflow is still available as [v3.0.1](https://github.com/kids-first/kf-rnaseq-workflow/tree/v3.0.1), and on CAVATICA, [revision 8](https://cavatica.sbgenomics.com/public/apps/cavatica/apps-publisher/kfdrc-rnaseq-workflow/8)

<p align="center">
  <img src="docs/kids_first_logo.svg" alt="Kids First repository logo" width="660px" />
</p>
<p align="center">
  <a href="https://github.com/kids-first/kf-rnaseq-workflow/blob/main/LICENSE"><img src="https://img.shields.io/github/license/kids-first/kf-rnaseq-workflow.svg?style=for-the-badge"></a>
</p>

## Introduction
This pipeline has an optional Cutadapt to trim adapters from the raw reads, alignment-to-FASTQ conversion if necessary, and passes the reads to STAR for alignment.
The alignment output is used by RSEM for gene expression abundance estimation and rMATS for differential alternative splicing events detection.
Additionally, Kallisto is used for quantification, but uses pseudoalignments to estimate the gene abundance from the raw data.
Fusion calling is performed using Arriba and STAR-Fusion detection tools on the STAR alignment outputs.
Filtering and prioritization of fusion calls is done by annoFuse.
Metrics for the workflow are generated by RNA-SeQC.
Junction files for the workflow are generated by rMATS.

If you would like to run this workflow using the CAVATICA public app, a basic primer on running public apps can be found [here](https://www.notion.so/d3b/Starting-From-Scratch-Running-Cavatica-af5ebb78c38a4f3190e32e67b4ce12bb).
Alternatively, if you would like to run it locally using `cwltool`, a basic primer on that can be found [here](https://www.notion.so/d3b/Starting-From-Scratch-Running-CWLtool-b8dbbde2dc7742e4aff290b0a878344d) and combined with app-specific info from the readme below.
This workflow is the current production workflow, equivalent to this [CAVATICA public app](https://cavatica.sbgenomics.com/public/apps#cavatica/apps-publisher/kfdrc-rnaseq-workflow).

### Cutadapt
Cutadapt v3.4: Cut adapter sequences from raw reads if needed.
 - [Github](https://github.com/marcelm/cutadapt)
 - [Publication](https://doi.org/10.14806/ej.17.1.200)
### STAR
STAR v2.7.10a: RNA-Seq raw data alignment.
 - [Github](https://github.com/alexdobin/STAR/tree/2.7.10a)
 - [Publication](https://doi.org/f4h523)
 - [README](docs/STAR_2.7.10a.md)
### [RSEM](docs/RSEM_1.3.1.md)
RSEM v1.3.1: Calculation of gene expression.
 - [Github](https://github.com/deweylab/RSEM/tree/v1.3.1)
 - [Publication](https://doi.org/10.1186/1471-2105-12-323)
 - [README](docs/RSEM_1.3.1.md)
### Kallisto
Kallisto v0.43.1: Raw data pseudoalignment to estimate gene abundance.
 - [Github](https://github.com/pachterlab/kallisto/tree/v0.43.1)
 - [Publication](https://doi.org/10.1038/nbt.3519)
### STAR-Fusion
STAR-Fusion v1.10.1: Fusion detection for `STAR` chimeric reads.
 - [Github](https://github.com/STAR-Fusion/STAR-Fusion/tree/STAR-Fusion-v1.10.1)
 - [Publication](https://doi.org/10.1101/120295)
 - [README](docs/STAR-Fusion_1.10.1.md)
### Arriba
Arriba v2.2.1 Fusion caller that uses `STAR` aligned reads and chimeric reads output.
 - [Github](https://github.com/suhrig/arriba)
 - [Publication](https://doi.org/10.1101/gr.257246.119)
 - [README](docs/ARRIBA_2.2.1.md)
### annoFuse
annoFuse 0.92.0 Filter and prioritize fusion calls.
 - [Github](https://github.com/d3b-center/annoFuse/releases/tag/v0.92.0)
 - [Publication](https://www.biorxiv.org/content/10.1101/839738v3)
 - [README](docs/D3B_ANNOFUSE.md)
### RNA-SeQC
RNA-SeQC v2.3.4 Generate metrics such as gene and transcript counts, sense/antisense mapping, mapping rates, etc.
 - [Github](https://github.com/broadinstitute/rnaseqc)
 - [Publication](https://doi.org/10.1093/bioinformatics/btab135)
### rMATS
rMATS turbo v4.1.2 Computational tool to detect differential alternative splicing events from RNA-Seq data.
 - [Github](https://github.com/Xinglab/rmats-turbo)
 - [Publication](https://doi.org/10.1038/s41596-023-00944-2)
 - [README](docs/D3B_RMATS.md)
### T1k
T1k v1.0.5 Genotype highly polymorphic genes (e.g. HLA) with bulk RNA-seq data.
 - [Github](https://github.com/mourisl/T1K)
 - [Publication](https://doi.org/10.1101/gr.277585.122)
 - [README](docs/T1K_README.md)

## Usage

### Runtime Estimates:
Based on a test set of five input BAMs, CAVATICA compute and storage estimates:
 - Typical 2 hour run time, 10 hours is a higher end possibility
 - Cost:
   - Pure spot instances with no terminations: $2.37 mean
   - Pure on-demand: $5.19 mean
   - Warning: If spot instance kill rate is high, especially for `c5.9xlarge` instance type, the cost could end up greater than on-demand
 - Storage:
   - Total output size 6GB mean
   - Storage estimate ~ $0.14 per month

## Inputs

### Reads Lists

Reads can be provided in both aligned/unaligned SAM/BAM/CRAM and FASTQ/FASTQ.GZ.
These reads are provided through a series of lists:
- `input_alignment_reads`
- `input_pe_reads`
- `input_pe_mates`
- `input_se_reads`

All SAM/BAM/CRAM files should be provided via the `input_alignment_reads` input.
All single end FASTQ files should be provided via the `input_se_reads` input.
For paired end FASTQ files:
- Provide the R1 FASTQ files via the `input_pe_reads` input
- Provide the R2 FASTQ files via the `input_pe_mates` input
- Warning! These lists need to be in the same order! The first file in the `input_pe_reads` list should be the pair of the first file of the `input_pe_mates` list, and so on!

### Read Group Strings and Additional Read Metadata

Users are encouraged to provide read group strings that correspond to their reads
files. For `input_alignment_reads` inputs, this workflow uses the read groups already
present in the SAM/BAM/CRAM header. The only exception to this is that the sample id
`SM` will be updated to whatever `sample_name` the user provides the workflow.

Users can also provide some strandedness and paired information:
- `is_paired_end`: Are the reads paired end?
- `wf_strand_param`: Describe the strandedness of the input data. 'default' for unstranded/auto, 'rf-stranded' if read1 in the FASTQ read pairs is reverse complement to the transcript, 'fr-stranded' if read1 same sense as transcript

If the user does not provide these values, the workflow will attempt to auto-detect
the values. Sometimes the software is incapable of making a determination. If the
user does not provide the value and the workflow is unable to guess the value, the
workflow will fail. Please check the error message to see the cause of the failure
and provide the necessary value.

Users should also provide any adapter information to the workflow. This information includes:
- `r1_adapter`: If the R1 reads still have adapters, supply the adapter sequence here
- `r2_adapter`: If the R2 reads still have adapters, supply the adapter sequence here
- `min_len`: If trimming adapters, what is the minimum length reads should have post trimming
- `quality_base`: Phred scale used for quality scores of the reads
- `quality_cutoff`: Quality trim cutoff, see https://cutadapt.readthedocs.io/en/v3.4/guide.html#quality-trimming for how 5' 3' is handled

At this time the workflow only accepts a single input for these options. If you
have multiple read groups with unique trimming needs, we recommend pre-trimming
the reads before running the workflow.

The workflow is designed to handle multiple read groups; however, not all tools
are capable of handling a mix of single and paired end inputs.

### samtools fastq

When providing one ore more `input_alignment_reads` inputs:
- `samtools_fastq_cores`: Cores for samtools fastq conversion
- `cram_reference`: If input is CRAM, provide the reference FASTA used to build the CRAM

At this time the workflow only accepts a single input for `cram_reference`. If
you have multiple CRAMs with different references, we recommend manually
converting them to BAM before running the workflow.

### STAR Align

STAR Align has many options available to the user at runtime. For the most part
users can leave these fields default. Kids First favors setting/overriding
defaults with "arriba-heavy" specified in [STAR docs](docs/STAR_2.7.10a.md),
however if it is not a tumor sample, then GTEx is preferred. Here are all of
the available options for STAR:
- `STARgenome`: TAR gzipped reference that will be unzipped at run time
- `runThreadN`: Adjust this value to change number of cores used.
- `twopassMode`: Enable two pass mode to detect novel splice events. Default is basic (on).
- `alignSJoverhangMin`: minimum overhang for unannotated junctions. ENCODE default used.
- `outFilterMismatchNoverLmax`: alignment will be output only if its ratio of mismatches to *mapped* length is less than or equal to this value
- `outFilterType`: type of filtering. Normal: standard filtering using only current alignment. BySJout (default): keep only those reads that contain junctions that passed filtering into SJ.out.tab.
- `outFilterScoreMinOverLread`: alignment will be output only if its score is higher than or equal to this value, normalized to read length (sum of mate's lengths for paired-end reads)
- `outFilterMatchNminOverLread`: alignment will be output only if the number of matched bases is higher than or equal to this value., normalized to the read length (sum of mates' lengths for paired-end reads)
- `outReadsUnmapped`: output of unmapped and partially mapped (i.e. mapped only one mate of a paired end read) reads in separate file(s). none (default): no output. Fastx: output in separate FASTA/FASTQ files, Unmapped.out.mate1/2.
- `limitSjdbInsertNsj`: maximum number of junction to be inserted to the genome on the fly at the mapping stage, including those from annotations and those detected in the 1st step of the 2-pass run
- `outSAMstrandField`: Cufflinks-like strand field flag. None: not used. intronMotif (default): strand derived from the intron motif. This option changes the output alignments: reads with inconsistent and/or non-canonical introns are filtered out.
- `outFilterIntronMotifs`: filter alignment using their motifs. None (default): no filtering. RemoveNoncanonical: filter out alignments that contain non-canonical junctions RemoveNoncanonicalUnannotated: filter out alignments that contain non-canonical unannotated junctions when using annotated splice junctions database. The annotated non-canonical junctions will be kept.
- `alignSoftClipAtReferenceEnds`: allow the soft-clipping of the alignments past the end of the chromosomes. Yes (default): allow. No: prohibit, useful for compatibility with Cufflinks
- `quantMode`: types of quantification requested. -: none. TranscriptomeSAM: output SAM/BAM alignments to transcriptome into a separate file GeneCounts: count reads per gene. Choices are additive, so default is 'TranscriptomeSAM GeneCounts'
- `outSAMtype`: type of SAM/BAM output. None: no SAM/BAM output. Otherwise, first word is output type (BAM or SAM), second is sort type (Unsorted or SortedByCoordinate)
- `outSAMunmapped`: output of unmapped reads in the SAM format. None: no output. Within (default): output unmapped reads within the main SAM file (i.e. Aligned.out.sam) Within KeepPairs: record unmapped mate for each alignment, and, in case of unsorted output, keep it adjacent to its mapped mate. Only affects multi-mapping reads
- `genomeLoad`: mode of shared memory usage for the genome file. In this context, the default value makes the most sense, the others are their as a courtesy.
- `chimMainSegmentMultNmax`: maximum number of multi-alignments for the main chimeric segment. =1 will prohibit multimapping main segments
- `outSAMattributes`: a string of desired SAM attributes, in the order desired for the output SAM. Tags can be listed in any combination/order. Please refer to the STAR manual, as there are numerous combinations: https://raw.githubusercontent.com/alexdobin/star_2.7.10a/master/doc/STARmanual.pdf
- `alignInsertionFlush`: how to flush ambiguous insertion positions. None (default): insertions not flushed. Right: insertions flushed to the right. STAR Fusion recommended (SF)
- `alignIntronMax`: maximum intron size. SF recommends 100000
- `alignMatesGapMax`: maximum genomic distance between mates, SF recommends 100000 to avoid readthru fusions within 100k
- `alignSJDBoverhangMin`: minimum overhang for annotated junctions. SF recommends 10
- `outFilterMismatchNmax`: maximum number of mismatches per pair, large number switches off this filter
- `alignSJstitchMismatchNmax`: maximum number of mismatches for stitching of the splice junctions. Value '5 -1 5 5' improves SF chimeric junctions, also recommended by arriba (AR)
- `alignSplicedMateMapLmin`: minimum mapped length for a read mate that is spliced. SF recommends 30
- `alignSplicedMateMapLminOverLmate`: alignSplicedMateMapLmin normalized to mate length. SF recommends 0, AR 0.5
- `chimJunctionOverhangMin`: minimum overhang for a chimeric junction. SF recommends 8, AR 10
- `chimMultimapNmax`: maximum number of chimeric multi-alignments. SF recommends 20, AR 50.
- `chimMultimapScoreRange`: the score range for multi-mapping chimeras below the best chimeric score. Only works with chimMultimapNmax > 1. SF recommends 3
- `chimNonchimScoreDropMin`: int>=0: to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be greater than this value. SF recommends 10
- `chimOutJunctionFormat`: formatting type for the Chimeric.out.junction file, value 1 REQUIRED for SF
- `chimOutType`: type of chimeric output. Args are additive, and defined as such - Junctions: Chimeric.out.junction. SeparateSAMold: output old SAM into separate Chimeric.out.sam file WithinBAM: output into main aligned BAM files (Aligned.*.bam). WithinBAM HardClip: hard-clipping in the CIGAR for supplemental chimeric alignments WithinBAM SoftClip:soft-clipping in the CIGAR for supplemental chimeric alignments
- `chimScoreDropMax`: max drop (difference) of chimeric score (the sum of scores of all chimeric segments) from the read length. AR recommends 30
- `chimScoreJunctionNonGTAG`: penalty for a non-GT/AG chimeric junction. default -1, SF recommends -4, AR -1
- `chimScoreSeparation`: int>=0: minimum difference (separation) between the best chimeric score and the next one. AR recommends 1
- `chimSegmentMin`: minimum length of chimeric segment length, if ==0, no chimeric output. REQUIRED for SF, 12 is their default, AR recommends 10
- `chimSegmentReadGapMax`: maximum gap in the read sequence between chimeric segments. AR recommends 3
- `outFilterMultimapNmax`: max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped. ENCODE value is default. AR recommends 50
- `peOverlapMMp`: maximum proportion of mismatched bases in the overlap area. SF recommends 0.1
- `peOverlapNbasesMin`: minimum number of overlap bases to trigger mates merging and realignment. Specify >0 value to switch on the 'merging of overlapping mates' algorithm. SF recommends 12,  AR recommends 10

These are the defaults set by the workflow:
- `runThreadN`: 36
- `twopassMode`: "Basic"
- `alignSJoverhangMin`: 8
- `outFilterMismatchNoverLmax`: 0.1
- `outFilterType`: "BySJout"
- `outFilterScoreMinOverLread`: 0.33
- `outFilterMatchNminOverLread`: 0.33
- `outReadsUnmapped`: "None"
- `limitSjdbInsertNsj`: 1200000
- `outSAMstrandField`: "intronMotif"
- `outFilterIntronMotifs`: "None"
- `alignSoftClipAtReferenceEnds`: "Yes"
- `quantMode`: TranscriptomeSAM GeneCounts
- `outSAMtype`: "BAM Unsorted"
- `outSAMunmapped`: "Within"
- `genomeLoad`: "NoSharedMemory"
- `chimMainSegmentMultNmax`: 1
- `outSAMattributes`: 'NH HI AS nM NM ch RG'
- `alignInsertionFlush`: "None"
- `alignIntronMax`: 1000000
- `alignMatesGapMax`: 1000000
- `alignSJDBoverhangMin`: 1
- `outFilterMismatchNmax`: 999
- `alignSJstitchMismatchNmax`: "5 -1 5 5"
- `alignSplicedMateMapLmin`: 0
- `alignSplicedMateMapLminOverLmate`: 0.5
- `chimJunctionOverhangMin`: 10
- `chimMultimapNmax`: 50
- `chimMultimapScoreRange`: 1
- `chimNonchimScoreDropMin`: 20
- `chimOutJunctionFormat`: 1
- `chimOutType`: "Junctions WithinBAM SoftClip"
- `chimScoreDropMax`: 30
- `chimScoreJunctionNonGTAG`: -1
- `chimScoreSeparation`: 1
- `chimSegmentMin`: 10
- `chimSegmentReadGapMax`: 3
- `outFilterMultimapNmax`: 50
- `peOverlapMMp`: 0.01
- `peOverlapNbasesMin`: 10

### Arriba
- `arriba_memory`: Mem intensive tool. Set in GB

### STAR Fusion
- `FusionGenome`: STAR-Fusion Cancer Transcriptome Analysis Toolkit (CTAT) Genome lib
- `compress_chimeric_junction`: If part of a workflow, recommend compressing this file as final output

### RNAseQC
- `RNAseQC_GTF`: GTF file from `gtf_anno` that has been collapsed GTEx-style

### kallisto
- `kallisto_idx`: Specialized index of a **transcriptome** FASTA file for kallisto

### RSEM:
- `RSEMgenome`: RSEM reference tar ball
- `estimate_rspd`: Set this option if you want to estimate the read start position distribution (RSPD) from data

### annoFuse:
- `sample_name`: Sample ID of the input reads. If not provided, will use reads1 file basename.
- `annofuse_col_num`: 0-based column number in file of fusion name.
- `fusion_annotator_ref`: Tar ball with fusion_annot_lib.idx and blast_pairs.idx from STAR-Fusion CTAT Genome lib. Can be same as FusionGenome, but only two files needed from that package

### rmats
- `rmats_variable_read_length`: Allow reads with lengths that differ from --readLength to be processed. --readLength will still be used to determine IncFormLen and SkipFormLen.
- `rmats_novel_splice_sites`: Select for novel splice site detection or unannotated splice sites. 'true' to detect or add this parameter, 'false' to disable denovo detection. Tool Default: false
- `rmats_stat_off`: Select to skip statistical analysis, either between two groups or on single sample group. 'true' to add this parameter. Tool default: false
- `rmats_allow_clipping`: Allow alignments with soft or hard clipping to be used.
- `rmats_threads`: Threads to allocate to RMATs.
- `rmats_ram`: GB of RAM to allocate to RMATs.

### T1k
- `run_t1k`: Set to false to disable T1k HLA typing
- `hla_rna_ref_seqs`: FASTA file containing the HLA allele reference sequences for RNA.
- `hla_rna_gene_coords`: FASTA file containing the coordinates of the HLA genes for RNA.

### Run:

1. Reads inputs:
   - For PE FASTQ input, please enter the reads 1 file in `reads1` and the reads 2 file in `reads2`.
   - For SE FASTQ input, enter the single ends reads file in `reads1` and leave `reads2` empty as it is optional.
   - For alignment input (SAM/BAM/CRAM), please enter the reads file in `reads1` and leave `reads2` empty as it is optional.
2. `r1_adapter` and `r2_adapter` are OPTIONAL:
   - If the input reads have already been trimmed, leave these as null and cutadapt step will simple pass on the FASTQ files to STAR.
   - If they do need trimming, supply the adapters and the cutadapt step will trim, and pass trimmed FASTQs along.
   - `min_len` if adapter is trimmed, currently set to min `20` bp. Change this as you see fit
   - `quality_base` set to Phred scale `33` by default if trimming. There was a weird time when `64` was used - change if different
   - `quality_cutoff` if adapter is trimmed and you want to set a min bp quality. A single value will apply to both paired ends, 2 values will allow you to assign a different one to each (unusual)
3. `wf_strand_param` is now *optional* as the workflow will try to determine strandedness for you. Note: if the workflow fails to detect a strandedness, it will fail. If you would like to override autodetect, it is a workflow convenience param so that, if you input the following, the equivalent will propagate to the four tools that use that parameter:
   - `default`: 'rsem_std': null, 'kallisto_std': null, 'rnaseqc_std': null, 'arriba_std': null. This means unstranded or auto in the case of arriba.
   - `rf-stranded`: 'rsem_std': 0, 'kallisto_std': 'rf-stranded', 'rnaseqc_std': 'rf', 'arriba_std': 'reverse'.  This means if read1 in the input FASTQ/BAM is reverse complement to the transcript that it maps to.
   - `fr-stranded`: 'rsem_std': 1, 'kallisto_std': 'fr-stranded', 'rnaseqc_std': 'fr', 'arriba_std': 'yes'. This means if read1 in the input FASTQ/BAM is the same sense (maps 5' to 3') to the transcript that it maps to.
4. Suggested STAR `outSAMattrRGline` format is `ID:sample_name LB:aliquot_id   PL:platform SM:BSID`:
   - For example, `ID:7316-242 LB:750189 PL:ILLUMINA SM:BS_W72364MN`
   - These `KEY:VALUE` fields can be separated by either a whitespace or tab
     character. Any unquoted whitespace will be automatically converted to a tab
     value by STAR. If you wish to include whitespaces in your `VALUE`, you must put
     double quotes around the `VALUE`. For example if you wanted a `DS` key with a
     `I love read groups` value, the entry would look like: `ID:xxx DS:"I love read
     groups"`. See the STAR documentation on `outSAMattrRGline` for complete details.
5. Suggested REFERENCE inputs are:
   - `reference_fasta`: [GRCh38.primary_assembly.genome.fa](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.primary_assembly.genome.fa.gz), will need to unzip
   - `gtf_anno`: [gencode.v39.primary_assembly.annotation.gtf](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.primary_assembly.annotation.gtf.gz), will need to unzip
   - `FusionGenome`: GRCh38_v39_CTAT_lib_Mar242022.CUSTOM.tar.gz. A custom library built using instructions from (https://github.com/STAR-Fusion/STAR-Fusion/wiki/installing-star-fusion#preparing-the-genome-resource-lib), using GENCODE 39 reference.
   - `RNAseQC_GTF`: gencode.v39.primary_assembly.rnaseqc.stranded.gtf OR gencode.v39.primary_assembly.rnaseqc.unstranded.gtf, built using `gtf_anno` and following build instructions [here](https://github.com/broadinstitute/rnaseqc#usage) and [here](https://github.com/broadinstitute/gtex-pipeline/tree/master/gene_model)
   - `RSEMgenome`: RSEM_GENCODE39.tar.gz, built using the `reference_fasta` and `gtf_anno`, following `GENCODE` instructions from [here](https://deweylab.github.io/RSEM/README.html), then creating a tar ball of the results.
   - `STARgenome`: STAR_2.7.10a_GENCODE39.tar.gz, created using the star_2.7.10a_genome_generate.cwl tool, using the `reference_fasta`, `gtf_anno`, and setting `sjdbOverhang` to 100
   - `kallisto_idx`: RSEM_GENCODE39.transcripts.kallisto.idx, built from RSEM GENCODE 39 transcript fasts, in `RSEMgenome` tar ball, following instructions from [here](https://pachterlab.github.io/kallisto/manual)
   - `hla_rna_ref_seqs`: hla_v3.43.0_gencode_v39_rna_seq.fa, created using https://github.com/mourisl/T1K/blob/master/t1k-build.pl with [hla.dat v3.43.0](http://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla.dat) and [GENCODE v39 primary assembly GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.primary_assembly.annotation.gtf.gz)
   - `hla_rna_gene_coords`: hla_v3.43.0_gencode_v39_rna_coord.fa, created using https://github.com/mourisl/T1K/blob/master/t1k-build.pl with [hla.dat v3.43.0](http://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla.dat) and [GENCODE v39 primary assembly GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.primary_assembly.annotation.gtf.gz)
6. rMATS requires the length of the reads in the sample. This workflow will attempt to estimate the read length based on a polling of reads. If the user wishes to override this value they can set `read_length_median` to their desired read length. Additionally, there is a `rmats_variable_read_length` boolean that users can set if their reads are not uniform in length. This workflow will poll the reads and set that value to true if it observes multiple read lengths. Like read length, user-provided input will override this guess.
7. While `output_basename`, `sample_name`, and the `*rg_str` inputs are optional, it is strongly recommended that the user provide these values for data quality purposes. If the user does not provide these values, the basename of the reads1 file will be substituted in their place.
   - `output_basename` and `sample_name` values will become `reads1.basename.split('.')[0]`
   - The STAR align `outSAMattrRGline` value will become:
      - For aligned reads inputs, we will use the RG line set in the BAM header (the `SM` value will be set to what we described above
      - For FASTQ inputs, `ID:reads1.basename.split('.')[0]_1 LB:reads1.basename.split('.')[0] SM:reads1.basename.split('.')[0] PL:Illumina`
   - Additionally for FASTQ inputs, if no `outSAMattrRGline` input is provided a disclaimer will be added to the `@RG` header line that reads: `DS:Values for this read group were auto-generated and may not reflect the true read group information.`

## Outputs
- `cutadapt_stats`: Cutadapt stats output, only if adapter is supplied.
- `STAR_sorted_genomic_cram`: STAR sorted and indexed genomic alignment CRAM
- `STAR_chimeric_junctions`: STAR chimeric junctions
- `STAR_gene_count`: STAR genecounts
- `STAR_junctions_out`: STARjunction reads
- `STAR_final_log`: STAR metricslog file of unique, multi-mapping, unmapped, and chimeric reads
- `STAR-Fusion_results`: STAR fusion detection from chimeric reads
- `arriba_fusion_results`: Fusion output from Arriba
- `arriba_fusion_viz`: pdf output from Arriba
- `RSEM_isoform`: RSEM isoform expression estimates
- `RSEM_gene`: RSEM gene expression estimates
- `RNASeQC_Metrics`: Metrics on mapping, intronic, exonic rates, count information, etc
- `RNASeQC_counts`: Contains gene tpm, gene read, and exon counts
- `kallisto_Abundance`: Gene abundance output from STAR genomic BAM file
- `annofuse_filtered_fusions_tsv`: Filtered fusions called by annoFuse.
- `rmats_filtered_alternative_3_prime_splice_sites_jc`: Alternative 3 prime splice sites JC.txt output from RMATs containing only those calls with 10 or more junction spanning read counts of support
- `rmats_filtered_alternative_5_prime_splice_sites_jc`: Alternative 5 prime splice sites JC.txt output from RMATs containing only those calls with 10 or more junction spanning read counts of support
- `rmats_filtered_mutually_exclusive_exons_jc`: Mutually exclusive exons JC.txt output from RMATs containing only those calls with 10 or more junction spanning read counts of support
- `rmats_filtered_retained_introns_jc`: Retained introns JC.txt output from RMATs containing only those calls with 10 or more junction spanning read counts of support
- `rmats_filtered_skipped_exons_jc`: Skipped exons JC.txt output from RMATs containing only those calls with 10 or more junction spanning read counts of support
- `t1k_genotype_tsv`: Genotyping results from T1k

### Reference build notes
 - STAR-Fusion reference built with command `/usr/local/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl --gtf gencode.v39.primary_assembly.annotation.gtf --annot_filter_rule ../AnnotFilterRule.pm --CPU 36 --fusion_annot_lib ../fusion_lib.Mar2021.dat.gz --genome_fa ../GRCh38.primary_assembly.genome.fa --output_dir GRCh38_v39_CTAT_lib_Mar242022.CUSTOM --human_gencode_filter --pfam_db current --dfam_db human 2> build.errs > build.out &`
 - fusion_annotator_ref built by placing GRCh38_v39_CTAT_lib_Mar242022.CUSTOM/fusion_annot_lib.idx and GRCh38_v39_CTAT_lib_Mar242022.CUSTOM/blast_pairs.idx into its own tar ball
 - kallisto index built using RSEM `RSEM_GENCODE39.transcripts.fa` file as transcriptome FASTA, using command: `kallisto index -i RSEM_GENCODE39.transcripts.kallisto.idx RSEM_GENCODE39.transcripts.fa`
 - RNA-SEQc reference built using [collapse GTF script](https://github.com/broadinstitute/gtex-pipeline/blob/master/gene_model/collapse_annotation.py)
   - Two references needed if data are stranded vs. unstranded
   - Flag `--collapse_only` used for stranded

# [Kids First STAR Diploid Beta](docs/STAR_2.7.11b_DIPLOID.md)
This is an alternative alignment and quantification method currently in beta phase.
It uses DNA variant calls from a patient to create a "personal genome" for improved alignment.
See doc linked in section header.
