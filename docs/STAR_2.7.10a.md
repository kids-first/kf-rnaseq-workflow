# STAR Aligner v2.7.10a
STAR is a parameter-rich RNAseq aligner.
This tool is based on this release: https://github.com/alexdobin/STAR/tree/2.7.10a.
Depending on the parameters used, run time is typically around 5-10 hours, with spot instance cost of $2.50 - $8.00

## Alignment tool
`tools/star_2.7.10a_align.cwl` <br>
### Required references:
 - `genomeDir`: Reference generated using [Reference generation tool](#reference-generation-tool). Will be tar gzipped, and tool will unpack at run time.

### Required inputs:
 - `readFilesIn1`: Fastq reads file, either uncompressed or gzipped
 - `readFilesIn2`: Required IF data is paired end. Same format as `readFilesIn1`
 - `outSAMattrRGline`: Suggested setting, with TABS SEPARATING THE TAGS, format is: `ID:sample_name LB:aliquot_id PL:platform SM:BSID`.
   For example ID: <br>
 `7316-242  LB:750189   PL:ILLUMINA SM:BS_W72364MN`
 - `outFileNamePrefix`: output files name prefix (including full or relative path).
 Can only be defined on the command line.
 Tool will add '.' after prefix to easily delineate between file name and suffix

### Optional inputs
Try not to panic, but there are many options.
We do set many defaults to make life easier, and should yield reasonable results.
Below is a table of params and values we set by default, as well as suggested changes based on downstream fusion callers, which should have little effect on expression calculations, but potentially valuable increases in sensitivity for fusion calling.

| Param/Description                | WF Default (arriba recommended)  | GTEx/Broad                                        | STAR-Fusion Recommended          |
|----------------------------------|----------------------------------|---------------------------------------------------|----------------------------------|
| runThreadN                       | 16                               | 16                                                | 16                               |
| twopassMode                      | Basic                            | Basic                                             | Basic                            |
| alignSJoverhangMin               | 8                                | 8                                                 | 8                                |
| outFilterMismatchNoverLmax       | 0.1                              | 0.1                                               | 0.1                              |
| outFilterType                    | BySJout                          | BySJout                                           | BySJout                          |
| outFilterScoreMinOverLread       | 0.33                             | 0.33                                              | 0.33                             |
| outFilterMatchNminOverLread      | 0.33                             | 0.33                                              | 0.33                             |
| outReadsUnmapped                 | None                             | None                                              | None                             |
| limitSjdbInsertNsj               | 1200000                          | 1200000                                           | 1200000                          |
| outSAMstrandField                | intronMotif                      | intronMotif                                       | intronMotif                      |
| outFilterIntronMotifs            | None                             | None                                              | None                             |
| alignSoftClipAtReferenceEnds     | Yes                              | Yes                                               | Yes                              |
| quantMode                        | TranscriptomeSAM   GeneCounts    | TranscriptomeSAM   GeneCounts                     | TranscriptomeSAM   GeneCounts    |
| outSAMtype                       | BAM   Unsorted                   | BAM   Unsorted                                    | BAM   Unsorted                   |
| outSAMunmapped                   | Within                           | Within                                            | Within                           |
| genomeLoad                       | NoSharedMemory                   | NoSharedMemory                                    | NoSharedMemory                   |
| chimMainSegmentMultNmax          | 1                                | 1                                                 | 1                                |
| outSAMattributes                 | NH   HI   AS   nM   NM   ch      | NH   HI   AS   nM   NM   ch                       | NH   HI   AS   nM   NM   ch      |
| alignInsertionFlush              | None                             | None                                              | Right                            |
| alignIntronMax                   | 1000000                          | 1000000                                           | 100000                           |
| alignMatesGapMax                 | 1000000                          | 1000000                                           | 100000                           |
| alignSJDBoverhangMin             | 1                                | 1                                                 | 10                               |
| outFilterMismatchNmax            | 999                              | 999                                               | 999                              |
| alignSJstitchMismatchNmax        | 5   -1   5   5                   | 0   -1   0   0                                    | 5   -1   5   5                   |
| alignSplicedMateMapLmin          | 0                                | 0                                                 | 30                               |
| alignSplicedMateMapLminOverLmate | 0.5                              | 0.66                                              | 0                                |
| chimJunctionOverhangMin          | 10                               | 15                                                | 8                                |
| chimMultimapNmax                 | 50                               | 0                                                 | 20                               |
| chimMultimapScoreRange           | 1                                | 1                                                 | 3                                |
| chimNonchimScoreDropMin          | 20                               | 20                                                | 10                               |
| chimOutJunctionFormat            | 1                                | 1                                                 | 1                                |
| chimOutType                      | Junctions   WithinBAM   SoftClip | Junctions   SeparateSAMold   WithinBAM   SoftClip | Junctions   WithinBAM   SoftClip |
| chimScoreDropMax                 | 30                               | 20                                                | 20                               |
| chimScoreJunctionNonGTAG         | -1                               | -1                                                | -4                               |
| chimScoreSeparation              | 1                                | 10                                                | 10                               |
| chimSegmentMin                   | 10                               | 15                                                | 12                               |
| chimSegmentReadGapMax            | 3                                | 0                                                 | 0                                |
| outFilterMultimapNmax            | 50                               | 20                                                | 20                               |
| peOverlapMMp                     | 0.01                             | 0.01                                              | 0                                |
| peOverlapNbasesMin               | 10                               | 0                                                 | 12                               |


### Outputs:
 - `log_progress_out`: Simple progress output. Can use to gauge speed and run time
 - `log_out`: Contains a summary of all params used and reference files
 - `log_final_out`: Overall summary of read mapping statistics
 - `genomic_bam_out`: UNSORTED read mapping to genomic coordinates
 - `junctions_out`: high confidence collapsed splice junctions in tab-delimited form
 - `transcriptome_bam_out`: Read mapping to transcriptome
 - `chimeric_sam_out`: Deprecated output. Incompatible with certain options, has chimeric read alignments
 - `chimeric_junctions`: Chimeric junctions output file. May be used for downstream tools for fusion analysis
 - `gene_counts`: STAR-generated read counts by gene

## Reference generation tool
`tools/star_2.7.10a_genome_generate.cwl` <br>
This tool is used to create the necessary genome indices for STAR.
You should only need to do this once for each gene model, fasta reference, and possibly read length of input files (STAR author says values of 100 should work for most cases)

### Required references:
 - `genomeDir`: Output dirname. Recommend STAR_{version}_GENCODE{version num}
 - `genome_fa`: Fasta file to index. Recommend from GENCODE, PRI assembly. Must unzip first if compressed
 - `gtf`: Matched GTF file to index. Recommend from GENCODE, PRI assembly

### Optional inputs:
 - `runThreadN`: Number of threads to use. Default 16
 - `sjdbOverhang`: Splice junction database overhang. Ideal value is read len minus 1, but default 100 ok for most cases

### Outputs:
 - `star_ref`: A tar gzipped package with all required files and indices to run STAR