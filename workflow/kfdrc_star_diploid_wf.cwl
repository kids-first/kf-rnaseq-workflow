cwlVersion: v1.2
class: Workflow
id: kfdrc-star-diploid-wf
label: KFDRC STAR Diploid Workflow
doc: |-
  # Kids First Kids First STAR Diploid Beta

  This is an alternative alignment and quantification method currently in **_beta phase_**.
  It uses a patient's DNA variant calls to create a "personal genome" (PG) for improved alignment.
  It is purported to have fewer multi-mapping/better unique mapping for potentially improved gene and isoform level quantification.
  It cannot be used in fusion calling
  The STAR Diploid mode has a known bug with an unknown cause manifests as a seg fault for up to 20% of normal sample inputs and 100% of tumor inputs. 

  ![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

  ## Introduction
  This pipeline runs the following steps:
  1. STAR Genome Generate per individual, or will skip if one has been generated already.
      -  Strip existing annotations as defined by the user. Since STAR requires an uncompressed vcf, this helps make the input file size smaller as only variant calls and `FORMAT` info are needed
      -  Filtering steps for input patient DNA variant calls in order to focus on higher quality calls.
      Recommended filtering criteria recommendations are still being established, but various scenarios and suggestions from our and partner institutions can be found in the inputs section.
      - Use a genome fasta and gtf - recommend fasta matches input DNA calls - in conjunction with filtered DNA VCF to create PG
  1. If needed, convert input bam reads to fastq
  1. If needed, cutadapt to remove any adapters
  1. Run STAR aligner 
      - Use PG as refs
      - Align input reads
  1. Run custom tool to filter bam for RSEM. Removes indels and soft-clipped reads
  1. RSEM quantification

  ### Cutadapt
  [Cutadapt v3.4](https://github.com/marcelm/cutadapt) Cut adapter sequences from raw reads if needed.
  ### STAR v2.7.11b_alpha_2024-03-29
  [STAR v2.7.11b_alpha_2024-03-29](https://doi.org/f4h523) RNA-Seq raw data alignment.
  ### [RSEM](docs/RSEM_1.3.1.md)

  ## INPUTS
  A brief note - many references and filtering requirements take a bit of up front work. In the [Filtering and Input Appendix](#filtering-and-input-appendix) section, we try to lay this out. The workflow has more inputs that this README will not cover, but can be tweaked by an advanced user.
  ### Common Required
  These are inputs used in multiple steps
   - `output_basename`: String to prepend to results files from STAR and RSEM
  ### STAR Genome Generate
  If a pre-existing PG does not exist, need the following inputs to create:
   - `input_vcf`: Currently not standardized, matched DNA variant calls. For INCLUDE dataset, trio calls were used when available, otherwise singe sample genotyping - both from GATK workflows
   - `strip_info`: Given that input vcf needs to be uncompressed, stripping `INFO` is a good way to reduce file size. Current recommended strip based on typical KF runs: `INFO/CLNDISDB,INFO/CLNDISDBINCL,INFO/CLNDN,INFO/CLNDNINCL,INFO/CLNHGVS,INFO/CLNREVSTAT,INFO/CLNSIG,INFO/CLNSIGCONF,INFO/CLNSIGINCL,INFO/CLNVC,INFO/CLNVCSO,INFO/CLNVI,INFO/CSQ,INFO/ClippingRankSum,INFO/DB,INFO/DP,INFO/DS,INFO/END,INFO/ExcessHet,INFO/FS,INFO/HaplotypeScore,INFO/InbreedingCoeff,INFO/Intervar,INFO/Intervar_STATUS,INFO/MLEAC,INFO/MLEAF,INFO/MQ,INFO/MQRankSum,INFO/NEGATIVE_TRAIN_SITE,INFO/OLD_VARIANT,INFO/POSITIVE_TRAIN_SITE,INFO/QD,INFO/RAW_MQ,INFO/ReadPosRankSum,INFO/SOR,INFO/VQSLOD,INFO/culprit,INFO/gnomad_3_1_1_AC,INFO/gnomad_3_1_1_AC_controls_and_biobanks,INFO/gnomad_3_1_1_AC_popmax,INFO/gnomad_3_1_1_AF,INFO/gnomad_3_1_1_AF_controls_and_biobanks,INFO/gnomad_3_1_1_AF_non_cancer,INFO/gnomad_3_1_1_AF_popmax,INFO/gnomad_3_1_1_AN,INFO/gnomad_3_1_1_AN_controls_and_biobanks,INFO/gnomad_3_1_1_AN_popmax,INFO/gnomad_3_1_1_nhomalt,INFO/gnomad_3_1_1_nhomalt_popmax,INFO/gnomad_3_1_1_primate_ai_score,INFO/gnomad_3_1_1_splice_ai_consequence`
    - `include_expression`: Filters DNA vcf for high quality variants. Current recommended:
      - Trio called VCF: `STRLEN(REF)<=50 && STRLEN(ALT)<=50 && FILTER="PASS" && GT="alt"`
      - Single sample VCF:  `STRLEN(REF)<=50 && STRLEN(ALT)<=50 && FILTER="PASS"`
   - `subtract_bed`: Recommend to filter regions from repeat and low complexity regions. Recommend obtaining repeat-masker bed file from UCSC, run bedtools sort + merge to simplify. Removes variant calls from `input_vcf` from notoriously difficult regions
   - `vcf_sample_name`: **If input is trio**, provide the patient sample name to ensure desired `include_expression` is applied to the specific patient
   - `genome_fa`: Should match input used for DNA. For KF/INCLUDE, recommend `Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta`.
   - `genomeTransformType`: `Diploid`, set by default
   - `gtf`: Recommend `PRI` assembly from [GENCODE version 45](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.primary_assembly.annotation.gtf.gz) for CFDE
   - `sjdbOverhang`: Default is 100. Normally fine as-is, but for PG should probably just set it to read length minus 1
   ### STAR aligner
    - `reads1`: BAM/CRAM/FASTQ reads input. If BAM/CRAM workflow will convert to FASTQ. If FASTQ, read 1 file would go here
    - `reads2`: If FASTQ and paired end, mates file, aka read 2 goes here
    - `cram_reference`: If input reads are CRAM, provide alignment FASTA reference used for it here
    - `outSAMattrRGline`: Output alignment read group. With **_tabs separating the tags_**, format is: ID:sample_name LB:aliquot_id PL:platform SM:BSID for example ID:7316-242 LB:750189 PL:ILLUMINA SM:BS_W72364MN
    - `genomeDir`: If pre-built tar-gzipped PG exists, provide here to skip STAR Genome step
  ### RSEM
   - `wf_strand_param`: Strandedness of input reads. Default is `rf-stranded`. Use 'default' for unstranded/auto, 'rf-stranded' if read1 in the fastq read pairs is reverse complement to the transcript, 'fr-stranded' if read1 same sense as transcript
   - `RSEMgenome`: RSEM reference tar ball.
  ## OUTPUTS
   - `star_ref`: If existing `genomeDir` tar ball was not provided, workflow will have created and provided the patient's PG. This can be re-used in the event a user wants to align another sample frm the patient and/or try different aligner parameters
   - `debug_log`: Log output from STAR GEnome Generate 
   - `STAR_sorted_genomic_cram`: Aligned reads to genome in CRAM format
   - `STAR_transcriptome_bam`: Typically not kept as it's seldom re-used, given that this is in beta phase, we'll keep this
   - `STAR_gene_count`: Gene counts from STAR
   - `STAR_junctions_out`: STAR junctions file
   - `STAR_final_log`: STAR metrics log file of unique, multi-mapping, unmapped, and chimeric reads
   - `RSEM_isoform` RSEM isoform expression estimates
   - `RSEM_gene`: RSEM gene expression estimates

  ## Filtering and Input Appendix
  ### Input creation
   - `genome_fa` was generated by first downloading the [complete hg38 fasta](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta) from Broad, then creating a chr1-22,X,Y,M chromosome list, and running the following commands to get the fasta file and index: 
     ```sh
     samtools faidx Homo_sapiens_assembly38.fasta -r chr_list.txt > Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
     samtools faidx Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
     ```
   - `subtract_bed` was generated by navigating to https://genome.ucsc.edu/cgi-bin/hgTables, then set the following:
     - **group**: `Repeats`
     - **track**: `RepeatMasker`
     - **output format**: `BED`
     - Result was then piped to `bedtools sort | bedtools merge | gzip > rpt_merge_sort.bed.gz` 
  ### VCF filtering
  #### Trio VCFS
  Two main sources of soft `FILTER` are added during creation using [this workflow](https://github.com/kids-first/kf-jointgenotyping-workflow/blob/v2.4.0/README.md):
   - Adding tranches during VQSR
   - Adding `lowGQ` filter `GQ < 20.0`
  These variants get removed during bcftools step of filtering on `PASS`
  ### Single sample VCFS
  These have even more `FILTER` added that are dropped on `PASS`:
   - Common
     - ExcessHet: `ExcessHet > 54.69`
     - QD2: `QD < 2.0`
     - QUAL30: `QUAL < 30.0`
   - SNP only:
     - SOR3_SNP: `SOR > 3.0`
     - FS60_SNP: `FS > 60.0`
     - MQ40_SNP: `MQ < 40.0`
     - MQRankSum-12.5_SNP: `MQRankSum < -12.5`
     - ReadPosRankSum-8_SNP: `ReadPosRankSum < -8.0`
   - INDEL only:
     - QD2: `QD < 2.0`
     - QUAL30: `QUAL < 30.0`
     - FS200_INDEL: `FS > 200.0`
     - ReadPosRankSum-20_INDEL: `ReadPosRankSum < -20.0`
  ### Broad joint cohort calls
  Provided by our collaborators at Broad, these are recommended for joint cohort calls:
   - Genotype-quality (GQ): filter out genotypes with a quality below this threshold. Default: 20 - 
   - Allelic-balance: filter out genotypes with allelic balance outside of [1-{allelic-balance}, {allelic-balance}]. Default: 0.8 - 
   - Missingness-threshold: filter out variants with missingness above this threshold. Default: 0.02 - 
   - HWE N-individuals: minimum number of individuals within a subpopulation to perform HWE test. Default: 100 - 
   - HWE: the autosomal (or non-pseudoautosomal X in females only) site failed Hardy-Weinberg equilibrium expectations in (sub)population(s). Minimum acceptable p-value for HWE: Default: 1e-8 - 
   - hmiss pval-threshold: minimum acceptable p-value for hmiss: Default: 1e-8 - 
   - Sample-blacklist: list of samples to exclude. - 
   - ExcessHet (excess of heterozygosity):filter out variants > 54.69 - 
   - LCR (low-complexity regions): filter out variants in low-complexity regions - 
   - Monomorphic alleles: filter out monomorphic alleles (0/0) - 
   - VQSR: Variant Quality Score Recalibration. It is a sophisticated filtering technique applied on the variant callset that uses machine learning to model the technical profile of variants in a training set and uses that to filter out probable artifacts from the callset. We remove variants with VQSR filter. - 
   - Use LCR-hg38-noHLA.interval_list to filter out calls in those regions


requirements:
- class: MultipleInputFeatureRequirement
- class: SubworkflowFeatureRequirement
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
  # Strip, subset and PASS vars
  input_vcf: {type: 'File?', secondaryFiles: ['.tbi']}
  strip_info: {type: 'string?', doc: "If given, remove previous annotation information based on INFO file, i.e. to strip VEP info,
      use INFO/ANN", default: "INFO/CLNDISDB,INFO/CLNDISDBINCL,INFO/CLNDN,INFO/CLNDNINCL,INFO/CLNHGVS,INFO/CLNREVSTAT,INFO/CLNSIG,INFO/CLNSIGCONF,INFO/CLNSIGINCL,INFO/CLNVC,INFO/CLNVCSO,INFO/CLNVI,INFO/CSQ,INFO/ClippingRankSum,INFO/DB,INFO/DP,INFO/DS,INFO/END,INFO/ExcessHet,INFO/FS,INFO/HaplotypeScore,INFO/InbreedingCoeff,INFO/Intervar,INFO/Intervar_STATUS,INFO/MLEAC,INFO/MLEAF,INFO/MQ,INFO/MQRankSum,INFO/NEGATIVE_TRAIN_SITE,INFO/OLD_VARIANT,INFO/POSITIVE_TRAIN_SITE,INFO/QD,INFO/RAW_MQ,INFO/ReadPosRankSum,INFO/SOR,INFO/VQSLOD,INFO/culprit,INFO/gnomad_3_1_1_AC,INFO/gnomad_3_1_1_AC_controls_and_biobanks,INFO/gnomad_3_1_1_AC_popmax,INFO/gnomad_3_1_1_AF,INFO/gnomad_3_1_1_AF_controls_and_biobanks,INFO/gnomad_3_1_1_AF_non_cancer,INFO/gnomad_3_1_1_AF_popmax,INFO/gnomad_3_1_1_AN,INFO/gnomad_3_1_1_AN_controls_and_biobanks,INFO/gnomad_3_1_1_AN_popmax,INFO/gnomad_3_1_1_nhomalt,INFO/gnomad_3_1_1_nhomalt_popmax,INFO/gnomad_3_1_1_primate_ai_score,INFO/gnomad_3_1_1_splice_ai_consequence"}
  output_basename: {type: 'string?', doc: "String to use as basename for outputs. Will use read1 file basename if null"}
  sample_name: {type: 'string?', doc: "Sample ID of the input reads. If not provided, will use reads1 file basename."}
  include_expression: {type: 'string?', doc: "Prefilter data frem VCF file before personal genome gen", default: STRLEN(REF)<=50 &&
      STRLEN(ALT)<=50 && FILTER="PASS"}
  subtract_bed: {type: 'File?', doc: "Supply if you want to remove regions for any reason, like low complexity or repeat mask, etc"}
  vcf_sample_name: {type: 'string?', doc: "csv string of samples if user wishes to apply filtering to and output specific samples"}
  # Genome gen vars
  genome_fa: {type: 'File?', doc: "Fasta file use for PG creation/input"}
  genomeTransformType: {type: ['null', {type: enum, name: genomeTransformType, symbols: ["None", "Haploid", "Diploid"]}], default: Diploid,
    doc: "type of genome transformation - None: no transformation. Haploid: eplace reference alleles with alternative alleles from
      VCF file (e.g. consensus allele) Diploid: create two haplotypes for each chromosome listed in VCF file, for genotypes 1—2, assumes
      perfect phasing (e.g. personal genome)"}
  gtf: {type: 'File?', doc: "Matched GTF file to index. Recommend from GENCODE, PRI assembly"}
  runThreadN: {type: 'int?', default: 32}
  memory: {type: 'int?', doc: "Mem in GB required. With no VCF, 60GB is fine, need more with VCF", default: 96}
  sjdbOverhang: {type: 'int?', default: 100, doc: "Ideal value is read len minus 1, but default 100 ok for most cases"}
  # Cutadapt 
  r1_adapter: {type: 'string?', doc: "Optional input. If the input reads have already been trimmed, leave these as null. If they do
      need trimming, supply the adapters."}
  r2_adapter: {type: 'string?', doc: "Optional input. If the input reads have already been trimmed, leave these as null. If they do
      need trimming, supply the adapters."}
  min_len: {type: 'int?', doc: "If you do not use this option, reads that have a length of zero (empty reads) are kept in the output",
    default: 20}
  quality_base: {type: 'int?', doc: "Phred scale used", default: 33}
  quality_cutoff: {type: 'int[]?', doc: "Quality trim cutoff, see https://cutadapt.readthedocs.io/en/v3.4/guide.html#quality-trimming
      for how 5' 3' is handled"}
  # STAR Diploid Align Vars
  reads1: {type: File, doc: "Input fastq file, gzipped or uncompressed OR alignment file"}
  reads2: {type: 'File?', doc: "If paired end, R2 reads files, gzipped or uncompressed"}
  samtools_fastq_cores: {type: 'int?', doc: "Num cores for align2fastq conversion, if input is an alignment file", default: 16}
  cram_reference: {type: 'File?', secondaryFiles: [.fai], doc: "If input align is cram and you are uncertain all contigs are registered
      at http://www.ebi.ac.uk/ena/cram/md5/, provide here"}
  outSAMattrRGline: {type: string, doc: "Suggested setting, with TABS SEPARATING THE TAGS, format is: ID:sample_name LB:aliquot_id
      PL:platform SM:BSID for example ID:7316-242 LB:750189 PL:ILLUMINA SM:BS_W72364MN"}
  genomeDir: {type: 'File?', doc: "Tar gzipped reference that will be unzipped at run time. Provide to skip genome generate"}
  twopassMode: {type: ['null', {type: enum, name: twopassMode, symbols: ["Basic", "None"]}], default: "Basic", doc: "Enable two pass
      mode to detect novel splice events. Default is basic (on)."}
  alignSJoverhangMin: {type: 'int?', default: 8, doc: "minimum overhang for unannotated junctions. ENCODE default used."}
  outFilterMismatchNoverLmax: {type: 'float?', default: 0.1, doc: "alignment will be output only if its ratio of mismatches to *mapped*
      length is less than or equal to this value"}
  outFilterType: {type: ['null', {type: enum, name: outFilterType, symbols: ["BySJout", "Normal"]}], default: "BySJout", doc: "type
      of filtering. Normal: standard filtering using only current alignment. BySJout (default): keep only those reads that contain
      junctions that passed filtering into SJ.out.tab."}
  outFilterScoreMinOverLread: {type: 'float?', default: 0.33, doc: "alignment will be output only if its score is higher than or equal
      to this value, normalized to read length (sum of mate's lengths for paired-end reads)"}
  outFilterMatchNminOverLread: {type: 'float?', default: 0.33, doc: "alignment will be output only if the number of matched bases
      is higher than or equal to this value., normalized to the read length (sum of mates' lengths for paired-end reads)"}
  outReadsUnmapped: {type: ['null', {type: enum, name: outReadsUnmapped, symbols: ["None", "Fastx"]}], default: "None", doc: "output
      of unmapped and partially mapped (i.e. mapped only one mate of a paired end read) reads in separate file(s). none (default):
      no output. Fastx: output in separate fasta/fastq files, Unmapped.out.mate1/2."}
  limitSjdbInsertNsj: {type: 'int?', default: 1200000, doc: "maximum number of junction to be inserted to the genome on the fly at
      the mapping stage, including those from annotations and those detected in the 1st step of the 2-pass run"}
  outSAMstrandField: {type: ['null', {type: enum, name: outSAMstrandField, symbols: ["intronMotif", "None"]}], default: "intronMotif",
    doc: "Cufflinks-like strand field flag. None: not used. intronMotif (default): strand derived from the intron motif. This option
      changes the output alignments: reads with inconsistent and/or non-canonical introns are filtered out."}
  outFilterIntronMotifs: {type: ['null', {type: enum, name: outFilterIntronMotifs, symbols: ["None", "RemoveNoncanonical", "RemoveNoncanonicalUnannotated"]}],
    default: "None", doc: "filter alignment using their motifs. None (default): no filtering. RemoveNoncanonical: filter out alignments
      that contain non-canonical junctions RemoveNoncanonicalUnannotated: filter out alignments that contain non-canonical unannotated
      junctions when using annotated splice junctions database. The annotated non-canonical junctions will be kept."}
  alignSoftClipAtReferenceEnds: {type: ['null', {type: enum, name: alignSoftClipAtReferenceEnds, symbols: ["Yes", "No"]}], default: "Yes",
    doc: "allow the soft-clipping of the alignments past the end of the chromosomes. Yes (default): allow. No: prohibit, useful for
      compatibility with Cufflinks"}
  quantMode: {type: ['null', {type: enum, name: quantMode, symbols: [TranscriptomeSAM GeneCounts, '-', TranscriptomeSAM, GeneCounts]}],
    default: TranscriptomeSAM GeneCounts, doc: "types of quantification requested. -: none. TranscriptomeSAM: output SAM/BAM alignments
      to transcriptome into a separate file GeneCounts: count reads per gene. Choices are additive, so default is 'TranscriptomeSAM
      GeneCounts'"}
  quantTranscriptomeSAMoutput: {type: ['null', {type: enum, name: quantTranscriptomeSAMoutput, symbols: [BanSingleEnd_BanIndels_ExtendSoftclip,
          BanSingleEnd, BanSingleEnd_ExtendSoftclip]}], default: BanSingleEnd_ExtendSoftclip, doc: "alignment filtering for TranscriptomeSAM
      output"}
  outSAMtype: {type: ['null', {type: enum, name: outSAMtype, symbols: ["BAM Unsorted", "None", "BAM SortedByCoordinate", "SAM Unsorted",
          "SAM SortedByCoordinate"]}], default: "BAM Unsorted", doc: "type of SAM/BAM output. None: no SAM/BAM output. Otherwise,
      first word is output type (BAM or SAM), second is sort type (Unsorted or SortedByCoordinate)"}
  outSAMunmapped: {type: ['null', {type: enum, name: outSAMunmapped, symbols: ["Within", "None", "Within KeepPairs"]}], default: "Within",
    doc: "output of unmapped reads in the SAM format. None: no output. Within (default): output unmapped reads within the main SAM
      file (i.e. Aligned.out.sam) Within KeepPairs: record unmapped mate for each alignment, and, in case of unsorted output, keep
      it adjacent to its mapped mate. Only affects multi-mapping reads"}
  genomeTransformOutput: {type: ['null', {type: enum, name: quantMode, symbols: [None, SAM, SJ, Quant, SAM SJ, SAM Quant, SAM SJ Quant,
          SJ Quant]}], default: "SAM SJ Quant", doc: "which output to transform back to original genome"}
  genomeLoad: {type: ['null', {type: enum, name: genomeLoad, symbols: ["NoSharedMemory", "LoadAndKeep", "LoadAndRemove", "LoadAndExit"]}],
    default: "NoSharedMemory", doc: "mode of shared memory usage for the genome file. In this context, the default value makes the
      most sense, the others are their as a courtesy."}
  chimMainSegmentMultNmax: {type: 'int?', doc: "maximum number of multi-alignments for the main chimeric segment. =1 will prohibit
      multimapping main segments"}
  outSAMattributes: {type: 'string?', default: "NH HI AS nM NM MD ha", doc: "a string of desired SAM attributes, in the order desired
      for the output SAM. Tags can be listed in any combination/order. Please refer to the STAR manual, as there are numerous combinations:
      https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf"}
  # fusion specific
  alignInsertionFlush: {type: ['null', {type: enum, name: alignInsertionFlush, symbols: ["None", "Right"]}], default: "None", doc: "how
      to flush ambiguous insertion positions. None (default): insertions not flushed. Right: insertions flushed to the right. STAR
      Fusion recommended (SF)"}
  alignIntronMax: {type: 'int?', default: 1000000, doc: "maximum intron size. SF recommends 100000"}
  alignMatesGapMax: {type: 'int?', default: 1000000, doc: "maximum genomic distance between mates, SF recommends 100000 to avoid readthru
      fusions within 100k"}
  alignSJDBoverhangMin: {type: 'int?', default: 1, doc: "minimum overhang for annotated junctions. SF recommends 10"}
  outFilterMismatchNmax: {type: 'int?', default: 999, doc: "maximum number of mismatches per pair, large number switches off this
      filter"}
  alignSJstitchMismatchNmax: {type: 'string?', default: "0 -1 0 0", doc: "maximum number of mismatches for stitching of the splice
      junctions. Value '5 -1 5 5' improves SF chimeric junctions, also recommended by arriba (AR)"}
  alignSplicedMateMapLmin: {type: 'int?', default: 0, doc: "minimum mapped length for a read mate that is spliced. SF recommends 30"}
  alignSplicedMateMapLminOverLmate: {type: 'float?', default: 0.66, doc: "alignSplicedMateMapLmin normalized to mate length. SF recommends
      0, AR 0.5"}
  chimJunctionOverhangMin: {type: 'int?', doc: "minimum overhang for a chimeric junction. SF recommends 8, AR 10"}
  chimMultimapNmax: {type: 'int?', default: 0, doc: "maximum number of chimeric multi-alignments. SF recommends 20, AR 50."}
  chimMultimapScoreRange: {type: 'int?', default: 1, doc: "the score range for multi-mapping chimeras below the best chimeric score.
      Only works with chimMultimapNmax > 1. SF recommends 3"}
  chimNonchimScoreDropMin: {type: 'int?', default: 20, doc: "int>=0: to trigger chimeric detection, the drop in the best non-chimeric
      alignment score with respect to the read length has to be greater than this value. SF recommends 10"}
  chimOutJunctionFormat: {type: 'int?', default: 1, doc: "formatting type for the Chimeric.out.junction file, value 1 REQUIRED for
      SF"}
  chimOutType: {type: ['null', {type: enum, name: chimOutType, symbols: ["Junctions SeparateSAMold WithinBAM SoftClip", "Junctions",
          "SeparateSAMold", "WithinBAM SoftClip", "WithinBAM HardClip", "Junctions SeparateSAMold", "Junctions WithinBAM SoftClip",
          "Junctions WithinBAM HardClip", "Junctions SeparateSAMold WithinBAM HardClip", "SeparateSAMold WithinBAM SoftClip", "SeparateSAMold
            WithinBAM HardClip"]}], doc: "type of chimeric output. Args are additive, and defined as such - Junctions: Chimeric.out.junction.
      SeparateSAMold: output old SAM into separate Chimeric.out.sam file WithinBAM: output into main aligned BAM files (Aligned.*.bam).
      WithinBAM HardClip: hard-clipping in the CIGAR for supplemental chimeric alignments WithinBAM SoftClip:soft-clipping in the
      CIGAR for supplemental chimeric alignments"}
  chimScoreDropMax: {type: 'int?', default: 20, doc: "max drop (difference) of chimeric score (the sum of scores of all chimeric segments)
      from the read length. AR recommends 30"}
  chimScoreJunctionNonGTAG: {type: 'int?', default: -1, doc: "penalty for a non-GT/AG chimeric junction. default -1, SF recommends
      -4, AR -1"}
  chimScoreSeparation: {type: 'int?', default: 10, doc: "int>=0: minimum difference (separation) between the best chimeric score and
      the next one. AR recommends 1"}
  chimSegmentMin: {type: 'int?', doc: "minimum length of chimeric segment length, if ==0, no chimeric output. REQUIRED for SF, 12
      is their default, AR recommends 10", default: 0}
  chimSegmentReadGapMax: {type: 'int?', default: 0, doc: "maximum gap in the read sequence between chimeric segments. AR recommends
      3"}
  outFilterMultimapNmax: {type: 'int?', default: 20, doc: "max number of multiple alignments allowed for a read: if exceeded, the
      read is considered unmapped. ENCODE value is default. AR recommends 50"}
  peOverlapMMp: {type: 'float?', default: 0.01, doc: "maximum proportion of mismatched bases in the overlap area. SF recommends 0.1"}
  peOverlapNbasesMin: {type: 'int?', default: 0, doc: "minimum number of overlap bases to trigger mates merging and realignment. Specify
      >0 value to switch on the 'merging of overlapping mates'algorithm. SF recommends 12,  AR recommends 10"}
  winAnchorMultimapNmax: {type: 'int?', default: 100, doc: "max number of loci anchors are allowed to map to"}
  wf_strand_param: {type: ['null', {type: 'enum', name: wf_strand_param, symbols: ["default", "rf-stranded", "fr-stranded"]}], doc: "use
      'default' for unstranded/auto, 'rf-stranded' if read1 in the fastq read pairs is reverse complement to the transcript, 'fr-stranded'
      if read1 same sense as transcript", default: "rf-stranded"}
  # RSEM inputs
  RSEMgenome: {type: 'File', doc: "RSEM reference tar ball", "sbg:suggestedValue": {class: File, path: 62853e7ad63f7c6d8d7ae5a5, name: RSEM_GENCODE39.tar.gz}}
  estimate_rspd: {type: 'boolean?', doc: "Set this option if you want to estimate the read start position distribution (RSPD) from
      data", default: true}
  #RNAseQC inputs
  RNAseQC_GTF: {type: 'File', doc: "gtf file from `gtf_anno` that has been collapsed
      GTEx-style", "sbg:suggestedValue": {class: File, path: 62853e7ad63f7c6d8d7ae5a3,
      name: gencode.v39.primary_assembly.rnaseqc.stranded.gtf}}
  RNAseQC_bed: { type: 'File?', doc: "BED file with intervals for estimating insert size distribution" }
  RNAseQC_legacy: { type: 'boolean?', doc: "If true, will output format compatible with version 1.1.9", default: false }
outputs:
  star_ref: {type: 'File?', outputSource: star_personal_genome_generate/star_ref}
  debug_log: {type: 'File?', outputSource: star_personal_genome_generate/debug_log}
  STAR_sorted_genomic_cram: {type: 'File', outputSource: samtools_bam_to_cram/output, doc: "STAR sorted and indexed genomic alignment
      cram"}
  STAR_transcriptome_bam: {type: File, outputSource: star_2-7-11b_diploid/transcriptome_bam_out}
  STAR_gene_count: {type: 'File', outputSource: star_2-7-11b_diploid/gene_counts, doc: "STAR genecounts"}
  STAR_junctions_out: {type: 'File', outputSource: star_2-7-11b_diploid/junctions_out, doc: "STARjunction reads"}
  STAR_final_log: {type: 'File', outputSource: star_2-7-11b_diploid/log_final_out, doc: "STAR metricslog file of unique, multi-mapping,
      unmapped, and chimeric reads"}
  RSEM_isoform: {type: 'File', outputSource: rsem/isoform_out, doc: "RSEM isoform expression estimates"}
  RSEM_gene: {type: 'File', outputSource: rsem/gene_out, doc: "RSEM gene expression estimates"}
  RNASeQC_Metrics: {type: 'File', outputSource: rna_seqc/Metrics, doc: "Metrics on
      mapping, intronic, exonic rates, count information, etc"}
  RNASeQC_counts: {type: 'File', outputSource: supplemental/RNASeQC_counts, doc: "Contains
      gene tpm, gene read, and exon counts"}
steps:
  basename_picker:
    run: ../tools/basename_picker.cwl
    in:
      root_name:
        source: reads1
        valueFrom: $(self.basename.split('.')[0])
      output_basename: output_basename
      sample_name: sample_name
      star_rg_line: outSAMattrRGline
    out: [outname, outsample, outrg]
  alignmentfile_pairedness:
    run: ../tools/alignmentfile_pairedness.cwl
    when: $(inputs.input_reads.basename.search(/.(b|cr|s)am$/) != -1)
    in:
      input_reads: reads1
      input_reference: cram_reference
    out: [is_paired_end]
  star_personal_genome_generate:
    run: personal_genome_input_wf.cwl
    when: $(inputs.input_genomeDir == null)
    in:
      input_vcf: input_vcf
      strip_info: strip_info
      output_basename: output_basename
      tool_name:
        valueFrom: "dna_in"
      include_expression: include_expression
      sample_name: vcf_sample_name
      subtract_bed: subtract_bed
      # Genome gen vars
      input_genomeDir: genomeDir
      genome_fa: genome_fa
      genomeTransformType: genomeTransformType
      gtf: gtf
      runThreadN: runThreadN
      memory: memory
      sjdbOverhang: sjdbOverhang
    out: [star_ref, debug_log]
  align2fastq:
    # Skip if input is FASTQ already
    run: ../tools/samtools_fastq.cwl
    when: $(inputs.input_reads_1.basename.search(/.(b|cr|s)am$/) != -1)
    in:
      input_reads_1: reads1
      SampleID: basename_picker/outname
      cores: samtools_fastq_cores
      is_paired_end: alignmentfile_pairedness/is_paired_end
      cram_reference: cram_reference
    out: [fq1, fq2]
  cutadapt_3-4:
    # Skip if no adapter given, get fastq from prev step if not null or wf input
    run: ../tools/cutadapter_3.4.cwl
    when: $(inputs.r1_adapter != null)
    in:
      readFilesIn1:
        source: [align2fastq/fq1, reads1]
        pickValue: first_non_null
      readFilesIn2:
        source: [align2fastq/fq2, reads2]
        pickValue: first_non_null
      r1_adapter: r1_adapter
      r2_adapter: r2_adapter
      min_len: min_len
      quality_base: quality_base
      quality_cutoff: quality_cutoff
      sample_name: basename_picker/outname
    out: [trimmedReadsR1, trimmedReadsR2, cutadapt_stats]
  star_2-7-11b_diploid:
    # will get fastq from first non-null in this order - cutadapt, align2fastq, wf input
    run: ../tools/star_2.7.11b_diploid_align.cwl
    in:
      outSAMattrRGline: basename_picker/outrg
      genomeDir:
        source: [star_personal_genome_generate/star_ref, genomeDir]
        pickValue: first_non_null
      readFilesIn1:
        source: [cutadapt_3-4/trimmedReadsR1, align2fastq/fq1, reads1]
        pickValue: first_non_null
      readFilesIn2:
        source: [cutadapt_3-4/trimmedReadsR2, align2fastq/fq2, reads2]
        pickValue: first_non_null
      outFileNamePrefix: basename_picker/outname
      runThreadN: runThreadN
      memory: memory
      twopassMode: twopassMode
      alignSJoverhangMin: alignSJoverhangMin
      outFilterMismatchNoverLmax: outFilterMismatchNoverLmax
      outFilterType: outFilterType
      outFilterScoreMinOverLread: outFilterScoreMinOverLread
      outFilterMatchNminOverLread: outFilterMatchNminOverLread
      outReadsUnmapped: outReadsUnmapped
      limitSjdbInsertNsj: limitSjdbInsertNsj
      outSAMstrandField: outSAMstrandField
      outFilterIntronMotifs: outFilterIntronMotifs
      alignSoftClipAtReferenceEnds: alignSoftClipAtReferenceEnds
      quantMode: quantMode
      quantTranscriptomeSAMoutput: quantTranscriptomeSAMoutput
      outSAMtype: outSAMtype
      outSAMunmapped: outSAMunmapped
      genomeTransformOutput: genomeTransformOutput
      genomeLoad: genomeLoad
      chimMainSegmentMultNmax: chimMainSegmentMultNmax
      outSAMattributes: outSAMattributes
      alignInsertionFlush: alignInsertionFlush
      alignIntronMax: alignIntronMax
      alignMatesGapMax: alignMatesGapMax
      alignSJDBoverhangMin: alignSJDBoverhangMin
      outFilterMismatchNmax: outFilterMismatchNmax
      alignSJstitchMismatchNmax: alignSJstitchMismatchNmax
      alignSplicedMateMapLmin: alignSplicedMateMapLmin
      alignSplicedMateMapLminOverLmate: alignSplicedMateMapLminOverLmate
      chimJunctionOverhangMin: chimJunctionOverhangMin
      chimMultimapNmax: chimMultimapNmax
      chimMultimapScoreRange: chimMultimapScoreRange
      chimNonchimScoreDropMin: chimNonchimScoreDropMin
      chimOutJunctionFormat: chimOutJunctionFormat
      chimOutType: chimOutType
      chimScoreDropMax: chimScoreDropMax
      chimScoreJunctionNonGTAG: chimScoreJunctionNonGTAG
      chimScoreSeparation: chimScoreSeparation
      chimSegmentMin: chimSegmentMin
      chimSegmentReadGapMax: chimSegmentReadGapMax
      outFilterMultimapNmax: outFilterMultimapNmax
      peOverlapMMp: peOverlapMMp
      peOverlapNbasesMin: peOverlapNbasesMin
      winAnchorMultimapNmax: winAnchorMultimapNmax
    out: [gene_counts, genomic_bam_out, junctions_out, log_final_out, log_out, log_progress_out, transcriptome_bam_out]
  strand_parse:
    run: ../tools/expression_parse_strand_param.cwl
    in:
      wf_strand_param: wf_strand_param
    out: [rsem_std, kallisto_std, rnaseqc_std, arriba_std]
  samtools_sort:
    run: ../tools/samtools_sort.cwl
    in:
      unsorted_bam: star_2-7-11b_diploid/genomic_bam_out
    out: [sorted_bam, sorted_bai]
  samtools_bam_to_cram:
    run: ../tools/samtools_bam_to_cram.cwl
    in:
      reference: genome_fa
      input_bam:
        source: [samtools_sort/sorted_bam, samtools_sort/sorted_bai]
        valueFrom: |
          ${
            var bundle = self[0];
            bundle.secondaryFiles = [self[1]];
            return bundle;
          }
    out: [output]
  rsem_compatibility:
    run: ../tools/rsem_compatibility.cwl
    in:
      input_bam: star_2-7-11b_diploid/transcriptome_bam_out
      output_basename: output_basename
    out: [rsem_compatible_bam]
  rsem:
    run: ../tools/rsem_calc_expression.cwl
    in:
      bam: rsem_compatibility/rsem_compatible_bam
      paired_end:
        # get value from tool, or test existence of reads2 fastq
        source: [alignmentfile_pairedness/is_paired_end, reads2]
        valueFrom: |
          ${ return self[0] != null ? self[0] : self[1] != null }
      estimate_rspd: estimate_rspd
      genomeDir: RSEMgenome
      outFileNamePrefix: basename_picker/outname
      strandedness: strand_parse/rsem_std
    out: [gene_out, isoform_out]
  rna_seqc:
    run: ../tools/rnaseqc_2.4.2.cwl
    in:
      aligned_sorted_reads: samtools_sort/sorted_bam
      collapsed_gtf: RNAseQC_GTF
      stranded: strand_parse/rnaseqc_std
      unpaired:
        source: [alignmentfile_pairedness/is_paired_end, reads2]
        valueFrom: |
          ${ return !(self[0] != null ? self[0] : self[1] != null) }
      bed: RNAseQC_bed
      legacy: RNAseQC_legacy
    out: [Metrics, Gene_TPM, Gene_count, Exon_count]
  supplemental:
    run: ../tools/supplemental_tar_gz.cwl
    in:
      outFileNamePrefix: basename_picker/outname
      Gene_TPM: rna_seqc/Gene_TPM
      Gene_count: rna_seqc/Gene_count
      Exon_count: rna_seqc/Exon_count
    out: [RNASeQC_counts]
$namespaces:
  sbg: https://sevenbridges.com
hints:
- class: "sbg:maxNumberOfParallelInstances"
  value: 3
sbg:license: Apache License 2.0
sbg:publisher: KFDRC
