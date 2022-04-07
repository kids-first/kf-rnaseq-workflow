cwlVersion: v1.2
class: Workflow
id: kfdrc-rnaseq-workflow-v4
label: Kids First DRC RNAseq Workflow Version 4
doc: "Latest and the greatest"

requirements:
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: SubworkflowFeatureRequirement

inputs:
  # many tool
  output_basename: { type: 'string', doc: "String to use as basename for outputs" }
  reads1: { type: File, doc: "Input fastq file, gzipped or uncompressed OR bam file" }
  reads2: { type: 'File?', doc: "If paired end, R2 reads files, gzipped or uncompressed" }

  # wf_strand_param: { type: [{type: 'enum', name: wf_strand_param, symbols: ["default",
  #         "rf-stranded", "fr-stranded"]}], doc: "use 'default' for unstranded/auto,\
  #     \ 'rf-stranded' if read1 in the fastq read pairs is reverse complement to the\
  #     \ transcript, 'fr-stranded' if read1 same sense as transcript" }
  # gtf_anno: { type: 'File', doc: "General transfer format (gtf) file with gene models corresponding to fasta reference" }
  # star_fusion_genome_untar_path: {type: 'string?', doc: "This is what the path will be\
  #     \ when genome_tar is unpackaged", default: "GRCh38_v39_CTAT_lib_Mar242022.CUSTOM"}

  # samtools fastq
  # samtools_fastq_cores: { type: 'int?', doc: "Num cores for bam2fastq conversion, if input is bam", default: 36 }
  # input_type: {type: [{type: 'enum', name: input_type, symbols: ["PEBAM", "SEBAM",
  #         "FASTQ"]}], doc: "Please select one option for input file type, PEBAM (paired-end\
  #     \ BAM), SEBAM (single-end BAM) or FASTQ."}
  # # cutadapt
  # r1_adapter: { type: 'string?', doc: "Optional input. If the input reads have already\
  #     \ been trimmed, leave these as null. If they do need trimming, supply the adapters." }
  # r2_adapter: { type: 'string?', doc: "Optional input. If the input reads have already\
  #     \ been trimmed, leave these as null. If they do need trimming, supply the adapters." }
  # STAR
  outSAMattrRGline: { type: string, doc: "Suggested setting, with TABS SEPARATING \
      THE TAGS, format is: ID:sample_name LB:aliquot_id PL:platform SM:BSID for \
      example ID:7316-242 LB:750189 PL:ILLUMINA SM:BS_W72364MN" }
  STARgenome: { type: File, doc: "Tar gzipped reference that will be unzipped at run time" }
  runThreadN: { type: 'int?', default: 16, doc: "Adjust this value to change number of cores used." }
  # twopassMode: { type: ['null', {type: enum, name: twopassMode, symbols: ["Basic", "None"]}], default: "Basic",
  # doc: "Enable two pass mode to detect novel splice events. Default is basic (on)." }
  # alignSJoverhangMin: { type: 'int?', default: 8, doc: "minimum overhang for unannotated junctions. ENCODE default used." }
  # outFilterMismatchNoverLmax: { type: 'float?', default: 0.1, doc: "alignment will be output only if its ratio of mismatches to *mapped* \
  # length is less than or equal to this value" }
  # outFilterType: { type: [ 'null', {type: enum, name: outFilterType, symbols: ["BySJout", "Normal"]}], default: "BySJout",
  # doc: "type of filtering. Normal: standard filtering using only current alignment. BySJout (default): keep only those reads that contain junctions \
  # that passed filtering into SJ.out.tab." }
  # outFilterScoreMinOverLread: { type: 'float?', default: 0.33, doc: "alignment will be output only if its score is higher than or equal to this value, \
  # normalized to read length (sum of mate's lengths for paired-end reads)" }
  # outFilterMatchNminOverLread: { type: 'float?', default: 0.33, doc: "alignment will be output only if the number of matched bases is higher than or \
  # equal to this value., normalized to the read length (sum of mates' lengths for paired-end reads)" }
  # outReadsUnmapped: { type: [ 'null', {type: enum, name: outReadsUnmapped, symbols: ["None", "Fastx"]}], default: "None",
  # doc: "output of unmapped and partially mapped (i.e. mapped only one mate of a paired end read) reads in separate file(s). \
  # none (default): no output. Fastx: output in separate fasta/fastq files, Unmapped.out.mate1/2." }
  # limitSjdbInsertNsj: { type: 'int?', default: 1200000, doc: "maximum number of junction to be inserted to the genome on the fly \
  # at the mapping stage, including those from annotations and those detected in the 1st step of the 2-pass run" }
  # outSAMstrandField: { type: [ 'null', {type: enum, name: outSAMstrandField, symbols: ["intronMotif", "None"]}], default: "intronMotif",
  # doc: "Cufflinks-like strand field flag. None: not used. intronMotif (default): strand derived from the intron motif. This option changes the output \
  # alignments: reads with inconsistent and/or non-canonical introns are filtered out." }
  # outFilterIntronMotifs: { type: [ 'null', {type: enum, name: outFilterIntronMotifs, symbols: ["None", "RemoveNoncanonical", "RemoveNoncanonicalUnannotated"]}],
  # default: "None",
  # doc: "filter alignment using their motifs. None (default): no filtering. RemoveNoncanonical: filter out alignments that contain non-canonical junctions \
  # RemoveNoncanonicalUnannotated: filter out alignments that contain non-canonical unannotated junctions when using annotated splice junctions database. \
  # The annotated non-canonical junctions will be kept." }
  # alignSoftClipAtReferenceEnds:  { type: [ 'null', {type: enum, name: alignSoftClipAtReferenceEnds, symbols: ["Yes", "No"]}], default: "Yes",
  # doc: "allow the soft-clipping of the alignments past the end of the chromosomes. Yes (default): allow. \
  # No: prohibit, useful for compatibility with Cufflinks" }
  # quantMode: { type: [ 'null', {type: enum, name: quantMode, symbols: [TranscriptomeSAM GeneCounts, -, TranscriptomeSAM, GeneCounts]}],
  # default: TranscriptomeSAM GeneCounts,
  # doc: "types of quantification requested. -: none. TranscriptomeSAM: output SAM/BAM alignments to transcriptome into a separate file \
  # GeneCounts: count reads per gene. Choices are additive, so default is 'TranscriptomeSAM GeneCounts'" }
  # outSAMtype: { type: [ 'null', {type: enum, name: outSAMtype, symbols: ["BAM Unsorted", "None", "BAM SortedByCoordinate", "SAM Unsorted", "SAM SortedByCoordinate"]}],
  # default: "BAM Unsorted",
  # doc: "type of SAM/BAM output. None: no SAM/BAM output. Otherwise, first word is output type (BAM or SAM), second is sort type (Unsorted or SortedByCoordinate)" }
  # outSAMunmapped: { type: [ 'null', {type: enum, name: outSAMunmapped, symbols: ["Within", "None", "Within KeepPairs"]}],
  # default: "Within",
  # doc: "output of unmapped reads in the SAM format. None: no output. Within (default): output unmapped reads within the main SAM file (i.e. Aligned.out.sam) \
  # Within KeepPairs: record unmapped mate for each alignment, and, in case of unsorted output, keep it adjacent to its mapped mate. Only affects \
  # multi-mapping reads" }
  # genomeLoad: { type: [ 'null', {type: enum, name: genomeLoad, symbols: ["NoSharedMemory", "LoadAndKeep", "LoadAndRemove", "LoadAndExit"]}],
  # default: "NoSharedMemory",
  # doc: "mode of shared memory usage for the genome file. In this context, the default value makes the most sense, the others are their as a courtesy." }
  # chimMainSegmentMultNmax: { type: 'int?', default: 1, doc: "maximum number of multi-alignments for the main chimeric segment. =1 will prohibit multimapping main segments" }
  # outSAMattributes: { type: 'string?', default: 'NH HI AS nM NM ch', doc: "a string of desired SAM attributes, in the order desired for the output SAM. Tags can be listed in any combination/order. \
  # Please refer to the STAR manual, as there are numerous combinations: https://raw.githubusercontent.com/alexdobin/star_2.7.10a/master/doc/STARmanual.pdf" }
  # alignInsertionFlush: { type: [ 'null', {type: enum, name: alignInsertionFlush, symbols: ["None", "Right"]}], default: "None",
  # doc: "how to flush ambiguous insertion positions. None (default): insertions not flushed. Right: insertions flushed to the right.
  # STAR Fusion recommended (SF)" }
  # alignIntronMax: { type: 'int?', default: 1000000, doc: "maximum intron size. SF recommends 100000" }
  # alignMatesGapMax: { type: 'int?', default: 1000000, doc: "maximum genomic distance between mates, SF recommends 100000 \
  # to avoid readthru fusions within 100k" }
  # alignSJDBoverhangMin: { type: 'int?', default: 1, doc: "minimum overhang for annotated junctions. SF recommends 10" }
  # outFilterMismatchNmax: { type: 'int?', default: 999,  doc: "maximum number of mismatches per pair, large number switches off this filter" }
  # alignSJstitchMismatchNmax: { type: 'string?', default: "0 -1 0 0", doc: "maximum number of mismatches for stitching of the splice junctions. \
  # Value '5 -1 5 5' improves SF chimeric junctions, also recommended by arriba (AR)" }
  # alignSplicedMateMapLmin: { type: 'int?', default: 0, doc: "minimum mapped length for a read mate that is spliced. SF recommends 30" }
  # alignSplicedMateMapLminOverLmate: { type: 'float?', default: 0.66,
  # doc: "alignSplicedMateMapLmin normalized to mate length. SF recommends 0, AR 0.5" }
  # chimJunctionOverhangMin: { type: 'int?', default: 15, doc: "minimum overhang for a chimeric junction. SF recommends 8, AR 10" }
  # chimMultimapNmax: { type: 'int?', default: 0, doc: "maximum number of chimeric multi-alignments. SF recommends 20, AR 50." }
  # chimMultimapScoreRange: { type: 'int?', default: 1, doc: "the score range for multi-mapping chimeras below the best chimeric \
  # score. Only works with chimMultimapNmax > 1. SF recommends 3" }
  # chimNonchimScoreDropMin: { type: 'int?', default: 20,
  # doc: "int>=0: to trigger chimeric detection, the drop in the best non-chimeric \
  # alignment score with respect to the read length has to be greater than this value. SF recommends 10" }
  # chimOutJunctionFormat: { type: 'int?', default: 1, doc: "formatting type for the Chimeric.out.junction file, value 1 REQUIRED for SF" }
  # chimOutType: { type: [ 'null', {type: enum, name: chimOutType, symbols: [
  #     "Junctions SeparateSAMold WithinBAM SoftClip",
  #     "Junctions", "SeparateSAMold",
  #     "WithinBAM SoftClip",
  #     "WithinBAM HardClip",
  #     "Junctions SeparateSAMold",
  #     "Junctions WithinBAM SoftClip",
  #     "Junctions WithinBAM HardClip",
  #     "Junctions SeparateSAMold WithinBAM HardClip",
  #     "SeparateSAMold WithinBAM SoftClip",
  #     "SeparateSAMold WithinBAM HardClip"
  #     ]}],
  # default: "Junctions SeparateSAMold WithinBAM SoftClip",
  # doc: "type of chimeric output. Args are additive, and defined as such - Junctions: Chimeric.out.junction. SeparateSAMold: output old SAM into separate Chimeric.out.sam file \
  # WithinBAM: output into main aligned BAM files (Aligned.*.bam). WithinBAM HardClip: hard-clipping in the CIGAR for supplemental chimeric alignments \
  # WithinBAM SoftClip:soft-clipping in the CIGAR for supplemental chimeric alignments" }
  # chimScoreDropMax: { type: 'int?', default: 20,
  # doc: "max drop (difference) of chimeric score (the sum of scores of all chimeric segments) from the read length. AR recommends 30" }
  # chimScoreJunctionNonGTAG: { type: 'int?', default: -1, doc: "penalty for a non-GT/AG chimeric junction. \
  # default -1, SF recommends -4, AR -1" }
  # chimScoreSeparation: { type: 'int?', default: 10,
  # doc: "int>=0: minimum difference (separation) between the best chimeric score and the next one. AR recommends 1" }
  # chimSegmentMin: { type: 'int?', default: 12, doc: "minimum length of chimeric segment length, if ==0, no chimeric output. \
  # REQUIRED for SF, 12 is their default, AR recommends 10" }
  # chimSegmentReadGapMax: { type: 'int?', default: 0, doc: "maximum gap in the read sequence between chimeric segments. AR recommends 3" }
  # outFilterMultimapNmax: { type: 'int?', default: 20, doc: "max number of multiple alignments allowed for \
  # a read: if exceeded, the read is considered unmapped. ENCODE value is default. AR recommends 50" }
  # peOverlapMMp: { type: 'float?', default: 0.01, doc: "maximum proportion of mismatched bases in the overlap area. SF recommends 0.1" }
  # peOverlapNbasesMin: { type: 'int?', default: 0,
  # doc: "minimum number of overlap bases to trigger mates merging and realignment. Specify >0 value to switch \
  # on the 'merging of overlapping mates'algorithm. SF recommends 12,  AR recommends 10" }

  # arriba
  # reference_fasta: {type: 'File', doc: "GRCh38.primary_assembly.genome.fa", "sbg:suggestedValue": {
  #     class: File, path: 5f500135e4b0370371c051b4, name: GRCh38.primary_assembly.genome.fa}}
  # arriba_memory: { type: 'int?', doc: "Mem intensive tool. Set in GB", default: 64 }
  # # STAR Fusion
  # FusionGenome: { type: 'File', doc: "STAR-Fusion CTAT Genome lib" }
  # compress_chimeric_junction: { type: 'boolean?', default: true,
  # doc: 'If part of a workflow, recommend compressing this file as final output' }
  # RNAseQC_GTF: {type: 'File', doc: "gtf file from `gtf_anno` that has been collapsed GTEx-style"}
  # # kallisto
  # kallisto_idx: { type: 'File', doc: "Specialized index of a **transcriptome** fasta file for kallisto" }
  # kallisto_avg_frag_len:  {type: 'int?', doc: "Optional input. Average fragment length\
  #     \ for Kallisto only if single end input." }
  # kallisto_std_dev: { type: 'long?', doc: "Optional input. Standard Deviation of the\
  #     \ average fragment length for Kallisto only needed if single end input." }
  # # RSEM
  # RSEMgenome: {type: 'File', doc: "RSEM_GENCODE27.tar.gz", "sbg:suggestedValue": {
  #     class: File, path: 5f500135e4b0370371c051be, name: RSEM_GENCODE27.tar.gz}}
  # paired_end: {type: 'boolean?', doc: "If input is paired-end, add this flag", default: true }
  # estimate_rspd: { type: 'boolean?',
  # doc: "Set this option if you want to estimate the read start position distribution (RSPD) from data", default: false }
  # annoFuse
  # sample_name: { type: 'string', doc: "Sample ID of the input reads" }
  # annofuse_col_num: {type: 'int?', doc: "column number in file of fusion name."}
  # rmats
  # rmats_read_length: {type: 'int', doc: "Input read length for sample reads."}
  # rmats_variable_read_length: {type: 'boolean?', doc: "Allow reads with lengths that\
  #     \ differ from --readLength to be processed. --readLength will still be used\
  #     \ to determine IncFormLen and SkipFormLen."}
  # rmats_novel_splice_sites: {type: 'boolean?', doc: "Select for novel splice site\
  #     \ detection or unannotated splice sites. 'true' to detect or add this parameter,\
  #     \ 'false' to disable denovo detection. Tool Default: false"}
  # rmats_stat_off: {type: 'boolean?', doc: "Select to skip statistical analysis, either\
  #     \ between two groups or on single sample group. 'true' to add this parameter.\
  #     \ Tool default: false"}
  # rmats_allow_clipping: {type: 'boolean?', doc: "Allow alignments with soft or hard\
  #     \ clipping to be used."}
  # rmats_threads: {type: 'int?', doc: "Threads to allocate to RMATs."}
  # rmats_ram: {type: 'int?', doc: "GB of RAM to allocate to RMATs."}
  # rmats_read_type: { type: [ 'null', {type: enum, name: rmats_read_type, symbols: ["single", "paired"]}], default: "paired",
  # doc: "Indicate whether input reads are single- or paired-end" }

outputs:
  # cutadapt_stats: {type: 'File?', outputSource: cutadapt/cutadapt_stats, doc: "Cutadapt\
  #     \ stats output, only if adapter is supplied."}
  STAR_transcriptome_bam: {type: 'File', outputSource: star_2.7.10a/transcriptome_bam_out,
    doc: "STAR bam of transcriptome reads"}
  # STAR_sorted_genomic_bam: {type: 'File', outputSource: samtools_sort/sorted_bam,
  #   doc: "STAR sorted alignment bam"}
  # STAR_sorted_genomic_bai: {type: 'File', outputSource: samtools_sort/sorted_bai,
  #   doc: "STAR index for sorted aligned bam"}
  # STAR_chimeric_junctions: {type: 'File?', outputSource: star_fusion_1.10.1/chimeric_junction_compressed,
  #   doc: "STAR chimeric junctions"}
  STAR_gene_count: {type: 'File', outputSource: star_2.7.10a/gene_counts, doc: "STAR gene\
     counts"}
  STAR_junctions_out: {type: 'File', outputSource: star_2.7.10a/junctions_out, doc: "STAR\
     junction reads"}
  STAR_final_log: {type: 'File', outputSource: star_2.7.10a/log_final_out, doc: "STAR metrics\
     log file of unique, multi-mapping, unmapped, and chimeric reads"}
  # STAR-Fusion_results: {type: 'File', outputSource: star_fusion_1.10.1/abridged_coding, doc: "STAR\
  #     \ fusion detection from chimeric reads"}
  # arriba_fusion_results: {type: 'File', outputSource: arriba_fusion_2.2.1/arriba_fusions,
  #   doc: "Fusion output from Arriba"}
  # arriba_fusion_viz: {type: 'File', outputSource: arriba_draw_2.2.1/arriba_pdf, doc: "pdf\
  #     \ output from Arriba"}
  # RSEM_isoform: {type: 'File', outputSource: rsem/isoform_out, doc: "RSEM isoform\
  #     \ expression estimates"}
  # RSEM_gene: {type: 'File', outputSource: rsem/gene_out, doc: "RSEM gene expression\
  #     \ estimates"}
  # RNASeQC_Metrics: {type: 'File', outputSource: rna_seqc/Metrics, doc: "Metrics on\
  #     \ mapping, intronic, exonic rates, count information, etc"}
  # RNASeQC_counts: {type: 'File', outputSource: supplemental/RNASeQC_counts, doc: "Contains\
  #     \ gene tpm, gene read, and exon counts"}
  # kallisto_Abundance: {type: 'File', outputSource: kallisto/abundance_out, doc: "Gene\
  #     \ abundance output from STAR genomic bam file"}
  # # annofuse_filtered_fusions_tsv: {type: 'File?', outputSource: annofuse/annofuse_filtered_fusions_tsv,
  # #   doc: "Filtered fusions called by annoFuse."}
  # rmats_filtered_alternative_3_prime_splice_sites_jc: {type: 'File', outputSource: rmats/filtered_alternative_3_prime_splice_sites_jc,
  #   doc: "Alternative 3 prime splice sites JC.txt output from RMATs containing only\
  #     \ those calls with 10 or more read counts of support"}
  # rmats_filtered_alternative_5_prime_splice_sites_jc: {type: 'File', outputSource: rmats/filtered_alternative_5_prime_splice_sites_jc,
  #   doc: "Alternative 5 prime splice sites JC.txt output from RMATs containing only\
  #     \ those calls with 10 or more read counts of support"}
  # rmats_filtered_mutually_exclusive_exons_jc: {type: 'File', outputSource: rmats/filtered_mutually_exclusive_exons_jc,
  #   doc: "Mutually exclusive exons JC.txt output from RMATs containing only those\
  #     \ calls with 10 or more read counts of support"}
  # rmats_filtered_retained_introns_jc: {type: 'File', outputSource: rmats/filtered_retained_introns_jc,
  #   doc: "Retained introns JC.txt output from RMATs containing only those calls with\
  #     \ 10 or more read counts of support"}
  # rmats_filtered_skipped_exons_jc: {type: 'File', outputSource: rmats/filtered_skipped_exons_jc,
  #   doc: "Skipped exons JC.txt output from RMATs containing only those calls with\
  #     \ 10 or more read counts of support"}

steps:
  # bam2fastq:
  #   # Skip if input is FASTQ already
  #   run: ../tools/samtools_fastq.cwl
  #   when: $(inputs.input_type != "FASTQ")
  #   in:
  #     input_reads_1: reads1
  #     input_reads_2: reads2
  #     SampleID: output_basename
  #     cores: samtools_fastq_cores
  #     input_type: input_type
  #   out: [fq1, fq2]

  # cutadapt:
  #   # Skip if no adapter given, get fastq from prev step if not null or wf input
  #   run: ../tools/cutadapter.cwl
  #   when: $(inputs.r1_adapter != null)
  #   in:
  #     readFilesIn1: 
  #       source: [bam2fastq/fq1, reads1]
  #       pickValue: first_non_null
  #     readFilesIn2:
  #       source: [bam2fastq/fq2, reads2]
  #       pickValue: first_non_null
  #     r1_adapter: r1_adapter
  #     r2_adapter: r2_adapter
  #     sample_name: output_basename
  #   out: [trimmedReadsR1, trimmedReadsR2, cutadapt_stats]

  star_2.7.10a:
    # will get fastq from first non-null in this order - cutadapt, bam2fastq, wf input
    run: ../tools/star_2.7.10a_align.cwl
    in:
      outSAMattrRGline: outSAMattrRGline
      genomeDir: STARgenome
      readFilesIn1: reads1
        # source: [cutadapt/trimmedReadsR1, bam2fastq/fq1, reads1]
        # pickValue: first_non_null
      readFilesIn2: reads2
      #   source: [cutadapt/trimmedReadsR1, bam2fastq/fq2, reads2]
      #   pickValue: first_non_null
      outFileNamePrefix: output_basename
      runThreadN: runThreadN
      # twopassMode: twopassMode
      # alignSJoverhangMin: alignSJoverhangMin
      # outFilterMismatchNoverLmax: outFilterMismatchNoverLmax
      # outFilterType: outFilterType
      # outFilterScoreMinOverLread: outFilterScoreMinOverLread
      # outFilterMatchNminOverLread: outFilterMatchNminOverLread
      # outReadsUnmapped: outReadsUnmapped
      # limitSjdbInsertNsj: limitSjdbInsertNsj
      # outSAMstrandField: outSAMstrandField
      # outFilterIntronMotifs: outFilterIntronMotifs
      # alignSoftClipAtReferenceEnds: alignSoftClipAtReferenceEnds
      # quantMode: quantMode
      # outSAMtype: outSAMtype
      # outSAMunmapped: outSAMunmapped
      # genomeLoad: genomeLoad
      # chimMainSegmentMultNmax: chimMainSegmentMultNmax
      # outSAMattributes: outSAMattributes
      # alignInsertionFlush: alignInsertionFlush
      # alignIntronMax: alignIntronMax
      # alignMatesGapMax: alignMatesGapMax
      # alignSJDBoverhangMin: alignSJDBoverhangMin
      # outFilterMismatchNmax: outFilterMismatchNmax
      # alignSJstitchMismatchNmax: alignSJstitchMismatchNmax
      # alignSplicedMateMapLmin: alignSplicedMateMapLmin
      # alignSplicedMateMapLminOverLmate: alignSplicedMateMapLminOverLmate
      # chimJunctionOverhangMin: chimJunctionOverhangMin
      # chimMultimapNmax: chimMultimapNmax
      # chimMultimapScoreRange: chimMultimapScoreRange
      # chimNonchimScoreDropMin: chimNonchimScoreDropMin
      # chimOutJunctionFormat: chimOutJunctionFormat
      # chimOutType: chimOutType
      # chimScoreDropMax: chimScoreDropMax
      # chimScoreJunctionNonGTAG: chimScoreJunctionNonGTAG
      # chimScoreSeparation: chimScoreSeparation
      # chimSegmentMin: chimSegmentMin
      # chimSegmentReadGapMax: chimSegmentReadGapMax
      # outFilterMultimapNmax: outFilterMultimapNmax
      # peOverlapMMp: peOverlapMMp
      # peOverlapNbasesMin: peOverlapNbasesMin
    out: [chimeric_junctions, chimeric_sam_out, gene_counts, genomic_bam_out, junctions_out,
      log_final_out, log_out, log_progress_out, transcriptome_bam_out]

  # samtools_sort:
  #   run: ../tools/samtools_sort.cwl
  #   in:
  #     unsorted_bam: star_2.7.10a/genomic_bam_out
  #     chimeric_sam_out: star_2.7.10a/chimeric_sam_out
  #   out: [sorted_bam, sorted_bai, chimeric_bam_out]

  # rmats:
  #   run: ../workflow/rmats_wf.cwl
  #   in:
  #     gtf_annotation: gtf_anno
  #     sample_1_bams:
  #       source: samtools_sort/sorted_bam
  #       valueFrom: |
  #         $([self])
  #     read_length: rmats_read_length
  #     variable_read_length: rmats_variable_read_length
  #     read_type: rmats_read_type
  #     strandedness:
  #       source: wf_strand_param
  #       valueFrom: |
  #         $(self == "rf-stranded" ? "fr-firststrand" : self == "fr-stranded" ? "fr-secondstrand" : "fr-unstranded")
  #     novel_splice_sites: rmats_novel_splice_sites
  #     stat_off: rmats_stat_off
  #     allow_clipping: rmats_allow_clipping
  #     output_basename: output_basename
  #     rmats_threads: rmats_threads
  #     rmats_ram: rmats_ram
  #   out: [filtered_alternative_3_prime_splice_sites_jc, filtered_alternative_5_prime_splice_sites_jc,
  #     filtered_mutually_exclusive_exons_jc, filtered_retained_introns_jc, filtered_skipped_exons_jc]

  # strand_parse:
  #   run: ../tools/expression_parse_strand_param.cwl
  #   in:
  #     wf_strand_param: wf_strand_param
  #   out: [rsem_std, kallisto_std, rnaseqc_std, arriba_std]

  # star_fusion_1.10.1:
  #   run: ../tools/star_fusion_1.10.1_call.cwl
  #   in:
  #     Chimeric_junction: star_2.7.10a/chimeric_junctions
  #     genome_tar: FusionGenome
  #     output_basename: output_basename
  #     genome_untar_path: star_fusion_genome_untar_path
  #     compress_chimeric_junction: compress_chimeric_junction
  #   out: [abridged_coding, chimeric_junction_compressed]

  # arriba_fusion_2.2.1:
  #   run: ../tools/arriba_fusion_2.2.1.cwl
  #   in:
  #     genome_aligned_bam:
  #       source: [samtools_sort/sorted_bam, samtools_sort/sorted_bai]
  #       valueFrom: |
  #         ${
  #           var bundle = self[0];
  #           bundle.secondaryFiles = [self[1]];
  #           return bundle;
  #         }
  #     memory: arriba_memory
  #     reference_fasta: reference_fasta
  #     gtf_anno: gtf_anno
  #     outFileNamePrefix: output_basename
  #     arriba_strand_flag: strand_parse/arriba_std
  #   out: [arriba_fusions]

  # arriba_draw_2.2.1:
  #   run: ../tools/arriba_draw_2.2.1.cwl
  #   in:
  #     fusions: arriba_fusion_2.2.1/arriba_fusions
  #     genome_aligned_bam:
  #       source: [samtools_sort/sorted_bam, samtools_sort/sorted_bai]
  #       valueFrom: |
  #         ${
  #           var bundle = self[0];
  #           bundle.secondaryFiles = [self[1]];
  #           return bundle;
  #         }
  #     gtf_anno: gtf_anno
  #   out: [arriba_pdf]

  # rsem:
  #   run: ../tools/rsem_calc_expression.cwl
  #   in:
  #     bam: star_2.7.10a/transcriptome_bam_out
  #     paired_end: paired_end
  #     estimate_rspd: estimate_rspd
  #     genomeDir: RSEMgenome
  #     outFileNamePrefix: output_basename
  #     strandedness: strand_parse/rsem_std
  #   out: [gene_out, isoform_out]

  # rna_seqc:
  #   run: ../tools/rnaseqc_2.4.2.cwl
  #   in:
  #     Aligned_sorted_bam: samtools_sort/sorted_bam
  #     collapsed_gtf: RNAseQC_GTF
  #     stranded: strand_parse/rnaseqc_std
  #     unpaired:
  #       source: rmats_read_type
  #       valueFrom: |
  #         $(self == "single" ? true : false)
  #   out: [Metrics, Gene_TPM, Gene_count, Exon_count]

  # supplemental:
  #   run: ../tools/supplemental_tar_gz.cwl
  #   in:
  #     outFileNamePrefix: output_basename
  #     Gene_TPM: rna_seqc/Gene_TPM
  #     Gene_count: rna_seqc/Gene_count
  #     Exon_count: rna_seqc/Exon_count
  #   out: [RNASeQC_counts]

  # kallisto:
  #   run: ../tools/kallisto_calc_expression.cwl
  #   in:
  #     transcript_idx: kallisto_idx
  #     strand: strand_parse/kallisto_std
  #     reads1:
  #       source: [cutadapt/trimmedReadsR1, bam2fastq/fq1, reads1]
  #       pickValue: first_non_null
  #     reads2:
  #       source: [cutadapt/trimmedReadsR2, bam2fastq/fq2, reads2]
  #       pickValue: first_non_null
  #     SampleID: output_basename
  #     avg_frag_len: kallisto_avg_frag_len
  #     std_dev: kallisto_std_dev
  #   out: [abundance_out]

  # annofuse:
  #   run: ../workflow/kfdrc_annoFuse_wf.cwl
  #   in:
  #     sample_name: sample_name
  #     FusionGenome: FusionGenome
  #     genome_untar_path: star_fusion_genome_untar_path
  #     rsem_expr_file: rsem/gene_out
  #     arriba_output_file: arriba_fusion_2.2.1/arriba_fusions
  #     star_fusion_output_file: star_fusion_1.10.1/abridged_coding
  #     col_num: annofuse_col_num
  #     output_basename: output_basename
  #   out: [annofuse_filtered_fusions_tsv]

$namespaces:
  sbg: https://sevenbridges.com
hints:
- class: "sbg:maxNumberOfParallelInstances"
  value: 3
"sbg:license": Apache License 2.0
"sbg:publisher": KFDRC
"sbg:categories":
- ALIGNMENT
- ANNOFUSE
- ARRIBA
- BAM
- CUTADAPT
- FASTQ
- KALLISTO
- PE
- RNASEQ
- RNASEQC
- RSEM
- SE
- STAR
"sbg:links":
- id: 'https://github.com/kids-first/kf-rnaseq-workflow/releases/tag/v4.0.0'
  label: github-release
