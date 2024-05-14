cwlVersion: v1.2
class: CommandLineTool
id: star_2-7-11b_diploid_alignReads
label: "STAR Diploid Aligner v2.7.11b"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/brownm28/star:2.7.11b'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.runThreadN)
    ramMin: ${ return inputs.memory * 1000 }

baseCommand: [tar, -I pigz, -xvf]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.genomeDir.path)
  - position: 2
    shellQuote: false
    valueFrom: >-
      && STAR
      --genomeDir ./$(inputs.genomeDir.nameroot.replace(".tar", ""))/
      --readFilesIn $(inputs.readFilesIn1.path) $(inputs.readFilesIn2 ? inputs.readFilesIn2.path : '')
      --readFilesCommand $(inputs.readFilesIn1.nameext == '.gz' ? 'zcat' : '-')
      --outFileNamePrefix $(inputs.outFileNamePrefix).
  - position: 4 
    shellQuote: false
    valueFrom: >-
      1>&2 && pigz *ReadsPerGene.out.tab *SJ.out.tab

inputs:
  outSAMattrRGline: { type: string, doc: "Suggested setting, with TABS SEPARATING \
      THE TAGS, format is: ID:sample_name LB:aliquot_id PL:platform SM:BSID for \
      example ID:7316-242 LB:750189 PL:ILLUMINA SM:BS_W72364MN",
      inputBinding: { position: 3, prefix: '--outSAMattrRGline', shellQuote: false }}
  genomeDir: { type: File, doc: "Tar gzipped reference that will be unzipped at run time" }
  readFilesIn1: { type: File, doc: "Input fastq file, gzipped or uncompressed" }
  readFilesIn2: { type: 'File?', doc: "If paired end, R2 reads files, gzipped or uncompressed" }
  outFileNamePrefix: { type: string, doc: "output files name prefix (including full or relative path). Can only be defined on the command line. \
  Tool will add '.' after prefix to easily delineate between file name and suffix" }
  runThreadN: { type: 'int?', default: 16, doc: "Adjust this value to change number of cores used.", inputBinding: { position: 3, prefix: '--runThreadN' } }
  memory: { type: 'int?', doc: "Mem in GB required. With no VCF, 60GB is fine, need more with VCF", default: 60}
  twopassMode: { type: ['null', {type: enum, name: twopassMode, symbols: ["Basic", "None"]}], default: "Basic",
  doc: "Enable two pass mode to detect novel splice events. Default is basic (on).", inputBinding: { position: 3, prefix: '--twopassMode' } }
  alignSJoverhangMin: { type: 'int?', default: 8, doc: "minimum overhang for unannotated junctions. ENCODE default used.",
  inputBinding: { position: 3, prefix: '--alignSJoverhangMin', shellQuote: false } }
  outFilterMismatchNoverLmax: { type: 'float?', default: 0.1, doc: "alignment will be output only if its ratio of mismatches to *mapped* \
  length is less than or equal to this value", inputBinding: { position: 3, prefix: '--outFilterMismatchNoverLmax' } }
  outFilterType: { type: [ 'null', {type: enum, name: outFilterType, symbols: ["BySJout", "Normal"]}], default: "BySJout",
  doc: "type of filtering. Normal: standard filtering using only current alignment. BySJout (default): keep only those reads that contain junctions \
  that passed filtering into SJ.out.tab.",
  inputBinding: { position: 3, prefix: '--outFilterType', shellQuote: false } }
  outFilterScoreMinOverLread: { type: 'float?', default: 0.33, doc: "alignment will be output only if its score is higher than or equal to this value, \
  normalized to read length (sum of mate's lengths for paired-end reads)", inputBinding: { position: 3, prefix: '--outFilterScoreMinOverLread' } }
  outFilterMatchNminOverLread: { type: 'float?', default: 0.33, doc: "alignment will be output only if the number of matched bases is higher than or \
  equal to this value., normalized to the read length (sum of mates' lengths for paired-end reads)",
  inputBinding: { position: 3, prefix: '--outFilterMatchNminOverLread' } }
  outReadsUnmapped: { type: [ 'null', {type: enum, name: outReadsUnmapped, symbols: ["None", "Fastx"]}], default: "None",
  doc: "output of unmapped and partially mapped (i.e. mapped only one mate of a paired end read) reads in separate file(s). \
  none (default): no output. Fastx: output in separate fasta/fastq files, Unmapped.out.mate1/2.",
  inputBinding: { position: 3, prefix: '--outReadsUnmapped', shellQuote: false } }
  limitSjdbInsertNsj: { type: 'int?', default: 1200000, doc: "maximum number of junction to be inserted to the genome on the fly \
  at the mapping stage, including those from annotations and those detected in the 1st step of the 2-pass run",
  inputBinding: { position: 3, prefix: '--limitSjdbInsertNsj' } }
  outSAMstrandField: { type: [ 'null', {type: enum, name: outSAMstrandField, symbols: ["intronMotif", "None"]}], default: "intronMotif",
  doc: "Cufflinks-like strand field flag. None: not used. intronMotif (default): strand derived from the intron motif. This option changes the output \
  alignments: reads with inconsistent and/or non-canonical introns are filtered out.",
  inputBinding: { position: 3, prefix: '--outSAMstrandField', shellQuote: false } }
  outFilterIntronMotifs: { type: [ 'null', {type: enum, name: outFilterIntronMotifs, symbols: ["None", "RemoveNoncanonical", "RemoveNoncanonicalUnannotated"]}],
  default: "None",
  doc: "filter alignment using their motifs. None (default): no filtering. RemoveNoncanonical: filter out alignments that contain non-canonical junctions \
  RemoveNoncanonicalUnannotated: filter out alignments that contain non-canonical unannotated junctions when using annotated splice junctions database. \
  The annotated non-canonical junctions will be kept.",
  inputBinding: { position: 3, prefix: '--outFilterIntronMotifs', shellQuote: false } }
  alignSoftClipAtReferenceEnds:  { type: [ 'null', {type: enum, name: alignSoftClipAtReferenceEnds, symbols: ["Yes", "No"]}], default: "Yes",
  doc: "allow the soft-clipping of the alignments past the end of the chromosomes. Yes (default): allow. \
  No: prohibit, useful for compatibility with Cufflinks",
  inputBinding: { position: 3, prefix: '--alignSoftClipAtReferenceEnds', shellQuote: false } }
  quantMode: { type: [ 'null', {type: enum, name: quantMode, symbols: [TranscriptomeSAM GeneCounts, -, TranscriptomeSAM, GeneCounts]}],
  default: TranscriptomeSAM GeneCounts,
  doc: "types of quantification requested. -: none. TranscriptomeSAM: output SAM/BAM alignments to transcriptome into a separate file \
  GeneCounts: count reads per gene. Choices are additive, so default is 'TranscriptomeSAM GeneCounts'",
  inputBinding: { position: 3, prefix: '--quantMode', shellQuote: false } }
  quantTranscriptomeSAMoutput: { type: [ 'null', {type: enum, name: quantTranscriptomeSAMoutput, symbols: [BanSingleEnd_BanIndels_ExtendSoftclip, BanSingleEnd, BanSingleEnd_ExtendSoftclip]}],
  default: BanSingleEnd_BanIndels_ExtendSoftclip,
  doc: "alignment filtering for TranscriptomeSAM output",
  inputBinding: { position: 3, prefix: '--quantTranscriptomeSAMoutput', shellQuote: false } }
  outSAMtype: { type: [ 'null', {type: enum, name: outSAMtype, symbols: ["BAM Unsorted", "None", "BAM SortedByCoordinate", "SAM Unsorted", "SAM SortedByCoordinate"]}],
  default: "BAM Unsorted",
  doc: "type of SAM/BAM output. None: no SAM/BAM output. Otherwise, first word is output type (BAM or SAM), second is sort type (Unsorted or SortedByCoordinate)",
  inputBinding: { position: 3, prefix: '--outSAMtype', shellQuote: false} }
  outSAMunmapped: { type: [ 'null', {type: enum, name: outSAMunmapped, symbols: ["Within", "None", "Within KeepPairs"]}],
  default: "Within",
  doc: "output of unmapped reads in the SAM format. None: no output. Within (default): output unmapped reads within the main SAM file (i.e. Aligned.out.sam) \
  Within KeepPairs: record unmapped mate for each alignment, and, in case of unsorted output, keep it adjacent to its mapped mate. Only affects \
  multi-mapping reads",
  inputBinding: { position: 3, prefix: '--outSAMunmapped', shellQuote: false } }
  genomeTransformOutput: { type: [ 'null', {type: enum, name: quantMode, symbols: [None, SAM, SJ, Quant, SAM SJ, SAM Quant, SAM SJ Quant, SJ Quant ]}],
  default: None,
  doc: "which output to transform back to original genome",
  inputBinding: { position: 3, prefix: '--genomeTransformOutput', shellQuote: false } }
  genomeLoad: { type: [ 'null', {type: enum, name: genomeLoad, symbols: ["NoSharedMemory", "LoadAndKeep", "LoadAndRemove", "LoadAndExit"]}],
  default: "NoSharedMemory",
  doc: "mode of shared memory usage for the genome file. In this context, the default value makes the most sense, the others are their as a courtesy.",
  inputBinding: { position: 3, prefix: '--genomeLoad', shellQuote: false } }
  chimMainSegmentMultNmax: { type: 'int?', doc: "maximum number of multi-alignments for the main chimeric segment. =1 will prohibit multimapping main segments",
  inputBinding: { position: 3, prefix: '--chimMainSegmentMultNmax' } }
  outSAMattributes: { type: 'string?', default: 'NH HI AS nM NM ch ha', doc: "a string of desired SAM attributes, in the order desired for the output SAM. Tags can be listed in any combination/order. \
  Please refer to the STAR manual, as there are numerous combinations: https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf",
  inputBinding: { position: 3, prefix: '--outSAMattributes', shellQuote: false } }
  # fusion specific
  alignInsertionFlush: { type: [ 'null', {type: enum, name: alignInsertionFlush, symbols: ["None", "Right"]}], default: "None",
  doc: "how to flush ambiguous insertion positions. None (default): insertions not flushed. Right: insertions flushed to the right.
  STAR Fusion recommended (SF)",
  inputBinding: { position: 3, prefix: '--alignInsertionFlush', shellQuote: false } }
  alignIntronMax: { type: 'int?', default: 1000000, doc: "maximum intron size. SF recommends 100000", inputBinding: { position: 3, prefix: '--alignIntronMax' } }
  alignMatesGapMax: { type: 'int?', default: 1000000, doc: "maximum genomic distance between mates, SF recommends 100000 \
  to avoid readthru fusions within 100k",
  inputBinding: { position: 3, prefix: '--alignMatesGapMax' } }
  alignSJDBoverhangMin: { type: 'int?', default: 1, doc: "minimum overhang for annotated junctions. SF recommends 10",
  inputBinding: { position: 3, prefix: '--alignSJDBoverhangMin' } }
  outFilterMismatchNmax: { type: 'int?', default: 999,  doc: "maximum number of mismatches per pair, large number switches off this filter",
  inputBinding: { position: 3, prefix: '--outFilterMismatchNmax' } }
  alignSJstitchMismatchNmax: { type: 'string?', default: "0 -1 0 0", doc: "maximum number of mismatches for stitching of the splice junctions. \
  Value '5 -1 5 5' improves SF chimeric junctions, also recommended by arriba (AR)",
  inputBinding: { position: 3, prefix: '--alignSJstitchMismatchNmax', shellQuote: false } }
  alignSplicedMateMapLmin: { type: 'int?', default: 0, doc: "minimum mapped length for a read mate that is spliced. SF recommends 30",
  inputBinding: { position: 3, prefix: '--alignSplicedMateMapLmin' } }
  alignSplicedMateMapLminOverLmate: { type: 'float?', default: 0.66,
  doc: "alignSplicedMateMapLmin normalized to mate length. SF recommends 0, AR 0.5",
  inputBinding: { position: 3, prefix: '--alignSplicedMateMapLminOverLmate' } }
  chimJunctionOverhangMin: { type: 'int?', doc: "minimum overhang for a chimeric junction. SF recommends 8, AR 10",
  inputBinding: { position: 3, prefix: '--chimJunctionOverhangMin' } }
  chimMultimapNmax: { type: 'int?', default: 0, doc: "maximum number of chimeric multi-alignments. SF recommends 20, AR 50.",
  inputBinding: { position: 3, prefix: '--chimMultimapNmax' } }
  chimMultimapScoreRange: { type: 'int?', default: 1, doc: "the score range for multi-mapping chimeras below the best chimeric \
  score. Only works with chimMultimapNmax > 1. SF recommends 3",
  inputBinding: { position: 3, prefix: '--chimMultimapScoreRange' } }
  chimNonchimScoreDropMin: { type: 'int?', default: 20,
  doc: "int>=0: to trigger chimeric detection, the drop in the best non-chimeric \
  alignment score with respect to the read length has to be greater than this value. SF recommends 10",
  inputBinding: { position: 3, prefix: '--chimNonchimScoreDropMin' } }
  chimOutJunctionFormat: { type: 'int?', default: 1, doc: "formatting type for the Chimeric.out.junction file, value 1 REQUIRED for SF",
  inputBinding: { position: 3, prefix: '--chimOutJunctionFormat' } }
  chimOutType: { type: [ 'null', {type: enum, name: chimOutType, symbols: [
      "Junctions SeparateSAMold WithinBAM SoftClip",
      "Junctions", "SeparateSAMold",
      "WithinBAM SoftClip",
      "WithinBAM HardClip",
      "Junctions SeparateSAMold",
      "Junctions WithinBAM SoftClip",
      "Junctions WithinBAM HardClip",
      "Junctions SeparateSAMold WithinBAM HardClip",
      "SeparateSAMold WithinBAM SoftClip",
      "SeparateSAMold WithinBAM HardClip"
      ]}],
  doc: "type of chimeric output. Args are additive, and defined as such - Junctions: Chimeric.out.junction. SeparateSAMold: output old SAM into separate Chimeric.out.sam file \
  WithinBAM: output into main aligned BAM files (Aligned.*.bam). WithinBAM HardClip: hard-clipping in the CIGAR for supplemental chimeric alignments \
  WithinBAM SoftClip:soft-clipping in the CIGAR for supplemental chimeric alignments",
  inputBinding: { position: 3, prefix: '--chimOutType', shellQuote: false } }
  chimScoreDropMax: { type: 'int?', default: 20,
  doc: "max drop (difference) of chimeric score (the sum of scores of all chimeric segments) from the read length. AR recommends 30",
  inputBinding: { position: 3, prefix: '--chimScoreDropMax' } }
  chimScoreJunctionNonGTAG: { type: 'int?', default: -1, doc: "penalty for a non-GT/AG chimeric junction. \
  default -1, SF recommends -4, AR -1",
  inputBinding: { position: 3, prefix: '--chimScoreJunctionNonGTAG' } }
  chimScoreSeparation: { type: 'int?', default: 10,
  doc: "int>=0: minimum difference (separation) between the best chimeric score and the next one. AR recommends 1",
  inputBinding: { position: 3, prefix: '--chimScoreSeparation' } }
  chimSegmentMin: { type: 'int?', doc: "minimum length of chimeric segment length, if ==0, no chimeric output. \
  REQUIRED for SF, 12 is their default, AR recommends 10",
  inputBinding: { position: 3, prefix: '--chimSegmentMin' } }
  chimSegmentReadGapMax: { type: 'int?', default: 0, doc: "maximum gap in the read sequence between chimeric segments. AR recommends 3",
  inputBinding: { position: 3, prefix: '--chimSegmentReadGapMax' } }
  outFilterMultimapNmax: { type: 'int?', default: 20, doc: "max number of multiple alignments allowed for \
  a read: if exceeded, the read is considered unmapped. ENCODE value is default. AR recommends 50",
  inputBinding: { position: 3, prefix: '--outFilterMultimapNmax' } }
  peOverlapMMp: { type: 'float?', default: 0.01, doc: "maximum proportion of mismatched bases in the overlap area. SF recommends 0.1",
  inputBinding: { position: 3, prefix: '--peOverlapMMp' } }
  peOverlapNbasesMin: { type: 'int?', default: 0,
  doc: "minimum number of overlap bases to trigger mates merging and realignment. Specify >0 value to switch \
  on the 'merging of overlapping mates'algorithm. SF recommends 12,  AR recommends 10",
  inputBinding: { position: 3, prefix: '--peOverlapNbasesMin' } }
  winAnchorMultimapNmax: { type: 'int?', default: 100,
  doc: "max number of loci anchors are allowed to map to",
  inputBinding: { position: 3, prefix: '--winAnchorMultimapNmax' } }

outputs:
  log_progress_out: { type: File, doc: "Simple progress output. Can use to gauge speed and run time", outputBinding: {glob: '*Log.progress.out'} }
  log_out: { type: File, doc: "Contains a summary of all params used and reference files", outputBinding: {glob: '*Log.out'} }
  log_final_out: { type: File, doc: "Overall summary of read mapping statistics", outputBinding: {glob: '*Log.final.out'} }
  genomic_bam_out: { type: File, doc: "UNSORTED read mapping to genomic coordinates", outputBinding: {glob: '*Aligned.out.bam'} }
  junctions_out: { type: File, doc: "high confidence collapsed splice junctions in tab-delimited form", outputBinding: {glob: '*SJ.out.tab.gz'} }
  transcriptome_bam_out: { type: File, doc: "Read mapping to transcriptome", outputBinding: {glob: '*Aligned.toTranscriptome.out.bam'} }
  gene_counts: { type: File, doc: "STAR-generated read counts by gene", outputBinding: {glob: '*ReadsPerGene.out.tab.gz'} }

$namespaces:
  sbg: https://sevenbridges.com

hints: [
  {
      "class": "sbg:SaveLogs",
      "value": "*Log.*"
  }
]
