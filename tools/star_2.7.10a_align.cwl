cwlVersion: v1.0
class: CommandLineTool
id: star_2.7.10a_alignReads
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'migbro/star:2.7.10a'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 16
    ramMin: 60000

baseCommand: [tar,, -I pigz, -xvf]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.genomeDir.path)
  - position: 2
    shellQuote: false
    valueFrom: >-
      && STAR
      --genomeDir ./$(inputs.genomeDir.nameroot.split('.')[0])/
      --readFilesIn $(inputs.readFilesIn1.path) $(inputs.readFilesIn2 ? inputs.readFilesIn2.path : '')
      --readFilesCommand $(inputs.readFilesIn1.nameext == 'gz' ? 'zcat' : '-')
      --outFileNamePrefix $(inputs.outFileNamePrefix).
      --outSAMattributes NH HI AS nM NM ch
  - position: 4 
    shellQuote: false
    valueFrom: >-
      && pigz *ReadsPerGene.out.tab *SJ.out.tab

inputs:
  outSAMattrRGline: { type: string, doc: "Suggested setting, with TABS SEPARATING\
      \ THE TAGS, format is: ID:sample_name LB:aliquot_id PL:platform SM:BSID for\
      \ example ID:7316-242 LB:750189 PL:ILLUMINA SM:BS_W72364MN",
      inputBinding: { position: 3, prefix: '--outSAMattrRGline' }}
  genomeDir: { type: File, doc: "Tar gzipped reference that will be nuipped at run time" }
  readFilesIn1: { type: File, doc: "Input fastq file, gzipped or uncompressed" }
  readFilesIn2: { type: 'File?', doc: "If paired end, R2 reads files, gzipped or uncompressed" }
  outFileNamePrefix: { type: string, doc: "output files name prefix (including full or relative path). Can only be defined on the command line. \
  Tool will add '.' after prefix to easily delineate between file name and suffix"}
  runThreadN: { type: 'int?', default: 16, doc: "Adjust this value to change number of cores used.", inputBinding: { position: 3, prefix: '--runThreadN' } }
  twopassMode: { type: [{type: enum, name: twopassMode, symbols: ["Basic", "None"]}],
  doc: "Enable two pass mode to detect novel splice events. Default is basic (on).", inputBinding: { position: 3, prefix: '--twopassMode' }}
  outFilterMultimapNmax: { type: 'int?', default: 20, doc: "max number of multiple alignments allowed for \
  a read: if exceeded, the read is considered unmapped. ENCODE value is default", inputBinding: { position: 3, prefix: '--outFilterMultimapNmax' } }
  alignSJoverhangMin: { type: 'int?', default: 8, doc: "minimum overhang for unannotated junctions. ENCODE default used.",
  inputBinding: { position: 3, prefix: '--alignSJoverhangMin' } }
  alignSJDBoverhangMin: { type: 'int?', default: 1, doc: "minimum overhang for annotated junctions.", inputBinding: { position: 3, prefix: '--alignSJDBoverhangMin' } }
  alignSJstitchMismatchNmax { type: 'string?', default: "5 -1 5 5", doc: "maximum number of mismatches for stitching of the splice junctions. \
  Default improves STAR Fusion chimeric junctions", inputBinding: { position: 3, prefix: '--alignSJstitchMismatchNmax' } }
  outFilterMismatchNmax: { type: 'int?', default: 999,  doc: "maximum number of mismatches per pair, large number switches off this filter",
  inputBinding: { position: 3, prefix: '--outFilterMismatchNmax' } }
  outFilterMismatchNoverLmax: { type: 'float?', default: 0.1, doc: "alignment will be output only if its ratio of mismatches to *mapped* \
  length is less than or equal to this value", inputBinding: { position: 3, prefix: '--outFilterMismatchNoverLmax' } }
  alignIntronMax: { type: 'int?', default: 100000, doc: "maximum intron size", inputBinding: { position: 3, prefix: '--alignIntronMax' } }
  chimSegmentReadGapMax: { type: 'int?', default: 0, doc: "maximum gap in the read sequence between chimeric segments",
  inputBinding: { position: 3, prefix: '--chimSegmentReadGapMax' } }
  chimOutJunctionFormat: { type: 'int?', default: 1, doc: "formatting type for the Chimeric.out.junction file, value 1 required for STAR Fusion",
  inputBinding: { position: 3, prefix: '--chimOutJunctionFormat' } }
  alignMatesGapMax: { type: 'int?', default: 100000, doc: "maximum genomic distance between mates, default set according to STAR Fusion \
  to avoid readthru fusions within 100k",
  inputBinding: { position: 3, prefix: '--alignMatesGapMax' } }
  outFilterType: { type: [{type: enum, name: outFilterType, symbols: ["BySJout", "Normal"]}],
  doc: "type of filtering. Normal: standard filtering using only current alignment. BySJout (default): keep only those reads that contain junctions \
  that passed filtering into SJ.out.tab.",
  inputBinding: { position: 3, prefix: '--outFilterType' } }
  outFilterScoreMinOverLread: { type: 'float?', default: 0.33, doc: "alignment will be output only if its score is higher than or equal to this value, \
   normalized to read length (sum of mate's lengths for paired-end reads)", inputBinding: { position: 3, prefix: '--outFilterScoreMinOverLread' } }
  outFilterMatchNminOverLread: { type: 'float?', default: 0.33, doc: "alignment will be output only if the number of matched bases is higher than or equal to this value., \
   normalized to the read length (sum of mates' lengths for paired-end reads)", inputBinding: { position: 3, prefix: '--outFilterMatchNminOverLread' } }
  outReadsUnmapped: { type: [{type: enum, name: outReadsUnmapped, symbols: ["None", "Fastx"]}],
  doc: "output of unmapped and partially mapped (i.e. mapped only one mate of a paired end read) reads in separate file(s). \
  none (default): no output. Fastx: output in separate fasta/fastq files, Unmapped.out.mate1/2.",
  inputBinding: { position: 3, prefix: '--outReadsUnmapped' } }
  limitSjdbInsertNsj: { type: 'int?', default: 1200000, doc: "maximum number of junction to be inserted to the genome on the fly \
  at the mapping stage, including those from annotations and those detected in the 1st step of the 2-pass run",
  inputBinding: { position: 3, prefix: '--limitSjdbInsertNsj' } }
  outSAMstrandField: { type: [{type: enum, name: outSAMstrandField, symbols: ["intronMotif", "None"]}],},
  doc: "Cufflinks-like strand field flag. None: not used. intronMotif (default): strand derived from the intron motif. This option changes the output \
  alignments: reads with inconsistent and/or non-canonical introns are filtered out.",
  inputBinding: { position: 3, prefix: '--outSAMstrandField' } }
  outFilterIntronMotifs: { type: [{type: enum, name: outFilterIntronMotifs, symbols: ["None", "RemoveNoncanonical", "RemoveNoncanonicalUnannotated"]}],},
  doc: "filter alignment using their motifs. None (default): no filtering. RemoveNoncanonical: filter out alignments that contain non-canonical junctions \
  RemoveNoncanonicalUnannotated: filter out alignments that contain non-canonical unannotated junctions when using annotated splice junctions database. \
  The annotated non-canonical junctions will be kept.",
  inputBinding: { position: 3, prefix: '--outFilterIntronMotifs' } }
  alignSoftClipAtReferenceEnds:  { type: [{type: enum, name: alignSoftClipAtReferenceEnds, symbols: ["Yes", "No"]}],},
  doc: "allow the soft-clipping of the alignments past the end of the chromosomes. Yes (default): allow. \
  No: prohibit, useful for compatibility with Cufflinks",
  inputBinding: { position: 3, prefix: '--alignSoftClipAtReferenceEnds' } }
  quantMode { type: [{type: enum, name: quantMode, symbols: ["TranscriptomeSAM GeneCounts", "-", "TranscriptomeSAM", "GeneCounts"]}],},
  doc: "types of quantification requested. -: none. TranscriptomeSAM: output SAM/BAM alignments to transcriptome into a separate file \
  GeneCounts: count reads per gene. Choices are additive, so default is 'TranscriptomeSAM GeneCounts'",
  inputBinding: { position: 3, prefix: '--quantMode' } }
  outSAMtype  { type: [{type: enum, name: outSAMtype, symbols: ["BAM Unsorted", "None", "BAM SortedByCoordinate", "SAM Unsorted", "SAM SortedByCoordinate"]}],},
  doc: "type of SAM/BAM output. None: no SAM/BAM output. Otherwise, first word is output type (BAM or SAM), second is sort type (Unsorted or SortedByCoordinate)",
  inputBinding: { position: 3, prefix: '--outSAMtype' } }
  outSAMunmapped { type: [{type: enum, name: outSAMunmapped, symbols: ["Within", "None", "Within KeepPairs"]}],},
  doc: "output of unmapped reads in the SAM format. None: no output. Within (default): output unmapped reads within the main SAM file (i.e. Aligned.out.sam) \
  Within KeepPairs: record unmapped mate for each alignment, and, in case of unsorted output, keep it adjacent to its mapped mate. Only affects \
  multi-mapping reads",
  inputBinding: { position: 3, prefix: '--outSAMunmapped' } }
  genomeLoad: { type: [{type: enum, name: genomeLoad, symbols: ["NoSharedMemory", "LoadAndKeep", "LoadAndRemove", "LoadAndExit"]}],},
  doc: "mode of shared memory usage for the genome file. In this context, the default value makes the most sense, the others are their as a courtesy.",
  inputBinding: { position: 3, prefix: '--genomeLoad' } }
  chimSegmentMin: { type: 'int?', default: 12, doc: "minimum length of chimeric segment length, if ==0, no chimeric output. Needed for Star Fusion, 12 is their default",
  inputBinding: { position: 3, prefix: '--chimSegmentMin' } }
  chimJunctionOverhangMin: {type: 'int?', default: 15, doc: "minimum overhang for a chimeric junction",
  inputBinding: { position: 3, prefix: '--chimJunctionOverhangMin' }
  chimOutType: { type: [{type: enum, name: chimOutType, symbols: [
      "Junctions SeparateSAMold WithinBAM SoftClip",
      "Junctions",
      "SeparateSAMold",
      "WithinBAM SoftClip",
      "WithinBAM HardClip",
      "Junctions SeparateSAMold",
      "Junctions WithinBAM SoftClip".
      "Junctions WithinBAM HardClip",
      "Junctions SeparateSAMold WithinBAM HardClip",
      "SeparateSAMold WithinBAM SoftClip",
      "SeparateSAMold WithinBAM HardClip"
      ]}],},
  doc: "type of chimeric output. Args are additive, and defined as such - Junctions: Chimeric.out.junction. SeparateSAMold: output old SAM into separate Chimeric.out.sam file \
  WithinBAM: output into main aligned BAM files (Aligned.*.bam). WithinBAM HardClip: hard-clipping in the CIGAR for supplemental chimeric alignments \
  WithinBAM SoftClip:soft-clipping in the CIGAR for supplemental chimeric alignments",
  inputBinding: { position: 3, prefix: '--chimOutType' }
  chimMainSegmentMultNmax: { type: 'int?', default: 1, doc: "maximum number of multi-alignments for the main chimeric segment. =1 will prohibit multimapping main segments",
  inputBinding: { position: 3, prefix: '--chimMainSegmentMultNmax' }
  outSAMattributes: { type: 'string?', default: 'NH HI AS nM NM ch', doc: "a string of desired SAM attributes, in the order desired for the output SAM. Tags can be listed in any combination/order. \
  Please refer to the STAR manual, as their are numerous combinations: https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf",
  inputBinding: { position: 3, prefix: '--outSAMattributes' }

outputs:
  log_progress_out: {type: File, outputBinding: {glob: '*Log.progress.out'}}
  log_out: {type: File, outputBinding: {glob: '*Log.out'}}
  log_final_out: {type: File, outputBinding: {glob: '*Log.final.out'}}
  genomic_bam_out: {type: File, outputBinding: {glob: '*Aligned.out.bam'}}
  junctions_out: {type: File, outputBinding: {glob: '*SJ.out.tab.gz'}}
  transcriptome_bam_out: {type: File, outputBinding: {glob: '*Aligned.toTranscriptome.out.bam'}}
  chimeric_sam_out: {type: File, outputBinding: {glob: '*Chimeric.out.sam'}}
  chimeric_junctions: {type: File, outputBinding: {glob: '*Chimeric.out.junction'}}
  gene_counts: {type: File, outputBinding: {glob: '*ReadsPerGene.out.tab.gz'}}
