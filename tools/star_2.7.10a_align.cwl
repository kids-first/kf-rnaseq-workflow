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
      --outFilterMultimapNmax 20
      --alignSJoverhangMin 8
      --alignSJDBoverhangMin 10 *
      --alignSJstitchMismatchNmax 5 -1 5 5 *
      --outFilterMismatchNmax 999
      --outFilterMismatchNoverLmax 0.1
      --alignIntronMax 100000
      --chimSegmentReadGapMax 3
      --chimOutJunctionFormat 1
      --alignMatesGapMax 100000 *
      --outFilterType BySJout
      --outFilterScoreMinOverLread 0.33
      --outFilterMatchNminOverLread 0.33
      --outReadsUnmapped None
      --limitSjdbInsertNsj 1200000
      --outFileNamePrefix $(inputs.outFileNamePrefix).
      --outSAMstrandField intronMotif
      --outFilterIntronMotifs None
      --alignSoftClipAtReferenceEnds Yes
      --quantMode TranscriptomeSAM GeneCounts
      --outSAMtype BAM Unsorted
      --outSAMunmapped Within
      --genomeLoad NoSharedMemory
      --chimSegmentMin 12 *
      --chimJunctionOverhangMin 12
      --chimOutType Junctions SeparateSAMold WithinBAM SoftClip
      --chimMainSegmentMultNmax 1
      --outSAMattributes NH HI AS nM NM ch &&
      gzip *ReadsPerGene.out.tab *SJ.out.tab

inputs:
  outSAMattrRGline: { type: string, doc: "Suggested setting, with TABS SEPARATING\
      \ THE TAGS, format is: ID:sample_name LB:aliquot_id PL:platform SM:BSID for\
      \ example ID:7316-242 LB:750189 PL:ILLUMINA SM:BS_W72364MN",
      inputBinding: { position: 3, prefix: '--outSAMattrRGline' }}
  genomeDir: { type: File, doc: "Tar gzipped reference that will be nuipped at run time" }
  readFilesIn1: { type: File, doc: "Input fastq file, gzipped or uncompressed" }
  readFilesIn2: { type: 'File?', doc: "If paired end, R2 reads files, gzipped or uncompressed" }
  runThreadN: { type: 'int?', default: 16, doc: "Adjust this value to change number of cores used.", inputBinding: { position: 3, prefix: '--runThreadN' } }
  twopassMode: { type: [{type: enum, name: twopassMode, symbols: ["Basic", "None"]}], doc: "Enable two pass mode to detect novel splice events. Default is basic (on).", inputBinding: { position: 3, prefix: '--twopassMode' }}
  outFilterMultimapNmax: { type: 'int?', default: 20, doc: "max number of multiple alignments allowed for \
  a read: if exceeded, the read is considered unmapped. ENCODE value is default", inputBinding: { position: 3, prefix: '--outFilterMultimapNmax' } }
  alignSJoverhangMin: { type: 'int?', default: 8, doc: "minimum overhang for unannotated junctions. ENCODE default used.", inputBinding: { position: 3, prefix: '--alignSJoverhangMin' } }
  alignSJDBoverhangMin: { type: 'int?', default: 1, doc: "minimum overhang for annotated junctions.", inputBinding: { position: 3, prefix: '--alignSJDBoverhangMin' } }
  alignSJstitchMismatchNmax { type: 'string?', default: "5 -1 5 5", doc: "maximum number of mismatches for stitching of the splice junctions. \
  Default improves STAR Fusion chimeric junctions", inputBinding: { position: 3, prefix: '--alignSJstitchMismatchNmax' }}
    --outFilterMismatchNmax 999
    --outFilterMismatchNoverLmax 0.1
    --alignIntronMax 100000 *
    --chimSegmentReadGapMax 3
    --chimOutJunctionFormat 1 *
    --alignMatesGapMax 100000
    --outFilterType BySJout
    --outFilterScoreMinOverLread 0.33
    --outFilterMatchNminOverLread 0.33
    --outReadsUnmapped None
    --limitSjdbInsertNsj 1200000
    --outFileNamePrefix $(inputs.outFileNamePrefix).
    --outSAMstrandField intronMotif
    --outFilterIntronMotifs None
    --alignSoftClipAtReferenceEnds Yes
    --quantMode TranscriptomeSAM GeneCounts
    --outSAMtype BAM Unsorted
    --outSAMunmapped Within
    --genomeLoad NoSharedMemory
    --chimSegmentMin 12
    --chimJunctionOverhangMin 8 *
    --chimOutType Junctions SeparateSAMold WithinBAM SoftClip
    --chimMainSegmentMultNmax 1
    --outSAMattributes NH HI AS nM NM ch
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
