cwlVersion: v1.0
class: CommandLineTool
id: star_alignReads
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/star:2.6.1d'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 16
    ramMin: 60000

baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      tar -xzf $(inputs.genomeDir.path) &&
      STAR
      --outSAMattrRGline $(inputs.outSAMattrRGline)
      --genomeDir ./$(inputs.genomeDir.nameroot.split('.')[0])/
      --readFilesIn $(inputs.readFilesIn1.path) $(inputs.readFilesIn2 ? inputs.readFilesIn2.path : '')
      --readFilesCommand $(inputs.readFilesIn1.nameext == '.gz' ? 'zcat' : '-')
      --runThreadN $(inputs.runThreadN)
      --twopassMode Basic
      --outFilterMultimapNmax 20
      --alignSJoverhangMin 8
      --alignSJDBoverhangMin 10
      --alignSJstitchMismatchNmax 5 -1 5 5
      --outFilterMismatchNmax 999
      --outFilterMismatchNoverLmax 0.1
      --alignIntronMax 100000
      --chimSegmentReadGapMax 3
      --chimOutJunctionFormat 1
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
      --chimJunctionOverhangMin 12
      --chimOutType Junctions SeparateSAMold WithinBAM SoftClip
      --chimMainSegmentMultNmax 1
      --outSAMattributes NH HI AS nM NM ch &&
      gzip *ReadsPerGene.out.tab *SJ.out.tab

inputs:
  outSAMattrRGline: string
  readFilesIn1: File
  readFilesIn2: File?
  genomeDir: File
  runThreadN: int
  outFileNamePrefix: string

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
