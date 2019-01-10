cwlVersion: v1.0
class: CommandLineTool
id: star_alignReads
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/star:latest'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: "$(inputs.runThreadN ? inputs.runThreadN : 16)"
    ramMin: 60000

baseCommand: [tar, -xzf]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.genomeDir.path) &&
      STAR --outSAMattrRGline $(inputs.outSAMattrRGline)
      --genomeDir ./$(inputs.genomeDir.nameroot.split('.')[0])/
      --readFilesIn $(inputs.readFilesIn1.path) $(inputs.readFilesIn2.path)
      --readFilesCommand zcat
      --runThreadN $(inputs.runThreadN)
      --twopassMode Basic
      --outFilterMultimapNmax 20
      --alignSJoverhangMin 8
      --alignSJDBoverhangMin 1
      --outFilterMismatchNmax 999
      --outFilterMismatchNoverLmax 0.1
      --alignIntronMin 20
      --alignIntronMax 1000000
      --alignMatesGapMax 1000000
      --outFilterType BySJout
      --outFilterScoreMinOverLread 0.33
      --outFilterMatchNminOverLread 0.33
      --limitSjdbInsertNsj 1200000
      --outFileNamePrefix $(inputs.outFileNamePrefix).
      --outSAMstrandField intronMotif
      --outFilterIntronMotifs None
      --alignSoftClipAtReferenceEnds Yes
      --quantMode TranscriptomeSAM GeneCounts
      --outSAMtype BAM Unsorted
      --outSAMunmapped Within
      --genomeLoad NoSharedMemory
      --chimSegmentMin 15
      --chimJunctionOverhangMin 15
      --chimOutType Junctions SeparateSAMold WithinBAM SoftClip
      --chimMainSegmentMultNmax 1
      --outSAMattributes NH HI AS nM NM ch

inputs:
  outSAMattrRGline: string
  readFilesIn1: File
  readFilesIn2: File
  genomeDir: File
  runThreadN: int
  outFileNamePrefix: string

outputs:
  log_progress_out: {type: File, outputBinding: {glob: '*Log.progress.out'}}
  log_out: {type: File, outputBinding: {glob: '*Log.out'}}
  log_final_out: {type: File, outputBinding: {glob: '*Log.final.out'}}
  genomic_bam_out: {type: File, outputBinding: {glob: '*Aligned.out.bam'}}
  junctions_out: {type: File, outputBinding: {glob: '*SJ.out.tab'}}
  transcriptome_bam_out: {type: File, outputBinding: {glob: '*Aligned.toTranscriptome.out.bam'}}
  chimeric_sam_out: {type: File, outputBinding: {glob: '*Chimeric.out.sam'}}
  chimeric_junctions: {type: File, outputBinding: {glob: '*Chimeric.out.junction'}}
  gene_counts: {type: File, outputBinding: {glob: '*ReadsPerGene.out.tab'}}
