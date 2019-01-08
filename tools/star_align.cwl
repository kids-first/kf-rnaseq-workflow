cwlVersion: v1.0
class: CommandLineTool
id: star_alignReads
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/star:latest'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 16
    ramMin: 60000

baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: |
      tar -zxf $(inputs.genomeDir.path)

      STAR \
      --runMode alignReads \
      --runThreadN 16 \
      --genomeDir ./STAR_index/ \
      --genomeLoad NoSharedMemory \
      --readFilesIn $(inputs.readFilesIn1.path) $(inputs.readFilesIn2.path) \
      --readFilesCommand zcat \
      --limitSjdbInsertNsj 1200000 \
      --outFileNamePrefix $(inputs.outFileNamePrefix) \
      --outSAMtype BAM Unsorted \
      --outSAMstrandField intronMotif \
      --outSAMattributes NH HI AS nM NM ch \
      --outSAMunmapped Within \
      --outSAMattrRGline $(inputs.outSAMattrRGline) \
      --outFilterType BySJout \
      --outFilterMultimapNmax 20 \
      --outFilterMismatchNmax 999 \
      --outFilterMismatchNoverLmax 0.1 \
      --outFilterScoreMinOverLread 0.33 \
      --outFilterMatchNminOverLread 0.33 \
      --alignIntronMin 20 \
      --alignIntronMax 1000000 \
      --alignMatesGapMax 1000000 \
      --alignSJDBoverhangMin 1 \
      --alignSJoverhangMin 8 \
      --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
      --chimSegmentMin 15 \
      --chimJunctionOverhangMin 15 \
      --chimMainSegmentMultNmax 1 \
      --quantMode TranscriptomeSAM GeneCounts \
      --twopassMode Basic 

inputs:
  genomeDir: File
  readFilesIn1: File[]
  readFilesIn2: File[]
  outSAMattrRGline: string
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
  star_pass1_genome: {type: Directory, outputBinding: {glob: '*_STARgenome'}}
  star_pass1: {type: Directory, outputBinding: {glob: '*_STARpass1'}}

