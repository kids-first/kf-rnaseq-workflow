# Dockers of kfdrc_RNAseq_workflow.cwl

TOOL|DOCKER
-|-
alignmentfile_pairedness.cwl|quay.io/biocontainers/pysam:0.22.0--py310h41dec4a_0
annoFuse.cwl|pgc-images.sbgenomics.com/d3b-bixu/annofuse:0.92.0
arriba_draw_2.2.1.cwl|pgc-images.sbgenomics.com/d3b-bixu/arriba:2.2.1
arriba_fusion_2.2.1.cwl|pgc-images.sbgenomics.com/d3b-bixu/arriba:2.2.1
awk_junction_filtering.cwl|ubuntu:20.04
bam_strandness.cwl|pgc-images.sbgenomics.com/d3b-bixu/stranded:1.1.0
basename_picker.cwl|None
build_reads_record.cwl|None
cutadapter_3.4.cwl|pgc-images.sbgenomics.com/d3b-bixu/cutadapt:3.4
expression_parse_strand_param.cwl|None
format_arriba_fusion_file.cwl|pgc-images.sbgenomics.com/d3b-bixu/annofuse:0.92.0
fusion_annotator.cwl|pgc-images.sbgenomics.com/d3b-bixu/fusionannotator:0.1.1
kallisto_calc_expression.cwl|images.sbgenomics.com/uros_sipetic/kallisto:0.43.1
rmats_both_bam.cwl|xinglab/rmats:v4.1.2
rnaseqc_2.4.2.cwl|gcr.io/broad-cga-aarong-gtex/rnaseqc:2.4.2
rsem_calc_expression.cwl|images.sbgenomics.com/uros_sipetic/rsem:1.3.1
samtools_bam_to_cram.cwl|pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9
samtools_fastq.cwl|pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9
samtools_head.cwl|staphb/samtools:1.20
samtools_readlength_bam.cwl|pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9
samtools_sort.cwl|pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9
samtools_split.cwl|staphb/samtools:1.20
star_2.7.10a_align.cwl|pgc-images.sbgenomics.com/d3b-bixu/star:2.7.10a
star_fusion_1.10.1_call.cwl|pgc-images.sbgenomics.com/d3b-bixu/star:fusion-1.10.1
supplemental_tar_gz.cwl|pgc-images.sbgenomics.com/d3b-bixu/ubuntu:18.04
t1k.cwl|pgc-images.sbgenomics.com/d3b-bixu/t1k:v1.0.5
