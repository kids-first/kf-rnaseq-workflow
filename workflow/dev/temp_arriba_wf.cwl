cwlVersion: v1.0
class: Workflow
id: temp_arriba_only_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  sample_name: string
  wf_strand_param: {type: ['null', string], doc: "use 'default' for unstranded/auto, rf_stranded if read1 in the fastq read pairs is reverse complement to the transcript, fr-stranded if read1 same sense as transcript"}
  reference_fasta: File
  gtf_anno: File
  genome_aligned_bam: File
  genome_aligned_bai: File
  chimeric_sam_out: File

outputs:
  arriba_fusion_results: {type: File, outputSource: arriba_fusion/arriba_fusions}
  arriba_fusion_viz: {type: File, outputSource: arriba_fusion/arriba_pdf}

steps:
  strand_parse:
    run: ../../tools/expression_parse_strand_param.cwl
    in:
      strand: wf_strand_param
    out:
      [
        rsem_std,
        kallisto_std,
        rnaseqc_std,
        arriba_std
      ]
  arriba_fusion:
    run: ../../tools/arriba_fusion.cwl
    in:
      genome_aligned_bam: genome_aligned_bam
      genome_aligned_bai: genome_aligned_bai
      chimeric_sam_out: chimeric_sam_out
      reference_fasta: reference_fasta
      gtf_anno: gtf_anno
      outFileNamePrefix: sample_name
      arriba_strand_flag: strand_parse/arriba_std
    out:
      [
        arriba_fusions,
        arriba_pdf
      ]
