cwlVersion: v1.2
class: Workflow
id: fusion-wf
label: Fusion Workflow
doc: |
  # Fusion Workflow

  This workflow is a sub workflow for calling RNA fusions. For full documentation on usage, inputs, outputs, and defaults please see doc section of Kids First DRC RNAseq Workflow.

requirements:
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: SubworkflowFeatureRequirement
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement
inputs:
  reference_fasta: {type: 'File', doc: "GRCh38.primary_assembly.genome.fa", "sbg:suggestedValue": {class: File, path: 5f500135e4b0370371c051b4,
      name: GRCh38.primary_assembly.genome.fa, secondaryFiles: [{class: File, path: 62866da14d85bc2e02ba52db, name: GRCh38.primary_assembly.genome.fa.fai}]},
    secondaryFiles: ['.fai']}
  output_basename: {type: 'string?', doc: "String to use as basename for outputs. Will use read1 file basename if null"}
  gtf_anno: {type: 'File', doc: "General transfer format (gtf) file with gene models corresponding to fasta reference", "sbg:suggestedValue": {
      class: File, path: 62853e7ad63f7c6d8d7ae5a4, name: gencode.v39.primary_assembly.annotation.gtf}}
  genome_untar_path: {type: 'string?', doc: "This is what the path will be when genome_tar is unpackaged", default: "GRCh38_v39_CTAT_lib_Mar242022.CUSTOM"}
  arriba_memory: {type: 'int?', doc: "Mem intensive tool. Set in GB", default: 64}
  FusionGenome: {type: 'File', doc: "STAR-Fusion CTAT Genome lib", "sbg:suggestedValue": {class: File, path: 62853e7ad63f7c6d8d7ae5a8,
      name: GRCh38_v39_CTAT_lib_Mar242022.CUSTOM.tar.gz}}
  compress_chimeric_junction: {type: 'boolean?', default: true, doc: 'If part of a workflow, recommend compressing this file as final
      output'}
  rsem_expr_file: {type: 'File'}
  annofuse_col_num: {type: 'int?', doc: "0-based column number in file of fusion name.", default: 30}
  fusion_annotator_ref: {type: 'File', doc: "Tar ball with fusion_annot_lib.idx and blast_pairs.idx from STAR-Fusion CTAT Genome lib.
      Can be same as FusionGenome, but only two files needed from that package", "sbg:suggestedValue": {class: 'File', path: '63cff818facdd82011c8d6fe',
      name: 'GRCh38_v39_fusion_annot_custom.tar.gz'}}
  sample_name: {type: 'string'}
  genome_aligned_bam: {type: 'File', secondaryFiles: [.bai]}
  arriba_strand_flag: { type: [ 'null', {type: enum, name: arriba_std, symbols: ["auto", "reverse", "yes"]}] }
  Chimeric_junction: {type: 'File'}

outputs:
  STAR-Fusion_results: {type: 'File?', outputSource: star_fusion_1-10-1/abridged_coding, doc: "STAR fusion detection from chimeric
      reads"}
  arriba_fusion_results: {type: 'File?', outputSource: arriba_fusion_2-2-1/arriba_fusions, doc: "Fusion output from Arriba"}
  arriba_fusion_viz: {type: 'File?', outputSource: arriba_draw_2-2-1/arriba_pdf, doc: "pdf output from Arriba"}
  annofuse_filtered_fusions_tsv: {type: 'File?', outputSource: annofuse/annofuse_filtered_fusions_tsv, doc: "Filtered fusions called
      by annoFuse."}
steps:
  star_fusion_1-10-1:
    run: ../tools/star_fusion_1.10.1_call.cwl
    in:
      Chimeric_junction: star_2-7-10a/chimeric_junctions
      genome_tar: FusionGenome
      output_basename: basename_picker/outname
      genome_untar_path: star_fusion_genome_untar_path
      compress_chimeric_junction: compress_chimeric_junction
    out: [abridged_coding, chimeric_junction_compressed]
  arriba_fusion_2-2-1:
    run: ../tools/arriba_fusion_2.2.1.cwl
    in:
      genome_aligned_bam:
        source: [samtools_sort/sorted_bam, samtools_sort/sorted_bai]
        valueFrom: |
          ${
            var bundle = self[0];
            bundle.secondaryFiles = [self[1]];
            return bundle;
          }
      memory: arriba_memory
      reference_fasta: reference_fasta
      gtf_anno: gtf_anno
      outFileNamePrefix: basename_picker/outname
      arriba_strand_flag: strand_parse/arriba_std
    out: [arriba_fusions]
  arriba_draw_2-2-1:
    run: ../tools/arriba_draw_2.2.1.cwl
    in:
      fusions: arriba_fusion_2-2-1/arriba_fusions
      genome_aligned_bam:
        source: [samtools_sort/sorted_bam, samtools_sort/sorted_bai]
        valueFrom: |
          ${
            var bundle = self[0];
            bundle.secondaryFiles = [self[1]];
            return bundle;
          }
      gtf_anno: gtf_anno
      memory: arriba_memory
    out: [arriba_pdf]
  annofuse:
    run: ../workflow/kfdrc_annoFuse_wf.cwl
    in:
      sample_name: basename_picker/outsample
      FusionGenome: fusion_annotator_ref
      genome_untar_path: star_fusion_genome_untar_path
      rsem_expr_file: rsem/gene_out
      arriba_output_file: arriba_fusion_2-2-1/arriba_fusions
      star_fusion_output_file: star_fusion_1-10-1/abridged_coding
      col_num: annofuse_col_num
      output_basename: basename_picker/outname
    out: [annofuse_filtered_fusions_tsv]
$namespaces:
  sbg: https://sevenbridges.com
hints:
- class: "sbg:maxNumberOfParallelInstances"
  value: 3
