cwlVersion: v1.2
class: CommandLineTool
id: star_fusion_1.10.1_gen_reference
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'migbro/star_fusion:1.10.1'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cores)
    ramMin: ${ return inputs.ram * 1000 }
  - class: InitialWorkDirRequirement
    listing: [$(inputs.ctat_source)]


baseCommand: [tar, -I pigz, -xvf]
arguments:
  - position: 2
    shellQuote: false
    valueFrom: >-
      && /STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl
      ${
          var arg = "--dfam_db ";
          for (var i = 0; i < inputs.dfam_db.length; i++) {
              if (inputs.dfam_db[i].nameext=='.hmm'){
                  arg += inputs.dfam_db[i].path;
                  return arg;
              }
          }
      }
  - position: 4
    shellQuote: false
    valueFrom: >-
      && tar -I pigz -cf $(inputs.output_dir).tar.gz $(inputs.output_dir)

inputs:
  ctat_source: { type: File, doc: "Resource file from https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/", inputBinding: { position: 1} }
  genome_fa: { type: File, doc: "Reference fasta file", inputBinding: { position: 3, prefix: "--genome_fa"} }
  fusion_annot_lib: { type: File, doc: "Can be extracted from ctat source, for from https://github.com/FusionAnnotator/CTAT_HumanFusionLib/releases",
  inputBinding: { position: 3, prefix: "--fusion_annot_lib"} }
  annot_filter_rule: {type: File, doc: "target AnnotFilterRule.pm, from https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB",
  inputBinding: { position: 3, prefix: "--annot_filter_rule"} }
  pfam_db: { type: File, doc: "pfam database obtained from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
  inputBinding: { position: 3, prefix: '--pfam_db' } }
  dfam_db: { type: 'File[]',
  doc: "DNA transposable element database (Dfam.hmm), required for repeat masking. Obtain from http://dfam.org/releases/Dfam_3.1/infrastructure/dfamscan/" }
  cores: { type: 'int?', doc: "Num cores to use", default: 16, inputBinding: { position: 3, prefix: '--CPU' } }
  ram: { type: 'int?', doc: "Memory in GB to use", default: 64}
  reference_gtf: { type: File, doc: "gene model definitions." , inputBinding: { position: 2, prefix: '--gtf' } }
  human_gencode_filter: { type: 'string?', doc: "flag for customized prep operations for human/gencode genome and annotation data.",
  default: "--human_gencode_filter" }
  output_dir: { type: string, doc: "Output dirname", inputBinding: { position: 3, prefix: "--output_dir"} }

outputs:
  star_fusion_reference:
    type: File
    outputBinding: 
      glob: '*tar.gz'
