cwlVersion: v1.2
class: CommandLineTool
id: star_fusion_1-10-1_gen_reference
doc: "Generate star fusion reference. Recommend to run locally as tool will take at least 12 hours to run"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/star:fusion-1.10.1'
  - class: ResourceRequirement
    coresMin: $(inputs.cpu)
    ramMin: $(inputs.ram * 1000)
  - class: InitialWorkDirRequirement
    listing: [$(inputs.ctat_source)]

baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      tar -I pigz -xvf
  - position: 10
    shellQuote: false
    valueFrom: >-
      && gunzip -c $(inputs.pfam_db.path) > Pfam-A.hmm
      && hmmpress Pfam-A.hmm
      && gunzip -c $(inputs.genome_fa.path) > genome.fa
      && gunzip -c $(inputs.reference_gtf.path) > annotations.gtf
  - position: 20
    shellQuote: false
    valueFrom: >-
      && /usr/local/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl --pfam_db Pfam-A.hmm --genome_fa genome.fa --gtf annotations.gtf
  - position: 30
    shellQuote: false
    valueFrom: >-
      $(inputs.resource_manifest ? "&& cp " + inputs.resource_manifest.path + " " + inputs.output_dir : "")
      && tar -I pigz -cf $(inputs.output_dir).tar.gz $(inputs.output_dir)

inputs:
  resource_manifest: { type: 'File?', doc: "Manifest of resources used to build the fusion reference." }
  ctat_source: { type: File, inputBinding: { position: 2 }, doc: "Resource file from https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/" }
  genome_fa: { type: File, doc: "GENCODE primary assembly reference fasta file" }
  fusion_annot_lib: { type: File, inputBinding: { position: 22, prefix: "--fusion_annot_lib"}, doc: "Can be extracted from ctat source, for from https://github.com/FusionAnnotator/CTAT_HumanFusionLib/releases" }
  annot_filter_rule: {type: File, inputBinding: { position: 22, prefix: "--annot_filter_rule"}, doc: "target AnnotFilterRule.pm, from https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB" }
  pfam_db: { type: 'File', doc: "Pfam-A.hmm.gz from from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/" }
  dfam_db: { type: 'File', secondaryFiles: [{"pattern":".h3f", "required": true},{"pattern":".h3i", "required": true},{"pattern":".h3m", "required": true},{"pattern":".h3p", "required": true}], inputBinding: { position: 22, prefix: '--dfam_db' }, doc: "DNA transposable element database (Dfam.hmm), required for repeat masking. Obtain from http://dfam.org/releases/Dfam_3.1/infrastructure/dfamscan/" }
  reference_gtf: { type: File, doc: "gene model definitions." }
  human_gencode_filter: { type: 'string?', default: "--human_gencode_filter", doc: "flag for customized prep operations for human/gencode genome and annotation data." }
  output_dir: { type: 'string', default: "star_fusion", inputBinding: { position: 22, prefix: "--output_dir"}, doc: "Name for prep_genome_lib output directory. Will also serve as basename for TAR.GZ file." }
  cpu: { type: 'int?', default: 16, inputBinding: { position: 22, prefix: '--CPU' }, doc: "CPUs to allocate to this task" }
  ram: { type: 'int?', default: 64, doc: "Memory in GB to allocate to this task" }

outputs:
  star_fusion_reference:
    type: File
    outputBinding: 
      glob: '$(inputs.output_dir).tar.gz'
