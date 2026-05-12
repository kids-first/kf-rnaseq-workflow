cwlVersion: v1.2
class: CommandLineTool
id: download_gencode
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.cpu)
  - class: DockerRequirement
    dockerPull: pgc-images.sbgenomics.com/danmiller/downloader:0.1.0 
  - class: InitialWorkDirRequirement
    listing:
      - entryname: download_gencode.sh
        entry:
          $include: ../scripts/download_gencode.sh
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      /bin/bash download_gencode.sh $(inputs.gencode_version) $(inputs.ctat_resource_version) $(inputs.ctat_fusion_version) $(inputs.hla_version) $(inputs.dfam_version) $(inputs.pfam_version)

inputs:
  gencode_version: { type: 'string', doc: "Version of GENCODE to download from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/" }
  ctat_resource_version: { type: 'string', doc: "Version of CTAT resource SOURCE file to download from https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/" }
  ctat_fusion_version: { type: 'string', doc: "Version of CTAT fusion dat.gz to download from https://github.com/FusionAnnotator/CTAT_HumanFusionLib/releases" }
  hla_version: { type: 'string', doc: "Version of HLA to download from https://github.com/ANHIG/IMGTHLA/releases" }
  dfam_version: { type: 'string', doc: "Version of DFAM to download from https://dfam.org/releases/" }
  pfam_version: { type: 'string', doc: "Version of PFAM to download from https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/" }
  cpu: { type: 'int?', doc: "Num processing threads to use", default: 1 }
  ram: { type: 'int?', doc: "Num GB memory to make available", default: 2 }
outputs:
  manifest: { type: 'File', outputBinding: { glob: "gencode_reference_manifest.txt" }}
  annot_filter_rule: { type: 'File', outputBinding: { glob: "AnnotFilterRule.pm" }}
  gencode_genome: { type: 'File', outputBinding: { glob: "*primary_assembly.genome.fa.gz" }}
  gencode_annotation: { type: 'File', outputBinding: { glob: "*annotation.gtf.gz" }}
  ctat_resource: { type: 'File', outputBinding: { glob: "*source.tar.gz" }}
  ctat_fusion: { type: 'File', outputBinding: { glob: "*dat.gz" }}
  hla: { type: 'File', outputBinding: { glob: "hla.dat" }}
  pfam: { type: 'File', outputBinding: { glob: "Pfam-A.hmm.gz" }}
  dfam:
    type: File
    secondaryFiles:
    - {"pattern": ".h3f", required: true}
    - {"pattern": ".h3i", required: true}
    - {"pattern": ".h3m", required: true}
    - {"pattern": ".h3p", required: true}
    outputBinding:
      glob: "homo_sapiens_dfam.hmm"
