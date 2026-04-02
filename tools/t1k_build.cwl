cwlVersion: v1.2
class: CommandLineTool
id: t1k-build
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.cpu)
  - class: DockerRequirement
    dockerPull: pgc-images.sbgenomics.com/d3b-bixu/t1k:v1.0.5
  - class: NetworkAccess
    networkAccess: true
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      perl /T1K-1.0.5/t1k-build.pl -o .

inputs:
  dat: { type: 'File?', inputBinding: { position: 2, prefix: '-d' }, doc: "EMBL-ENA dat file" }
  seq: { type: 'File?', inputBinding: { position: 2, prefix: '-f' }, doc: "plain gene sequence file" }
  download: { type: 'string?', inputBinding: { position: 2, prefix: '--download' }, doc: "IPD-IMGT/HLA or IPD-KIR or user-specified dat file download link" }
  genome_annot: { type: 'File?', inputBinding: { position: 2, prefix: '-g' }, doc: "genome annotation GTF. can be GZIP." }
  target: { type: 'string?', inputBinding: { position: 2, prefix: '--target' }, doc: "gene name keyword" }
  prefix: { type: 'string?', default: "hla", inputBinding: { position: 2, prefix: '--prefix' }, doc: "file prefix" }
  ignore_partial: { type: 'boolean?', inputBinding: { position: 2, prefix: '--ignore-partial' }, doc: "ignore partial allele at all" }
  cpu: { type: 'int?', doc: "Num processing threads to use", default: 2 }
  ram: { type: 'int?', doc: "Num GB memory to make available", default: 4 }
outputs:
  rna_ref_seqs: { type: 'File', outputBinding: { glob: '*rna_seq.fa' } }
  rna_gene_coords: { type: 'File', outputBinding: { glob: '*rna_coord.fa' } }
  dna_ref_seqs: { type: 'File', outputBinding: { glob: '*dna_seq.fa' } }
  dna_gene_coords: { type: 'File', outputBinding: { glob: '*dna_coord.fa' } }
