cwlVersion: v1.2
class: CommandLineTool
id: t1k_genotyper
doc: "Run T1K genotyper 'The ONE genotyper for Kir and HLA'"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.threads)
    https://platform.illumina.com/rdf/ica/resources:tier: economy
    https://platform.illumina.com/rdf/ica/resources:type: standard
    https://platform.illumina.com/rdf/ica/resources:size: large
  - class: DockerRequirement
    dockerPull: pgc-images.sbgenomics.com/d3b-bixu/t1k:v1.0.5
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      $(inputs.cram_reference != null && inputs.bam != null && inputs.bam.basename.search(/.cram$/) != -1 ? "perl /opt/seq_cache_populate.pl -root .cache/hts-ref " + inputs.cram_reference.path + " && export REF_CACHE=.cache/hts-ref/%2s/%2s/%s &&" : "")
  - position: 10
    shellQuote: false
    valueFrom: >-
      run-t1k
  - position: 20
    shellQuote: false
    valueFrom: >-
      && sed -i '1s/^/gene_name\tnum_diff_alleles\tallele_1\tabundance_1\tquality_1\tallele_2\tabundance_2\tquality_2\tsecondary_alleles\n/' *_genotype.tsv

inputs:
  reference: { type: 'File', inputBinding: { position: 12 , prefix: "-f" }, doc: "reference sequence file" }
  cram_reference: { type: 'File?', secondaryFiles: [{pattern: ".fai", required: true}], doc: "FAI-indexed FASTA file used to compress CRAM input." }
  reads: { type: 'File?', doc: "read1 fastq reads if paired-end", inputBinding: { position: 12 , prefix: "-1" } }
  mates: { type: 'File?', doc: "mates (read2) fastq reads if paired-end", inputBinding: { position: 12 , prefix: "-2" } }
  single_end: { type: 'File?', doc: "fastq reads if single-end", inputBinding: { position: 12 , prefix: "-u" } }
  interleaved: { type: 'File?', doc: "fastq reads if interleaved", inputBinding: { position: 12 , prefix: "-i" } }
  bam: { type: 'File?', secondaryFiles: [{pattern: ".bai", required: false}, {pattern: "^.bai", required: false}, {pattern: ".crai", required: false}, {pattern: "^.crai", required: false}], doc: "Indexed BAM/CRAM/SAM input", inputBinding: { position: 12 , prefix: "-b" } }
  output_basename: { type: 'string?', doc: "Prefix string for output file names. Default inferred from input", inputBinding: { position: 12, prefix: "-o"} }
  preset: {type: ['null', {type: enum, name: preset, symbols: ["hla", "hla-wgs", "kir-wgs", "kir-wes"]}], default: "hla", doc: "If paired-end, read orientation", inputBinding: { position: 12 , prefix: "--preset"} }
  stage: {type: ['null', {type: enum, name: stage, symbols: ["0", "1", "2"]}], inputBinding: { position: 12 , prefix: "--stage" }, doc: "start genotyping on specified stage; 0: start from beginning (candidate read extraction). 1: start from genotype with candidate reads. 2: start from post analysis" }

  gene_coordinates: { type: 'File?', inputBinding: { position: 12 , prefix: "-c" }, doc: "gene coordinate file (required when providing BAM/-b input)" }
  min_align_sim: { type: 'float?', inputBinding: { position: 12 , prefix: "-s" }, doc: "minimum alignment similarity" }
  frac: { type: 'float?', inputBinding: { position: 12 , prefix: "--frac" }, doc: "filter if abundance is less than the frac of dominant allele" }
  cov: { type: 'float?', inputBinding: { position: 12 , prefix: "--cov" }, doc: "filter genes with average coverage less than the specified value" }
  cross_gene_rate: { type: 'float?', inputBinding: { position: 12 , prefix: "--crossGeneRate" }, doc: "the effect from other gene's expression" }
  allele_digit_units: { type: 'int?', inputBinding: { position: 12 , prefix: "--alleleDigitUnits" }, doc: "the number of units in genotyping result" }
  allele_delimiter: { type: 'string?', inputBinding: { position: 12 , prefix: "--alleleDelimiter" }, doc: "the delimiter character for digit unit" }
  allele_whitelist: { type: 'string?', inputBinding: { position: 12 , prefix: "--alleleWhitelist" }, doc: "only consider read aligned to the listed allele series" }
  barcode_bam: { type: 'string?', inputBinding: { position: 12 , prefix: "--barcode" }, doc: "For BAM/-b inputs, name of field for barcode" }
  barcode_fastq: { type: 'File?', inputBinding: { position: 12 , prefix: "--barcode" }, doc: "For -1 -2/-u inputs, file containing barcodes" }
  barcode_range: { type: 'string?', inputBinding: { position: 12 , prefix: "--barcodeRange" }, doc: "start, end(-1 for length-1), strand in a barcode is the true barcode. Format: INT INT CHAR (e.g.: 0 -1 +)" }
  barcode_whitelist: { type: 'File?', inputBinding: { position: 12 , prefix: "--barcodeWhitelist" }, doc: "barcode whitelist file" }
  read1_range: { type: 'string?', inputBinding: { position: 12 , prefix: "--read1Range" }, doc: "start, end(-1 for length-1) in -1/-u files for genomic sequence. Format: INT INT (e.g.: 0 -1)" }
  read2_range: { type: 'string?', inputBinding: { position: 12 , prefix: "--read2Range" }, doc: "start, end(-1 for length-1) in -2 files for genomic sequence. Format: INT INT (e.g.: 0 -1)" }
  mate_id_suffix_len: { type: 'int?', inputBinding: { position: 12 , prefix: "--mateIdSuffixLen" }, doc: "the suffix length in read id for mate" }
  abnormal_unmap_flag: { type: 'boolean?', inputBinding: { position: 12 , prefix: "--abnormalUnmapFlag" }, doc: "Set if the flag in BAM for the unmapped read-pair is nonconcordant" }
  relax_intron_align: { type: 'boolean?', inputBinding: { position: 12 , prefix: "--relaxIntronAlign" }, doc: "Set to allow one more mismatch in intronic alignment" }
  no_extraction: { type: 'boolean?', inputBinding: { position: 12 , prefix: "--noExtraction" }, doc: "Set to directly use the files from provided -1 -2/-u for genotyping thus skipping extraction" }
  skip_post_analysis: { type: 'boolean?', inputBinding: { position: 12 , prefix: "--skipPostAnalysis" }, doc: "Set to skip post analysis and only conduct genotyping" }
  output_read_assignment: { type: 'boolean?', inputBinding: { position: 12 , prefix: "--outputReadAssignment" }, doc: "Set to output the allele assignment for each read to prefix_assign.tsv file" }

  threads: { type: 'int?', doc: "Num processing threads to use", default: 8, inputBinding: { position: 12, prefix: "-t" } }
  ram: { type: 'int?', doc: "Num GB memory to make available", default: 32 }
outputs:
  aligned_fasta: { type: 'File[]', outputBinding: { glob: '*_aligned_*.fa' } }
  allele_tsv: { type: File, outputBinding: { glob: '*_allele.tsv' } }
  allele_vcf: { type: 'File?', outputBinding: { glob: '*_allele.vcf' } }
  candidate_fastqs: { type: 'File[]', outputBinding: { glob: '*_candidate_*.fq' } }
  genotype_tsv: { type: File, outputBinding: { glob: '*_genotype.tsv' } }
  read_assignments: { type: 'File?', outputBinding: { glob: '*_assign.tsv' } }
