class: CommandLineTool
cwlVersion: v1.0
id: rmats-both-bam
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'xinglab/rmats:v4.1.2'
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.threads)
  - class: InitialWorkDirRequirement
    listing:
    - entryname: sample_1.txt
      entry: |
        $(inputs.sample_1.map(function(e){return e.path}).join())
    - entryname: sample_2.txt
      entry: |
        $(inputs.sample_2 ? inputs.sample_2.map(function(e){return e.path}).join() : '')
baseCommand: []
arguments:
  - position: 3
    shellQuote: false
    valueFrom: >-
      && for i in ./$(inputs.output_directory)/*.txt;
      do cp $i $(inputs.output_directory).`basename $i`; done
inputs:
  gtf_annotation: { type: 'File', inputBinding: { position: 2, prefix: '--gtf' }, doc: "Input gtf annotation file." }
  sample_1: { type: 'File[]', inputBinding: { position: 2, prefix: '--b1', valueFrom: 'sample_1.txt' }, doc: "Input sample 1 bam file." }
  sample_2: { type: 'File[]?', inputBinding: { position: 2, prefix: '--b2', valueFrom: 'sample_2.txt' }, doc: "Input sample 2 bam file." }
  output_directory: { type: 'string', inputBinding: { position: 2, prefix: '--od' }, doc: "Name of output directory" }
  temp_directory: { type: 'string?', default: 'temp',  inputBinding: { position: 2, prefix: '--tmp' }, doc: "Name of temporary directory" }
  read_type:
    type:
    - "null"
    - type: enum
      symbols:
      - paired
      - single
      name: read_type
    inputBinding:
      position: 2
      prefix: '-t'
    doc: "Select one option for input read type either paired or single. Tool default: paired"
  strandedness:
    type:
    - "null"
    - type: enum
      symbols:
        - fr-unstranded
        - fr-firststrand
        - fr-secondstrand
      name: strandedness
    inputBinding:
      position: 2
      prefix: '--libType'
    doc: "Select one option for input strandedness. Tool default: fr-unstranded"
  task:
    type:
    - "null"
    - type: enum
      symbols:
        - prep
        - post
        - both
      name: task
    inputBinding:
      position: 2
      prefix: '--task'
    doc: >-
      Specify which step(s) of rMATS to run. prep: preprocess BAMs and generate
      a .rmats file. post: load .rmats file(s) into memory, detect and count
      alternative splicing events, and calculate P value (if not --statoff).
      both: prep + post. Tool default: both
  read_length: { type: 'int', inputBinding: { position: 2, prefix: '--readLength' }, doc: "Input read length for sample reads." }
  variable_read_length: { type: 'boolean?', inputBinding: { position: 2, prefix: '--variable-read-length' }, doc: "Allow reads with lengths that differ from --readLength to be processed. --readLength will still be used to determine IncFormLen and SkipFormLen.", default: true }
  anchor_length: { type: 'int?', inputBinding: { position: 2, prefix: '--anchorLength' }, doc: "The anchor length. Tool default: 1" }
  tophat_anchor_length: { type: 'int?', inputBinding: { position: 2, prefix: '--tophatAnchor' }, doc: "The 'anchor length' or 'overhang length' used in the aligner. At least 'anchor length' NT must be mapped to each end of a given junction. (Only if using fastq). Tool default: 6." }
  star_indicies: { type: 'Directory?', inputBinding: { position: 2, prefix: '--bi' }, doc: "The directory name of the STAR binary indices (name of the directory that contains the SA file). (Only if using fastq)" }
  cutoff_splice_diff: { type: 'float?', doc: "The cutoff used in the null hypothesis test for differential splicing. Tool default is 0.0001 for 0.01% difference. Valid: 0 <= cutoff < 1. Does not apply to the paired stats model" }
  stat_off: { type: 'boolean?', inputBinding: { position: 2, prefix: '--statoff' }, doc: "Select to skip statistical analysis, either between two groups or on single sample group. 'true' to add this parameter. Tool default: false" }
  paired_stats: { type: 'boolean?', inputBinding: { position: 2, prefix: '--paired-stats' }, doc: "Use the paired stats model" }
  novel_splice_sites: { type: 'boolean?', inputBinding: { position: 2, prefix: '--novelSS' }, doc: "Select for novel splice site detection or unannotated splice sites. 'true' to detect or add this parameter, 'false' to disable denovo detection. Tool Default: false" }
  maximum_exon_length: { type: 'int?', inputBinding: { position: 2, prefix: '--mel' }, doc: "Maximum Exon Length. Only impacts --novelSS behavior. Tool Default: 500" }
  minimum_intron_length: { type: 'int?', inputBinding: { position: 2, prefix: '--mil' }, doc: "Minimum Intron Length. Only impacts --novelSS behavior. Tool Default: 50" }
  allow_clipping: { type: 'boolean?', inputBinding: { position: 2, prefix: '--allow-clipping' }, doc: "Allow alignments with soft or hard clipping to be used." }
  fixed_event_set: { type: 'Directory?', inputBinding: { position: 2, prefix: '--fixed-event-set' }, doc: "A directory containing fromGTF.[AS].txt files to be used instead of detecting a new set of events." }
  threads: { type: 'int?', default: 2, inputBinding: { position: 2, prefix: '--nthread' }, doc: "The number of threads. The optimal number of threads should be equal to the number of CPU cores." }
  statistical_threads: { type: 'int?', inputBinding: { position: 2, prefix: '--tstat' }, doc: "The number of threads for the statistical model. If not set then the value of thread is used" }
  ram: { type: 'int?', default: 4, doc: "GB of RAM to allocate to this task." }
outputs:
  alternative_3_prime_splice_sites_jc: { type: 'File', outputBinding: { glob: '*.A3SS.*JC.txt' } }
  alternative_5_prime_splice_sites_jc: { type: 'File', outputBinding: { glob: '*.A5SS.*JC.txt' } }
  mutually_exclusive_exons_jc: { type: 'File', outputBinding: { glob: '*.MXE.*JC.txt' } }
  retained_introns_jc: { type: 'File', outputBinding: { glob: '*.RI.*JC.txt' } }
  skipped_exons_jc: { type: 'File', outputBinding: { glob: '*.SE.*JC.txt' } }
  temp_read_outcomes: { type: File, outputBinding: { glob: 'temp/*_read_outcomes_by_bam.txt'} }
  summary_file: { type: File, outputBinding: { glob: '*summary.txt' }}
  fromGTF:
    type: 'File[]?'
    outputBinding:
      glob: '*fromGTF*'
      outputEval: |
        ${
          if (inputs.novel_splice_sites) {
            return self;
          }
        }
