cwlVersion: v1.0
class: Workflow
id: kfdrc_flagstat_qc
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  input_bam: 'File[]'

outputs:
  flagstats: {type: 'File[]', outputSource: samtools_flagstat/flagstat}

steps:
  samtools_flagstat:
    run: ../tools/samtools_flagstat.cwl
    in:
      input_bam: input_bam
    scatter: input_bam
    out: [flagstat]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:AWSInstanceType'
    value: c5.9xlarge;ebs-gp2;850
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4