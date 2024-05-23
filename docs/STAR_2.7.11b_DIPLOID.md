# Kids First Kids First STAR Diploid Beta

This is an alternative alignment and quantification method currently in **_beta phase_**.
It uses a patient's DNA variant calls to create a "personal genome" (PG) for improved alignment.
It is purported to have fewer multi-mapping/better unique mapping for potentially improved gene and isoform level quantification.
It cannot be used in fusion calling
The STAR Diploid mode has a known bug with an unknown cause manifests as a seg fault for up to 20% of normal sample inputs and 100% of tumor inputs. 

![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

## Introduction
This pipeline runs the following steps:
1. STAR Genome Generate per individual, or will skip if one has been generated already.
    -  Strip existing annotations as defined by the user. Since STAR requires an uncompressed VCF, this helps make the input file size smaller as only variant calls and `FORMAT` info are needed
    -  Filtering steps for input patient DNA variant calls in order to focus on higher quality calls.
    Recommended filtering criteria recommendations are still being established, but various scenarios and suggestions from our and partner institutions can be found in the inputs section.
    - Use a genome FASTA and GTF - recommend FASTA matches input DNA calls - in conjunction with filtered DNA VCF to create PG
1. If needed, convert input BAM reads to FASTQ
1. If needed, Cutadapt to remove any adapters
1. Run STAR aligner 
    - Use PG as refs
    - Align input reads
1. Run custom tool to filter BAM for RSEM. Removes indels and soft-clipped reads
1. RSEM quantification

### Cutadapt
[Cutadapt v3.4](https://github.com/marcelm/cutadapt) Cut adapter sequences from raw reads if needed.
### STAR v2.7.11b_alpha_2024-03-29
[STAR v2.7.11b_alpha_2024-03-29](https://doi.org/f4h523) RNA-Seq raw data alignment.
### [RSEM](docs/RSEM_1.3.1.md)

## INPUTS
A brief note - many references and filtering requirements take a bit of up front work. In the [Filtering and Input Appendix](#filtering-and-input-appendix) section, we try to lay this out. The workflow has more inputs that this README will not cover, but can be tweaked by an advanced user.
### Common Required
These are inputs used in multiple steps
 - `output_basename`: String to prepend to results files from STAR and RSEM
### STAR Genome Generate
If a pre-existing PG does not exist, need the following inputs to create:
 - `input_vcf`: Currently not standardized, matched DNA variant calls. For INCLUDE dataset, trio calls were used when available, otherwise singe sample genotyping - both from GATK workflows
 - `strip_info`: Given that input vcf needs to be uncompressed, stripping `INFO` is a good way to reduce file size. Current recommended strip based on typical KF runs: `INFO/CLNDISDB,INFO/CLNDISDBINCL,INFO/CLNDN,INFO/CLNDNINCL,INFO/CLNHGVS,INFO/CLNREVSTAT,INFO/CLNSIG,INFO/CLNSIGCONF,INFO/CLNSIGINCL,INFO/CLNVC,INFO/CLNVCSO,INFO/CLNVI,INFO/CSQ,INFO/ClippingRankSum,INFO/DB,INFO/DP,INFO/DS,INFO/END,INFO/ExcessHet,INFO/FS,INFO/HaplotypeScore,INFO/InbreedingCoeff,INFO/Intervar,INFO/Intervar_STATUS,INFO/MLEAC,INFO/MLEAF,INFO/MQ,INFO/MQRankSum,INFO/NEGATIVE_TRAIN_SITE,INFO/OLD_VARIANT,INFO/POSITIVE_TRAIN_SITE,INFO/QD,INFO/RAW_MQ,INFO/ReadPosRankSum,INFO/SOR,INFO/VQSLOD,INFO/culprit,INFO/gnomad_3_1_1_AC,INFO/gnomad_3_1_1_AC_controls_and_biobanks,INFO/gnomad_3_1_1_AC_popmax,INFO/gnomad_3_1_1_AF,INFO/gnomad_3_1_1_AF_controls_and_biobanks,INFO/gnomad_3_1_1_AF_non_cancer,INFO/gnomad_3_1_1_AF_popmax,INFO/gnomad_3_1_1_AN,INFO/gnomad_3_1_1_AN_controls_and_biobanks,INFO/gnomad_3_1_1_AN_popmax,INFO/gnomad_3_1_1_nhomalt,INFO/gnomad_3_1_1_nhomalt_popmax,INFO/gnomad_3_1_1_primate_ai_score,INFO/gnomad_3_1_1_splice_ai_consequence`
  - `include_expression`: Filters DNA vcf for high quality variants. Current recommended:
    - Trio called VCF: `STRLEN(REF)<=50 && STRLEN(ALT)<=50 && FILTER="PASS" && GT="alt"`
    - Single sample VCF:  `STRLEN(REF)<=50 && STRLEN(ALT)<=50 && FILTER="PASS"`
 - `subtract_bed`: Recommend to filter regions from repeat and low complexity regions. Recommend obtaining repeat-masker bed file from UCSC, run bedtools sort + merge to simplify. Removes variant calls from `input_vcf` from notoriously difficult regions
 - `vcf_sample_name`: **If input is trio**, provide the patient sample name to ensure desired `include_expression` is applied to the specific patient
 - `genome_fa`: Should match input used for DNA. For KF/INCLUDE, recommend `Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta`.
 - `genomeTransformType`: `Diploid`, set by default
 - `gtf`: Recommend `PRI` assembly from [GENCODE version 45](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.primary_assembly.annotation.gtf.gz) for CFDE
 - `sjdbOverhang`: Default is 100. Normally fine as-is, but for PG should probably just set it to read length minus 1
 ### STAR aligner
  - `reads1`: BAM/CRAM/FASTQ reads input. If BAM/CRAM workflow will convert to FASTQ. If FASTQ, read 1 file would go here
  - `reads2`: If FASTQ and paired end, mates file, aka read 2 goes here
  - `cram_reference`: If input reads are CRAM, provide alignment FASTA reference used for it here
  - `outSAMattrRGline`: Output alignment read group. With **_tabs separating the tags_**, format is: ID:sample_name LB:aliquot_id PL:platform SM:BSID for example ID:7316-242 LB:750189 PL:ILLUMINA SM:BS_W72364MN
  - `genomeDir`: If pre-built tar-gzipped PG exists, provide here to skip STAR Genome step
### RSEM
 - `wf_strand_param`: Strandedness of input reads. Default is `rf-stranded`. Use 'default' for unstranded/auto, 'rf-stranded' if read1 in the fastq read pairs is reverse complement to the transcript, 'fr-stranded' if read1 same sense as transcript
 - `RSEMgenome`: RSEM reference tar ball.
## OUTPUTS
 - `star_ref`: If existing `genomeDir` tar ball was not, workflow will have created and provided the patient's PG. This can be re-used in the event a user wants to align another sample frm the patient and/or try different aligner parameters
 - `debug_log`: Log output from STAR GEnome Generate 
 - `STAR_sorted_genomic_cram`: Aligned reads to genome in CRAM format
 - `STAR_transcriptome_bam`: Typically not kept as it's seldom re-used, given that this is in beta phase, we'll keep this
 - `STAR_gene_count`: Gene counts from STAR
 - `STAR_junctions_out`: STAR junctions file
 - `STAR_final_log`: STAR metrics log file of unique, multi-mapping, unmapped, and chimeric reads
 - `RSEM_isoform` RSEM isoform expression estimates
 - `RSEM_gene`: RSEM gene expression estimates

## Filtering and Input Appendix
### Input creation
 - `genome_fa` was generated by first downloading the [complete hg38 FASTA](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta) from Broad, then creating a chr1-22,X,Y,M chromosome list, and running the following commands to get the FASTA file and index: 
   ```sh
   samtools faidx Homo_sapiens_assembly38.fasta -r chr_list.txt > Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
   samtools faidx Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
   ```
 - `subtract_bed` was generated by navigating to https://genome.ucsc.edu/cgi-bin/hgTables, then set the following:
   - **group**: `Repeats`
   - **track**: `RepeatMasker`
   - **output format**: `BED`
   - Result was then piped to `bedtools sort | bedtools merge | gzip > rpt_merge_sort.bed.gz` 
### VCF filtering
#### Trio VCFS
Two main sources of soft `FILTER` are added during creation using [this workflow](https://github.com/kids-first/kf-jointgenotyping-workflow/blob/v2.4.0/README.md):
 - Adding tranches during VQSR
 - Adding `lowGQ` filter `GQ < 20.0`
These variants get removed during bcftools step of filtering on `PASS`
### Single sample VCFS
These have even more `FILTER` added that are dropped on `PASS`:
 - Common
   - ExcessHet: `ExcessHet > 54.69`
   - QD2: `QD < 2.0`
   - QUAL30: `QUAL < 30.0`
 - SNP only:
   - SOR3_SNP: `SOR > 3.0`
   - FS60_SNP: `FS > 60.0`
   - MQ40_SNP: `MQ < 40.0`
   - MQRankSum-12.5_SNP: `MQRankSum < -12.5`
   - ReadPosRankSum-8_SNP: `ReadPosRankSum < -8.0`
 - INDEL only:
   - QD2: `QD < 2.0`
   - QUAL30: `QUAL < 30.0`
   - FS200_INDEL: `FS > 200.0`
   - ReadPosRankSum-20_INDEL: `ReadPosRankSum < -20.0`
### Broad joint cohort calls
Provided by our collaborators at Broad, these are recommended for joint cohort calls:
 - Genotype-quality (GQ): filter out genotypes with a quality below this threshold. Default: 20 - 
 - Allelic-balance: filter out genotypes with allelic balance outside of [1-{allelic-balance}, {allelic-balance}]. Default: 0.8 - 
 - Missingness-threshold: filter out variants with missingness above this threshold. Default: 0.02 - 
 - HWE N-individuals: minimum number of individuals within a subpopulation to perform HWE test. Default: 100 - 
 - HWE: the autosomal (or non-pseudoautosomal X in females only) site failed Hardy-Weinberg equilibrium expectations in (sub)population(s). Minimum acceptable p-value for HWE: Default: 1e-8 - 
 - hmiss pval-threshold: minimum acceptable p-value for hmiss: Default: 1e-8 - 
 - Sample-blacklist: list of samples to exclude. - 
 - ExcessHet (excess of heterozygosity):filter out variants > 54.69 - 
 - LCR (low-complexity regions): filter out variants in low-complexity regions - 
 - Monomorphic alleles: filter out monomorphic alleles (0/0) - 
 - VQSR: Variant Quality Score Recalibration. It is a sophisticated filtering technique applied on the variant callset that uses machine learning to model the technical profile of variants in a training set and uses that to filter out probable artifacts from the callset. We remove variants with VQSR filter. - 
 - Use LCR-hg38-noHLA.interval_list to filter out calls in those regions