# T1K Genotyper
[The ONE genotyper for Kir and HLA](https://github.com/mourisl/T1K/).
Currently implemented as a tool, the [usage guide](https://github.com/mourisl/T1K/tree/v1.0.5#usage) explains all the possible tool options.

## Inputs
### Reference
- `output_basename`: Prefix string for output file names. Default inferred from input
- `reference`: FASTA file containing reference sequences for target contigs
- `gene_coordinates`: FASTA file containing gene coordinate information for target contigs. Required when providing a BAM/CRAM for input reads

### Input Reads

The reads for T1K can be delivered as either BAM/CRAM or FASTQ. The [T1K documentation](https://github.com/mourisl/T1K/tree/v1.0.5#inputoutput) outlines the various ways in which these files can be delivered.

#### Single End FASTQ

There is only one option for providing unmapped single end reads:
- `single_end`: Port for `-u` input, FASTQ reads if single-end (omit if paired-end or interleaved)

#### Paired End FASTQ

Unmapped paired end reads can be provided either as separate reads and mates files or as a single interleaved file:
- Paired files:
   - `reads`: Port for `-1` input, read1 FASTQ reads if paired-end (omit if single end or interleaved)
   - `mates`: Port for `-2` input, read2 FASTQ mates if paired-end (omit if single end or interleaved)
- Interleaved file:
   - `interleaved`: Port for `-i` input, FASTQ reads if interleaved (omit if single end, or read1 and mates are separate files)

#### BAM/CRAM

Finally, there is the option to provide the single or paired end reads through BAM/CRAM files. As mentioned above, providing an aligned file necessitates the `gene_coordinates` input:
- `bam`: Port for `-b` input, provide BAM/CRAMs here (requires `gene_coordinates` and, if providing a CRAM, `cram_reference`)
- `cram_reference`: If providing a CRAM input for `bam`, this input is required

### Presets

T1K can be used for RNA-seq, WES/WXS, or WGS analysis. T1k provides presets for HLA and KIR processing. These presets modify the `-s` and `--relaxIntronAlign` parameters:
- `hla`: HLA genotyping on RNA or WES data
   - `-s 0.97` for genotyper step
   - `-s 0.97` for analyzer step
- `hla-wgs`: HLA genotyping on WGS data
   - `-s 0.97` for genotyper step
   - `-s 0.97` for analyzer step
   - `-s 0.97` for FASTQ extraction step
- `kir-wes`: KIR genotyping on RNA or WES data
   - `--relaxIntronAlign` for genotyper step
   - `--relaxIntronAlign` for analyzer step
- `kir-wgs`: KIR genotyping on WGS data. Sets the following:
   - `-s 0.9` and `--relaxIntronAlign` for genotyper step
   - `-s 0.9` and `--relaxIntronAlign` for analyzer step

For more information, see the [T1K documentation](https://github.com/mourisl/T1K/tree/v1.0.5#preset-parameters).

### Tooling Options

All tool options are available in this tool. The options are too numerous to detail here. See the T1K documentation for information when and how to use these parameters.

### Resource Management
- `threads`: Number CPU cores to make available
- `ram`: GB of memory to make available

## Outputs
- `aligned_fasta`
- `allele_tsv`
- `allele_vcf`
- `candidate_fastqs`
- `genotype_tsv`
- `read_assignments`

Please see the [T1K documentation](https://github.com/mourisl/T1K/tree/v1.0.5#inputoutput) for descriptions of the outputs.
