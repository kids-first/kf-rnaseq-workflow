## Gene Expression Abundance Estimation
We used STAR [@doi:10/f4h523] v2.6.1d to align paired-end RNA-seq reads.
This output was used for all subsequent RNA analysis. The reference we used was that of ENSEMBL's GENCODE 27 [@url:https://www.gencodegenes.org/human/release_27.html], "Comprehensive gene annotation."
We used RSEM [@doi:10/cwg8n5] v1.3.1 for transcript- and gene-level quantification.
We also added a second method of quantification using kallisto [@doi:10.1038/nbt.3519] v0.43.1.
This method differs in that it uses pseudoaligments using fastq reads directly to the aforementioned GENCODE 27 reference.

## RNA Fusion Calling
We set up [Arriba v1.1.0](https://github.com/suhrig/arriba/) and STAR-Fusion 1.5.0 [@doi:10.1101/120295] fusion detection tools using CWL on CAVATICA.
For both these tools we used aligned BAM and chimeric SAM files from STAR as inputs and GRCh38_gencode_v27 GTF for gene annotation.
We ran STAR-Fusion with default parameters and annotated all fusion calls with GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz provided in the STAR-fusion release. 
For Arriba, we used a blacklist file (blacklist_hg38_GRCh38_2018-11-04.tsv.gz) from the Arriba release tarballs to remove recurrent fusion artifacts and transcripts present in healthy tissue.