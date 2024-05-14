#!/bin/bash
# Author: Abhishek Choudhary

set -euo pipefail
IFS=$'\n\t'

input_bam=$1
prefix=$2
threads=$3

suffix=Aligned.toTranscriptome_noIDS.out.bam
output_bam=$prefix.$suffix

samtools view $input_bam -@ $threads | awk '{print $1}' | sort | uniq > readids_all
samtools view $input_bam -@ $threads | awk '$6 ~ "I|D|S"' | awk '{print $1}' | sort | uniq > readids_flagIDS
comm -23 readids_all readids_flagIDS > readids_flagnonIDS

samtools view -@ $threads -N readids_flagnonIDS -hb -o $output_bam $input_bam
