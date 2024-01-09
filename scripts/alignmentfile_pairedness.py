#!/usr/bin/env python3

import pysam
import argparse
import pathlib

def parse_command() -> argparse.Namespace:
    """Function to parse the command line input
    Args:
        None
    Return:
        Namespace: returns the args as a standard Namespace object
    """
    parser = argparse.ArgumentParser(description='Determine if input BAM,CRAM,SAM file is paried or single end')
    parser.add_argument('--input_reads',
        required=True,
        help="Path to the BAM,CRAM,SAM file. If providing a CRAM, you must also provide an input_reference.")
    parser.add_argument('--input_reference',
        help="For CRAM only, provide the reference file used when making the input_reads")
    parser.add_argument('--max_reads',
        default=200_000,
        help="The max number of reads to examine to make PAIRED/SINGLE determination. default=%default")
    parser.add_argument('--threads',
        default=1,
        help="For BAM/CRAM decompression, provide the number of threads. default=%default")
    args = parser.parse_args()
    exit = False
    reasons = []
    if not args.input_reads.endswith(('.bam','.cram','.sam')):
        exit = True
        reasons.append(f"input_reads must be a BAM, CRAM, or SAM file")
    if args.input_reads.endswith('.cram') and args.input_reference is None:
        exit = True
        reasons.append(f"input_reference is required when providing a CRAM file")
    if not pathlib.Path(args.input_reads).is_file():
        exit = True
        reasons.append(f"input_reads is not a file")
    if args.input_reference is not None and not pathlib.Path(args.input_reference).is_file():
        exit = True
        reasons.append(f"input_reference is not a file")
    if exit:
        raise Exception("{}{}".format(chr(10),chr(10).join(reasons)))
    else:
        return args

def count_reads(insam: pysam.AlignmentFile, max_count: int) -> tuple[int, ...]:
    """
    Count paired in single reads from an alignment file up to max_count reads.
    Args:
      insam: (pysam.AlignmentFile) File from which to count reads
      max_count: (int) Maximum number of reads to count
    Return:
      tuple[int, ...]: returns tuple of counts paired, single, and total
    """
    total_count = paired_count = single_count = 0
    for aligned_read in insam:
        if aligned_read.is_paired:
            paired_count += 1
        else:
            single_count += 1
        total_count += 1
        if total_count >= max_count: break
    return paired_count, single_count, total_count

def main():
    args = parse_command()
    samfile = pysam.AlignmentFile(args.input_reads, 'r', reference_filename=args.input_reference)
    paired_count, single_count, total_count = count_reads(samfile, args.max_reads)
    if paired_count/total_count > 0.9:
        print("ReadType:PAIRED")
    elif paired_count/total_count < 0.1:
        print("ReadType:SINGLE")
    else:
        print("ReadType:MIXED")

if __name__ == '__main__':
    main()
