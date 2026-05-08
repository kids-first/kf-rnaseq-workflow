#!/usr/bin/env python3
"""
Fix or add RG tags based on Illumina-style QNAMEs and update @RG header lines.

Illumina QNAME assumptions
--------------------------
Supported formats:

Pre-CASAVA / legacy:
    <flowcell>:<lane>:...

CASAVA >= 1.8:
    <instrument>:<run>:<flowcell>:<lane>:<tile>:<x>:<y>:<index>

Derived RG ID:
    "<flowcell>_<lane>"

Any QNAME not matching these assumptions raises ValueError.
"""

import argparse
import json
import pysam
from typing import Optional


def first_not_none(*values):
    """Return the first non-None value from a list of values.

    Args:
        *values: Arbitrary positional arguments.

    Returns:
        The first value that is not None, or None if all values are None.
    """
    return next((v for v in values if v is not None), None)


def rg_from_qname(qname: str) -> str:
    """Derive a read group (RG) ID from an Illumina-style QNAME.

    Supported formats include legacy Illumina and CASAVA >= 1.8 read names.

    Args:
        qname: The read query name (QNAME).

    Returns:
        A read group ID of the form "<flowcell>_<lane>".

    Raises:
        ValueError: If the QNAME does not match known Illumina formats.
    """
    parts = qname.split(":")

    # Legacy Illumina: <flowcell>:<lane>:...
    if 2 <= len(parts) < 6:
        flowcell, lane = parts[0], parts[1]
        return f"{flowcell}_{lane}"

    # CASAVA >=1.8:
    # <instrument>:<run>:<flowcell>:<lane>:<tile>:<x>:<y>
    if len(parts) >= 7:
        flowcell, lane = parts[2], parts[3]
        return f"{flowcell}_{lane}"

    raise ValueError(f"Unrecognized Illumina QNAME format: {qname!r}")


def pick_rg_defaults_from_header(header_dict: dict):
    """Extract default SM, LB, and PL values from existing RG header lines.

    The first occurrence of each field across all @RG entries is used.

    Args:
        header_dict: BAM/CRAM header dictionary.

    Returns:
        A tuple (sm, lb, pl), where each element may be None.
    """

    """Pick first SM/LB/PL found across existing @RG header lines."""
    sm = lb = pl = None
    for rg in header_dict.get("RG", []):
        if sm is None and rg.get("SM"):
            sm = rg["SM"]
        if lb is None and rg.get("LB"):
            lb = rg["LB"]
        if pl is None and rg.get("PL"):
            pl = rg["PL"]
        if sm and lb and pl:
            break
    return sm, lb, pl


def set_rg_field(rg: dict, key: str, value: Optional[str], force: bool) -> None:
    """Set a field in an @RG header dictionary entry.

    Args:
        rg: The read group dictionary to update.
        key: Header field name (e.g., "SM", "LB", "PL").
        value: Value to set, or None to skip.
        force: If True, overwrite existing value; otherwise set only if absent.
    """
    if value is None:
        return
    if force:
        rg[key] = value
    else:
        rg.setdefault(key, value)



def ensure_rg_header_line(
    header_dict: dict,
    rg_id: str,
    sm: str,
    lb: str,
    pl: str,
    *,
    force_sm: bool = False,
    force_lb: bool = False,
    force_pl: bool = False,
) -> None:
    """Ensure that an @RG header line exists and is populated consistently.

    If an RG with the given ID exists, fields are updated; otherwise a new
    @RG line is created.

    Args:
        header_dict: BAM/CRAM header dictionary.
        rg_id: Read group ID.
        sm: Sample name value.
        lb: Library value.
        pl: Platform value.
        force_sm: Whether to overwrite existing SM.
        force_lb: Whether to overwrite existing LB.
        force_pl: Whether to overwrite existing PL.
    """
    rg_list = header_dict.setdefault("RG", [])

    for rg in rg_list:
        if rg.get("ID") == rg_id:
            set_rg_field(rg, "SM", sm, force_sm)
            set_rg_field(rg, "LB", lb, force_lb)
            set_rg_field(rg, "PL", pl, force_pl)
            return

    new_rg = {"ID": rg_id}
    if sm:
        new_rg["SM"] = sm
    if lb:
        new_rg["LB"] = lb
    if pl:
        new_rg["PL"] = pl
    rg_list.append(new_rg)


def open_in_bam(path: str, reference: str, threads: int) -> pysam.AlignmentFile:
    """Open a BAM or CRAM file for reading.

    Args:
        path: Path to input BAM or CRAM.
        reference: Reference FASTA path (required for CRAM).
        threads: Number of HTSlib threads.

    Returns:
        A pysam.AlignmentFile opened for reading.

    Raises:
        ValueError: If CRAM input is used without a reference.
    """

    if path.lower().endswith(".cram"):
        if not reference:
            raise ValueError("CRAM input requires --reference")
        return pysam.AlignmentFile(
            path, "rc",
            reference_filename=reference,
            threads=threads,
        )
    return pysam.AlignmentFile(path, "rb", threads=threads)


def open_out_bam(path: str, header: dict, reference: str, threads: int) -> pysam.AlignmentFile:
    """Open a BAM or CRAM file for writing.

    Args:
        path: Output file path.
        header: BAM/CRAM header dictionary.
        reference: Reference FASTA path (required for CRAM).
        threads: Number of HTSlib threads.

    Returns:
        A pysam.AlignmentFile opened for writing.

    Raises:
        ValueError: If CRAM output is requested without a reference.
    """
    if path.lower().endswith(".cram"):
        if not reference:
            raise ValueError("CRAM output requires --reference")
        return pysam.AlignmentFile(
            path, "wc",
            header=header,
            reference_filename=reference,
            threads=threads,
        )
    return pysam.AlignmentFile(path, "wb", header=header, threads=threads)


def format_stats(stats: dict, observed_rg_ids: set) -> str:
    """Format processing statistics for human-readable output.

    Args:
        stats: Dictionary of statistics counters.
        observed_rg_ids: Set of observed read group IDs.

    Returns:
        A formatted string suitable for printing.
    """
    return (
        f"Total reads: {stats['total']}\n"
        f"RG added (missing): {stats['rg_missing_added']}\n"
        f"RG fixed (mismatch): {stats['rg_mismatch_fixed']}\n"
        f"RG already OK: {stats['rg_already_ok']}\n"
        f"QNAME parse failed: {stats['qname_parse_failed']}\n"
        f"Distinct RG IDs observed: {len(observed_rg_ids)}\n"
    )


def write_stats_json(path: str, stats: dict, observed_rg_ids: set) -> None:
    """Write processing statistics to a JSON file.

    Args:
        path: Output JSON file path.
        stats: Dictionary of statistics counters.
        observed_rg_ids: Set of observed read group IDs.
    """
    out = dict(stats)
    out["distinct_rg_ids"] = len(observed_rg_ids)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)


def single_pass_fix_rg(
    in_bam: pysam.AlignmentFile,
    out_path: str,
    header_dict: dict,
    reference: str | None,
    threads: int,
    rg_parser,
    sm: str,
    lb: str,
    pl: str,
    force_sm: bool,
    force_lb: bool,
    force_pl: bool,
    buffer_limit: int,
    *,
    fail_on_unknown_qname: bool = False,
) -> None:
    """Fix or add RG tags in a single streaming pass with bounded buffering.

    This function processes reads from an input BAM/CRAM exactly once while
    lazily discovering read group (RG) IDs and constructing the output header.
    Reads are buffered until RG discovery has stabilized, ensuring that no read
    is written before its RG ID is declared in the header (BAM/CRAM compliance).

    The output header is finalized and written only after:
      - at least ``buffer_limit`` reads have been buffered, and
      - no new RG IDs have been observed for ``buffer_limit`` consecutive reads.

    If RG discovery never stabilizes, all reads are buffered until end-of-file,
    preserving correctness at the cost of higher memory usage.

    Args:
        in_bam: Open pysam AlignmentFile positioned at the start of reads.
        out_path: Output BAM/CRAM path.
        header_dict: Mutable BAM/CRAM header dictionary.
        reference: Reference FASTA path (required for CRAM).
        threads: HTSlib thread count.
        rg_parser: Function mapping QNAME → RG ID.
        sm: Sample name for @RG.
        lb: Library for @RG.
        pl: Platform for @RG.
        force_sm: Whether to overwrite existing SM values.
        force_lb: Whether to overwrite existing LB values.
        force_pl: Whether to overwrite existing PL values.
        buffer_limit: RG stabilization/buffering threshold.
        fail_on_unknown_qname: Whether to error on unknown QNAME formats.

    Raises:
        ValueError: If RG parsing fails and fail_on_unknown_qname is True.
    """
    observed_rg_ids = set()
    read_buffer = []
    out_bam = None
    reads_since_new_rg = 0

    def finalize_header_and_begin_writing():
        """Finalize the output header and begin streaming reads to disk.

        This function marks the transition from *buffering mode* to *write mode*
        in single-pass processing. It finalizes the BAM/CRAM header, opens the
        output file, and flushes all buffered reads in a safe, spec-compliant
        order.
        """
        header_dict["RG"].sort(key=lambda rg: rg.get("ID", ""))
        # out_bam is defined in the enclosing scope and we need to assign to it here rather than create a local variable
        nonlocal out_bam
        out_bam = open_out_bam(out_path, header_dict, reference, threads)
        for r in read_buffer:
            out_bam.write(r)
        read_buffer.clear()

    for read in in_bam.fetch(until_eof=True):
        try:
            rg_id = rg_parser(read.query_name)
        except ValueError:
            if fail_on_unknown_qname:
                raise
            rg_id = None

        # fail fast on late RG discovery
        if out_bam is not None and rg_id is not None and rg_id not in observed_rg_ids:
            raise RuntimeError(
                f"Late discovery of RG {rg_id!r} after header finalization. "
                "Try increasing --rg-buffer-reads or disable --single-pass for this input."
            )

        if rg_id:
            if rg_id not in observed_rg_ids:
                observed_rg_ids.add(rg_id)
                reads_since_new_rg = 0
                ensure_rg_header_line(
                    header_dict,
                    rg_id,
                    sm,
                    lb,
                    pl,
                    force_sm=force_sm,
                    force_lb=force_lb,
                    force_pl=force_pl,
                )
            else:
                reads_since_new_rg += 1

            read.set_tag("RG", rg_id, value_type="Z")

        if out_bam is None:
            read_buffer.append(read)

            if (
                len(read_buffer) >= buffer_limit
                and reads_since_new_rg >= buffer_limit
            ):
                finalize_header_and_begin_writing()
        else:
            out_bam.write(read)

    if out_bam is None:
        finalize_header_and_begin_writing()

    out_bam.close()


def main() -> None:
    """Entry point for RG normalization of BAM/CRAM files.

    This function parses command-line arguments, inspects the input file header,
    determines RG defaults, and dispatches processing in either two-pass or
    single-pass mode.

    Two-pass mode (default):
      1. First pass scans all reads to discover RG IDs and construct a complete
         output header.
      2. Second pass rewrites reads, fixing or adding RG tags and emitting a
         summary of processing statistics.

    Single-pass mode:
      - Uses bounded buffering with RG discovery stabilization.
      - Writes output header only when safe.
      - Skips statistics collection and RG header pruning.
      - Exits early after processing.

    Raises:
        ValueError: If CRAM input/output is requested without a reference, or
            if --fail-on-unknown-qname is set and a QNAME cannot be parsed.
    """
    ap = argparse.ArgumentParser(
        description="Fix/add RG tags based on Illumina QNAMEs and update @RG headers."
    )
    ap.add_argument("-i", "--input", required=True, help="Input BAM/CRAM")
    ap.add_argument("-o", "--output", required=True, help="Output BAM/CRAM")
    ap.add_argument("-r", "--reference", help="Reference FASTA (required for CRAM)")

    ap.add_argument("--sm", help="Value for @RG SM (default: inherit → SAMPLE)")
    ap.add_argument("--lb", help="Value for @RG LB (default: inherit → LIB1)")
    ap.add_argument("--pl", help="Value for @RG PL (default: inherit → ILLUMINA)")

    ap.add_argument("-p", "--hts-threads", type=int, default=6)
    ap.add_argument("--fail-on-unknown-qname", action="store_true")
    ap.add_argument(
        "--prune-rg-header",
        action="store_true",
        help="Remove @RG header lines not observed in reads",
    )
    ap.add_argument(
        "--stats-out",
        help="Write processing statistics to this file (JSON format).",
    )
    ap.add_argument(
        "--single-pass",
        action="store_true",
        help="Use single-pass RG discovery with bounded buffering. Disables stats output"
    )
    ap.add_argument(
        "--rg-buffer-reads",
        type=int,
        default=100000,
        help=(
            "Maximum reads to buffer during RG discovery in single-pass mode. "
            "Increasing this value reduces the chance of late RG discovery at the "
            "cost of higher memory usage."
        ),
    )
    args = ap.parse_args()

    observed_rg_ids = set()
    stats = {
        "total": 0,
        "rg_missing_added": 0,
        "rg_mismatch_fixed": 0,
        "rg_already_ok": 0,
        "qname_parse_failed": 0,
    }

    # ---------- Pass 1: discover RGs and build header ----------
    # --- Header-only open ---
    with open_in_bam(args.input, args.reference, args.hts_threads) as tmp_bam:
        header_dict = tmp_bam.header.to_dict()
        header_dict.setdefault("RG", [])

    hdr_sm, hdr_lb, hdr_pl = pick_rg_defaults_from_header(header_dict)

    sm = first_not_none(args.sm, hdr_sm, "SAMPLE")
    lb = first_not_none(args.lb, hdr_lb, "LIB1")
    pl = first_not_none(args.pl, hdr_pl, "ILLUMINA")

    force_sm = args.sm is not None
    force_lb = args.lb is not None
    force_pl = args.pl is not None

    # --- Read stream open ---
    with open_in_bam(args.input, args.reference, args.hts_threads) as in_bam:
        if args.single_pass:
            with open_in_bam(args.input, args.reference, args.hts_threads) as in_bam:
                single_pass_fix_rg(
                    in_bam=in_bam,
                    out_path=args.output,
                    header_dict=header_dict,
                    reference=args.reference,
                    threads=args.hts_threads,
                    rg_parser=rg_from_qname,
                    sm=sm,
                    lb=lb,
                    pl=pl,
                    force_sm=force_sm,
                    force_lb=force_lb,
                    force_pl=force_pl,
                    buffer_limit=args.rg_buffer_reads,
                    fail_on_unknown_qname=args.fail_on_unknown_qname,
                )
            return

        for read in in_bam.fetch(until_eof=True):
            stats["total"] += 1
            try:
                rg_id = rg_from_qname(read.query_name)
            except ValueError:
                stats["qname_parse_failed"] += 1
                if args.fail_on_unknown_qname:
                    raise
                continue
            observed_rg_ids.add(rg_id)

        if args.prune_rg_header:
            header_dict["RG"] = [
                rg for rg in header_dict["RG"]
                if rg.get("ID") in observed_rg_ids
            ]

        for rg_id in sorted(observed_rg_ids):
            ensure_rg_header_line(
                header_dict, rg_id, sm, lb, pl,
                force_sm=force_sm,
                force_lb=force_lb,
                force_pl=force_pl,
            )

    header_dict["RG"].sort(key=lambda rg: rg.get("ID", ""))

    # ---------- Pass 2: write output ----------
    with open_in_bam(args.input, args.reference, args.hts_threads) as in_bam, \
         open_out_bam(args.output, header_dict, args.reference, args.hts_threads) as out_bam:

        for read in in_bam.fetch(until_eof=True):
            try:
                desired = rg_from_qname(read.query_name)
            except ValueError:
                if args.fail_on_unknown_qname:
                    raise
                out_bam.write(read)
                continue

            if not read.has_tag("RG"):
                read.set_tag("RG", desired, value_type="Z")
                stats["rg_missing_added"] += 1
            else:
                current = read.get_tag("RG")
                if current != desired:
                    read.set_tag("RG", desired, value_type="Z")
                    stats["rg_mismatch_fixed"] += 1
                else:
                    stats["rg_already_ok"] += 1

            out_bam.write(read)

    summary = format_stats(stats, observed_rg_ids)
    print(summary)

    if args.stats_out:
        stats_path = args.stats_out + ".json" if not args.stats_out.lower().endswith(".json") else args.stats_out
        write_stats_json(stats_path, stats, observed_rg_ids)


if __name__ == "__main__":
    main()