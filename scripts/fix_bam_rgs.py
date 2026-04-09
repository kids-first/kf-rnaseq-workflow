#!/usr/bin/env python3

import argparse
import pysam


def rg_from_qname(qname: str) -> str:
    parts = qname.split(":")
    if len(parts) < 2:
        raise ValueError(f"QNAME not Illumina-like (expected >=4 ':'-fields): {qname!r}")
    flowcell = parts[0] if len(parts) < 6 else parts[2]
    lane = parts[1] if len(parts) < 6 else parts[3]
    return f"{flowcell}_{lane}"


def pick_rg_defaults_from_header(header_dict: dict):
    """
    If the input header already has @RG lines, pick SM/LB/PL defaults from them.
    We take the first non-empty value we find for each field.
    """
    sm = lb = pl = None
    for rg in header_dict.get("RG", []) or []:
        if sm is None and rg.get("SM"):
            sm = rg.get("SM")
        if lb is None and rg.get("LB"):
            lb = rg.get("LB")
        if pl is None and rg.get("PL"):
            pl = rg.get("PL")
        if sm and lb and pl:
            break
    return sm, lb, pl


def ensure_rg_header_line(header_dict: dict, rg_id: str, sm: str, lb: str, pl: str,
                          force_sm: bool = False, force_lb: bool = False, force_pl: bool = False) -> None:
    """
    Ensure an @RG entry with ID exists. Fill SM/LB/PL if missing.
    If force_* is True, overwrite that field.
    """
    rg_list = header_dict.setdefault("RG", [])
    for rg in rg_list:
        if rg.get("ID") == rg_id:
            if force_sm and sm is not None:
                rg["SM"] = sm
            else:
                if sm is not None:
                    rg.setdefault("SM", sm)

            if force_lb and lb is not None:
                rg["LB"] = lb
            else:
                if lb is not None:
                    rg.setdefault("LB", lb)

            if force_pl and pl is not None:
                rg["PL"] = pl
            else:
                if pl is not None:
                    rg.setdefault("PL", pl)
            return

    new_rg = {"ID": rg_id}
    if sm is not None:
        new_rg["SM"] = sm
    if lb is not None:
        new_rg["LB"] = lb
    if pl is not None:
        new_rg["PL"] = pl
    rg_list.append(new_rg)


def main():
    ap = argparse.ArgumentParser(
        description="Fix/add RG tags based on QNAME and update @RG header lines."
    )
    ap.add_argument("-i", "--input", required=True, help="Input BAM/CRAM")
    ap.add_argument("-o", "--output", required=True, help="Output BAM")
    ap.add_argument("-r", "--reference", default=None, help="Reference FASTA (required for CRAM; recommended for CRAM output)")

    # Use None defaults so we can inherit from header when present
    ap.add_argument("--sm", default=None, help="Value for @RG SM (default: inherit from input header, else SAMPLE)")
    ap.add_argument("--lb", default=None, help="Value for @RG LB (default: inherit from input header, else LIB1)")
    ap.add_argument("--pl", default=None, help="Value for @RG PL (default: inherit from input header, else ILLUMINA)")

    ap.add_argument("-p", "--hts_threads", type=int, default=6, help="BAM compression/decompression threads")
    ap.add_argument("--fail-on-unknown-qname", action="store_true",
                    help="If set, abort when a QNAME cannot be parsed to derive RG.")
    ap.add_argument("--keep-existing-header-rg", action="store_true",
                    help="If set, keep all existing @RG lines; otherwise add missing ones only.")
    args = ap.parse_args()

    if args.input.lower().endswith(".cram"):
        if not args.reference:
            raise ValueError("CRAM input requires reference")
        in_bam = pysam.AlignmentFile(args.input, "rc", reference_filename=args.reference, threads=args.hts_threads)
    else:
        in_bam = pysam.AlignmentFile(args.input, "rb", threads=args.hts_threads)
    header_dict = in_bam.header.to_dict()

    # Derive defaults from header if present, but only for args the user did not set
    hdr_sm, hdr_lb, hdr_pl = pick_rg_defaults_from_header(header_dict)

    sm = args.sm if args.sm is not None else (hdr_sm if hdr_sm is not None else "SAMPLE")
    lb = args.lb if args.lb is not None else (hdr_lb if hdr_lb is not None else "LIB1")
    pl = args.pl if args.pl is not None else (hdr_pl if hdr_pl is not None else "ILLUMINA")

    force_sm = args.sm is not None
    force_lb = args.lb is not None
    force_pl = args.pl is not None

    if not args.keep_existing_header_rg:
        header_dict.setdefault("RG", [])

    observed_rg_ids = set()
    stats = {"total": 0, "rg_missing_added": 0, "rg_mismatch_fixed": 0, "rg_already_ok": 0, "qname_parse_failed": 0}

    try:
        # Pass 1: collect RG IDs
        for read in in_bam.fetch(until_eof=True):
            stats["total"] += 1
            try:
                desired = rg_from_qname(read.query_name)
            except Exception:
                stats["qname_parse_failed"] += 1
                if args.fail_on_unknown_qname:
                    raise
                continue
            observed_rg_ids.add(desired)

        # Update header with observed RG IDs (fill/override SM/LB/PL as requested)
        for rg_id in sorted(observed_rg_ids):
            ensure_rg_header_line(
                header_dict, rg_id, sm, lb, pl,
                force_sm=force_sm, force_lb=force_lb, force_pl=force_pl
            )
        in_bam.close()

        # Pass 2: write
        if args.input.lower().endswith(".cram"):
            in_bam = pysam.AlignmentFile(args.input, "rc", reference_filename=args.reference, threads=args.hts_threads)
        else:
            in_bam = pysam.AlignmentFile(args.input, "rb", threads=args.hts_threads)
        out_bam = pysam.AlignmentFile(args.output, "wb", header=header_dict, threads=args.hts_threads)

        for read in in_bam.fetch(until_eof=True):
            try:
                desired = rg_from_qname(read.query_name)
            except Exception:
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

        out_bam.close()
        in_bam.close()

    finally:
        try:
            in_bam.close()
        except Exception:
            pass

    print(
        "Done.\n"
        f"Total reads: {stats['total']}\n"
        f"RG added (missing): {stats['rg_missing_added']}\n"
        f"RG fixed (mismatch): {stats['rg_mismatch_fixed']}\n"
        f"RG already OK: {stats['rg_already_ok']}\n"
        f"QNAME parse failed: {stats['qname_parse_failed']}\n"
        f"Distinct RG IDs in header: {len(observed_rg_ids)}"
    )


if __name__ == "__main__":
    main()
