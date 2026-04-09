#!/usr/bin/env python3

import argparse
import pysam


def first_not_none(*values):
    return next((v for v in values if v is not None), None)


def rg_from_qname(qname: str) -> str:
    parts = qname.split(":")
    if len(parts) < 2:
        raise ValueError(f"QNAME not Illumina-like (expected >=2 ':'-fields): {qname!r}")
    flowcell = parts[0] if len(parts) < 6 else parts[2]
    lane = parts[1] if len(parts) < 6 else parts[3]
    return f"{flowcell}_{lane}"


def pick_rg_defaults_from_header(header_dict: dict):
    sm = lb = pl = None
    for rg in header_dict.get("RG", []):
        if sm is None and rg.get("SM"):
            sm = rg.get("SM")
        if lb is None and rg.get("LB"):
            lb = rg.get("LB")
        if pl is None and rg.get("PL"):
            pl = rg.get("PL")
        if sm and lb and pl:
            break
    return sm, lb, pl


def set_rg_field(rg: dict, key: str, value: str | None, force: bool) -> None:
    if value is None:
        return
    if force:
        rg[key] = value
    else:
        rg.setdefault(key, value)


def ensure_rg_header_line(header_dict: dict, rg_id: str, sm: str, lb: str, pl: str,
                          force_sm: bool = False, force_lb: bool = False, force_pl: bool = False) -> None:
    # If RG already exists, update its info
    rg_list = header_dict.setdefault("RG", [])
    for rg in rg_list:
        if rg.get("ID") != rg_id:
            continue

        set_rg_field(rg, "SM", sm, force_sm)
        set_rg_field(rg, "LB", lb, force_lb)
        set_rg_field(rg, "PL", pl, force_pl)
        return

    # If RG does not exist, add it to the list
    new_rg = {"ID": rg_id}
    if sm is not None:
        new_rg["SM"] = sm
    if lb is not None:
        new_rg["LB"] = lb
    if pl is not None:
        new_rg["PL"] = pl
    rg_list.append(new_rg)


def open_in_bam(path, reference, threads):
    if path.lower().endswith(".cram"):
        if not reference:
            raise ValueError("CRAM input requires reference")
        return pysam.AlignmentFile(path, "rc", reference_filename=reference, threads=threads)
    return pysam.AlignmentFile(path, "rb", threads=threads)


def main():
    ap = argparse.ArgumentParser(
        description="Fix/add RG tags based on QNAME and update @RG header lines."
    )
    ap.add_argument("-i", "--input", required=True, help="Input BAM/CRAM")
    ap.add_argument("-o", "--output", required=True, help="Output BAM")
    ap.add_argument("-r", "--reference", default=None, help="Reference FASTA (required for CRAM; recommended for CRAM output)")

    ap.add_argument("--sm", default=None, help="Value for @RG SM (default: inherit from input header, else SAMPLE)")
    ap.add_argument("--lb", default=None, help="Value for @RG LB (default: inherit from input header, else LIB1)")
    ap.add_argument("--pl", default=None, help="Value for @RG PL (default: inherit from input header, else ILLUMINA)")

    ap.add_argument("-p", "--hts_threads", type=int, default=6, help="BAM compression/decompression threads")
    ap.add_argument("--fail-on-unknown-qname", action="store_true",
                    help="If set, abort when a QNAME cannot be parsed to derive RG.")
    args = ap.parse_args()

    observed_rg_ids = set()
    stats = {"total": 0, "rg_missing_added": 0, "rg_mismatch_fixed": 0, "rg_already_ok": 0, "qname_parse_failed": 0}

    # Pass 1: collect RG IDs + build updated header
    with open_in_bam(args.input, args.reference, args.hts_threads) as in_bam:
        # Get header and make sure it has RG fields
        header_dict = in_bam.header.to_dict()
        header_dict.setdefault("RG", [])

        # Get SM, LB, and PL from existing RGs
        hdr_sm, hdr_lb, hdr_pl = pick_rg_defaults_from_header(header_dict)

        # Pick defaults for new RGs
        sm = first_not_none(args.sm, hdr_sm, "SAMPLE")
        lb = first_not_none(args.lb, hdr_lb, "LIB1")
        pl = first_not_none(args.pl, hdr_pl, "ILLUMINA")

        # If SM, LB, and PR are set via args, we will overwrite existing RG fields
        force_sm = args.sm is not None
        force_lb = args.lb is not None
        force_pl = args.pl is not None

        # collect RG IDs from input
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

        # Keep only RG header lines observed in input
        existing_rg_lines = header_dict.get("RG", [])
        header_dict["RG"] = [
            rg for rg in existing_rg_lines
            if rg.get("ID") in observed_rg_ids
        ]

        # Cleanup or add observed RGs to header
        for rg_id in sorted(observed_rg_ids):
            ensure_rg_header_line(
                header_dict, rg_id, sm, lb, pl,
                force_sm=force_sm, force_lb=force_lb, force_pl=force_pl
            )

    # Pass 2: write output
    with open_in_bam(args.input, args.reference, args.hts_threads) as in_bam, \
         pysam.AlignmentFile(args.output, "wb", header=header_dict, threads=args.hts_threads) as out_bam:

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
