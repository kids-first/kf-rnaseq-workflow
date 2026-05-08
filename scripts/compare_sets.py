#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

SAMPLE_COL = "Sample"
LEFT_COL = "LeftBreakpoint"
RIGHT_COL = "RightBreakpoint"


def safe_name(s):
    """Sanitize a string for safe use in filenames.

    Replaces characters that are not alphanumeric, underscore, dot,
    or hyphen with underscores.

    Args:
        s (Any): Input value to sanitize (will be converted to string).

    Returns:
        str: A filename-safe version of the input string.
    """
    return re.sub(r'[^A-Za-z0-9_.-]', '_', str(s))


def load_files_with_key(paths):
    """Load multiple TSV files and add a breakpoint-pair key column.

    Each input file is read as a tab-separated file with all columns
    treated as strings. A new column named ``key`` is created by
    concatenating ``LeftBreakpoint`` and ``RightBreakpoint`` with
    a pipe character ("|").

    Args:
        paths (Iterable[str | Path]): Paths to TSV files to load.

    Returns:
        pandas.DataFrame: Concatenated DataFrame containing all rows
        from the input files, with an added ``key`` column. If no
        files are provided, an empty DataFrame with ``Sample`` and
        ``key`` columns is returned.

    Raises:
        ValueError: If any input file is missing the required
        ``LeftBreakpoint`` or ``RightBreakpoint`` columns.
    """
    dfs = []
    for p in paths:
        df = pd.read_csv(p, sep="\t", dtype=str)
        # ensure breakpoint columns exist
        if LEFT_COL not in df.columns or RIGHT_COL not in df.columns:
            raise ValueError(
                f"File {p} missing required column(s): {LEFT_COL}, {RIGHT_COL}"
            )
        df = df.copy()
        df["key"] = df[LEFT_COL].str.cat(df[RIGHT_COL], sep="|")
        dfs.append(df)

    if not dfs:
        return pd.DataFrame(columns=[SAMPLE_COL, "key"])

    return pd.concat(dfs, ignore_index=True)


def write_unique_rows(out_dir, sample, set_label, rows_df):
    """Write rows unique to a comparison set to a per-sample TSV file.

    Output files are named as:
        ``<sample>_<set_label>.unique.tsv``

    Filenames are sanitized to ensure filesystem safety. If the
    DataFrame is empty, no file is written.

    Args:
        out_dir (str | Path): Output directory for TSV files.
        sample (str): Sample identifier used in the output filename.
        set_label (str): Label indicating which comparison set
            (e.g., SetA or SetB) the rows belong to.
        rows_df (pandas.DataFrame): DataFrame containing rows
            unique to the given set and sample.

    Returns:
        None
    """
    if rows_df.empty:
        return

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fname = f"{safe_name(sample)}_{safe_name(set_label)}.unique.tsv"
    out_path = out_dir / fname
    rows_df.to_csv(out_path, sep="\t", index=False)


def main():
    """CLI entry point for comparing two sets of TSV files.

    This function:
      1. Parses command-line arguments.
      2. Loads TSV files for two comparison sets (A and B).
      3. Compares breakpoint pairs (LeftBreakpoint|RightBreakpoint)
         per sample.
      4. Computes counts of:
         - Entries unique to set A
         - Entries shared between A and B
         - Entries unique to set B
      5. Generates a stacked bar plot (count or percent).
      6. Optionally writes per-sample TSVs of unique entries.

    Returns:
        None
    """
    ap = argparse.ArgumentParser(
        description=(
            "Compare two sets of TSVs (A vs B) by "
            "LeftBreakpoint|RightBreakpoint per Sample, "
            "plot counts and optionally dump unique entries."
        )
    )
    ap.add_argument("-a", "--file-a", nargs="+", required=True,
                    help="One or more TSV files belonging to set A")
    ap.add_argument("-b", "--file-b", nargs="+", required=True,
                    help="One or more TSV files belonging to set B")

    ap.add_argument("--label-a", default="SetA",
                    help="Label to use for set A in the legend and filenames")
    ap.add_argument("--label-b", default="SetB",
                    help="Label to use for set B in the legend and filenames")

    ap.add_argument("--title",
                    help="Plot title")
    ap.add_argument("--tool", default="Tool",
                    help="Tool in comparison")
    ap.add_argument("--figsize", nargs=2, type=float, default=[10, 4],
                    metavar=("W", "H"), help="Figure size in inches")

    ap.add_argument("-o", "--out", required=True,
                    help="Output figure path (e.g., overlap.png, overlap.pdf)")

    ap.add_argument("--dump-unique-dir", default=None,
                    help=(
                        "Directory where per-sample TSVs unique to A "
                        "and B will be written (optional)"
                    ))

    ap.add_argument("--percent", action="store_true",
                    help="Stack the bars by percentage rather than raw count")

    args = ap.parse_args()

    if not args.title:
        args.title = (
            f"{args.label_a} vs {args.label_b} Fusion {args.tool} "
            f"{'Percent' if args.percent else 'Count'} Comparison"
        )

    A = load_files_with_key(args.file_a)
    B = load_files_with_key(args.file_b)

    samples = sorted(
        set(A[SAMPLE_COL].dropna() if SAMPLE_COL in A else []) |
        set(B[SAMPLE_COL].dropna() if SAMPLE_COL in B else [])
    )

    col_unique_a = f"Unique to {args.label_a}"
    col_unique_b = f"Unique to {args.label_b}"
    col_shared = "Shared"

    rows = []
    for s in samples:
        Arows = A.loc[A[SAMPLE_COL] == s].copy() if not A.empty else pd.DataFrame(columns=A.columns)
        Brows = B.loc[B[SAMPLE_COL] == s].copy() if not B.empty else pd.DataFrame(columns=B.columns)

        Akeys = set(Arows["key"].dropna())
        Bkeys = set(Brows["key"].dropna())

        uniqueA_keys = Akeys - Bkeys
        uniqueB_keys = Bkeys - Akeys
        shared_keys = Akeys & Bkeys

        rows.append({
            "Sample": s,
            col_unique_a: len(uniqueA_keys),
            col_shared: len(shared_keys),
            col_unique_b: len(uniqueB_keys),
        })

        if args.dump_unique_dir:
            uniqueA_rows = Arows[Arows["key"].isin(uniqueA_keys)].drop(columns=["key"], errors="ignore")
            write_unique_rows(args.dump_unique_dir, s, args.label_a, uniqueA_rows)

            uniqueB_rows = Brows[Brows["key"].isin(uniqueB_keys)].drop(columns=["key"], errors="ignore")
            write_unique_rows(args.dump_unique_dir, s, args.label_b, uniqueB_rows)

    counts = pd.DataFrame(rows)

    pivot = counts.set_index("Sample")[[col_unique_a, col_shared, col_unique_b]]

    if args.percent:
        totals = pivot.sum(axis=1)
        pivot = pivot.div(totals.where(totals != 0, 1), axis=0) * 100.0
        ylabel = "Percent (%)"
    else:
        ylabel = "Count (by LeftBreakpoint|RightBreakpoint)"

    samples_order = pivot.index.tolist()
    categories = pivot.columns.tolist()

    plt.figure(figsize=args.figsize)
    ax = plt.gca()
    bottom = None
    colors = sns.color_palette("tab10", n_colors=len(categories))

    for i, cat in enumerate(categories):
        vals = pivot[cat].values
        if bottom is None:
            ax.bar(samples_order, vals, label=cat, color=colors[i])
            bottom = vals.copy()
        else:
            ax.bar(samples_order, vals, bottom=bottom, label=cat, color=colors[i])
            bottom = bottom + vals

    ax.set_xlabel("Sample")
    ax.set_ylabel(ylabel)
    ax.set_title(args.title)
    ax.legend()

    if len(samples_order) > 5:
        plt.xticks(rotation=45, ha="right")

    plt.tight_layout()
    plt.savefig(args.out, dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()