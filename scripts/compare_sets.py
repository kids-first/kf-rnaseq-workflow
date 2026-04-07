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
    """Sanitize a string for use in filenames."""
    return re.sub(r'[^A-Za-z0-9_.-]', '_', str(s))

def load_files_with_key(paths):
    dfs = []
    for p in paths:
        df = pd.read_csv(p, sep="\t", dtype=str)
        # ensure breakpoint columns exist
        if LEFT_COL not in df.columns or RIGHT_COL not in df.columns:
            raise ValueError(f"File {p} missing required column(s): {LEFT_COL}, {RIGHT_COL}")
        df = df.copy()
        df["key"] = df[LEFT_COL].astype(str) + "|" + df[RIGHT_COL].astype(str)
        dfs.append(df)
    if not dfs:
        return pd.DataFrame(columns=[SAMPLE_COL, "key"])
    return pd.concat(dfs, ignore_index=True)

def write_unique_rows(out_dir, sample, set_label, rows_df):
    """
    Write rows_df to:
      out_dir / f"{safe(sample)}_{safe(set_label)}.unique.tsv"
    If rows_df is empty, nothing is written.
    """
    if rows_df.empty:
        return
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fname = f"{safe_name(sample)}_{safe_name(set_label)}.unique.tsv"
    out_path = out_dir / fname
    rows_df.to_csv(out_path, sep="\t", index=False)

def main():
    ap = argparse.ArgumentParser(
        description="Compare two sets of TSVs (A vs B) by LeftBreakpoint|RightBreakpoint per Sample, plot counts and optionally dump unique entries."
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
                    help="Output figure path (e.g., overlap.png, overlap.pdf, overlap.svg)")

    ap.add_argument("--dump-unique-dir", default=None,
                    help="Directory where per-sample TSVs unique to A and B will be written (optional)")

    ap.add_argument("--percent", action="store_true",
                    help="Stack the bars by percentage rather than raw count")

    args = ap.parse_args()

    if not args.title:
        args.title = f"{args.label_a} vs {args.label_b} Fusion {args.tool} {'Percent' if args.percent else 'Count'} Comparison"

    A = load_files_with_key(args.file_a)
    B = load_files_with_key(args.file_b)

    samples = sorted(set(A.get(SAMPLE_COL, pd.Series(dtype=str)).dropna()) |
                     set(B.get(SAMPLE_COL, pd.Series(dtype=str)).dropna()))

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

        # counts for plot
        rows.append({
            "Sample": s,
            col_unique_a: len(uniqueA_keys),
            col_shared: len(shared_keys),
            col_unique_b: len(uniqueB_keys),
        })

        # dump unique rows if requested — remove the key column for output
        if args.dump_unique_dir:
            uniqueA_rows = Arows[Arows["key"].isin(uniqueA_keys)].drop(columns=["key"], errors='ignore')
            write_unique_rows(args.dump_unique_dir, s, args.label_a, uniqueA_rows)
            uniqueB_rows = Brows[Brows["key"].isin(uniqueB_keys)].drop(columns=["key"], errors='ignore')
            write_unique_rows(args.dump_unique_dir, s, args.label_b, uniqueB_rows)

    counts = pd.DataFrame(rows)
    plot_df = counts.melt(id_vars="Sample", var_name="Category", value_name="Count")

    # Create stacked bar chart
    pivot = counts.set_index("Sample")[[col_unique_a, col_shared, col_unique_b]]

    if args.percent:
        totals = pivot.sum(axis=1)
        # avoid divide-by-zero; rows with total 0 become all 0
        pivot = pivot.div(totals.where(totals != 0, 1), axis=0) * 100.0
        ylabel = "Percent (%)"
    else:
        ylabel = "Count (by LeftBreakpoint|RightBreakpoint)"

    samples_order = pivot.index.tolist()
    categories = pivot.columns.tolist()  # order: unique A, shared, unique B

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
    ax.set_ylabel("Count (by LeftBreakpoint|RightBreakpoint)")
    ax.set_title(args.title)
    ax.legend()
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    plt.savefig(args.out, dpi=300, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    main()
