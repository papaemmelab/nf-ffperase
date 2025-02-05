#!/usr/bin/env python3
"""
merge_picard.py

Merge multiple Picard artifact metrics (pre-adapter and bait-bias) by:
1) Optionally copying original Picard metric files (removing header).
2) Searching for partial metrics in `picard_outdir/tmpPicard`.
3) Summing relevant columns across all partial files.
4) Computing ERROR_RATE, QSCORE, etc.
5) Removing the temporary directory.

Example usage:
    python collect_picard.py --dir /path/to/picard_metrics
"""
from os.path import join, isdir
from glob import glob

import math
import shutil
import argparse
import pandas as pd

pd.options.display.float_format = "{:.2f}".format


PICARD_CONTEXT_COLS = ["SAMPLE_ALIAS", "LIBRARY", "REF_BASE", "ALT_BASE", "CONTEXT"]

PICARD_PRE_ADAPTER_COLS = [
    "PRO_REF_BASES",
    "PRO_ALT_BASES",
    "CON_REF_BASES",
    "CON_ALT_BASES",
]

PICARD_BAIT_BIAS_COLS = [
    "FWD_CXT_REF_BASES",
    "FWD_CXT_ALT_BASES",
    "REV_CXT_REF_BASES",
    "REV_CXT_ALT_BASES",
]

def get_picard_metrics(picard_dir):
    """
    Merges Picard metrics for each interval:
      - Sums up partial pre_adapter files from tmpPicard.
      - Sums up partial bait_bias files from tmpPicard.
      - Computes ERROR_RATE, QSCORE, etc.
      - Writes final merged files to picard_outdir
      - Removes picard_outdir/tmpPicard

    Output:
      - picard_outdir/pre_adapter_metrics.tsv
      - picard_outdir/bait_bias_metrics.tsv
    """
    aggregate = True
    artifacts = join(picard_dir, "tmpPicard")
    if not isdir(artifacts):
        aggregate = False
        artifacts = picard_dir

    # 1) Merge Pre-adapter artifact
    pre_adapter_files = glob(join(artifacts, "*pre_adapter_detail_metrics"))
    if not pre_adapter_files:
        raise FileNotFoundError(
            f"No *pre_adapter_detail_metrics files found in {artifacts}"
        )
    
    if not aggregate:
        # Do not compute metrics if these are provided
        outfile = pd.read_csv(pre_adapter_files[0], sep="\t", skiprows=6)
    else:
        # Use the first file as a baseline
        df0 = pd.read_csv(pre_adapter_files[0], sep="\t", skiprows=6)
        base_counts = df0[PICARD_PRE_ADAPTER_COLS].copy()
        outfile = df0[PICARD_CONTEXT_COLS].copy()

        # Sum up additional files
        for fpath in pre_adapter_files[1:]:
            dfi = pd.read_csv(fpath, sep="\t", skiprows=6)
            base_counts = base_counts.add(dfi[PICARD_PRE_ADAPTER_COLS])

        # Join contextual columns with the summed columns
        outfile = outfile.join(base_counts)

        # Compute ERROR_RATE and QSCORE
        def compute_pre_error_rate(row):
            numerator = float(row["PRO_ALT_BASES"] - row["CON_ALT_BASES"])
            denominator = float(
                row["PRO_ALT_BASES"]
                + row["PRO_REF_BASES"]
                + row["CON_ALT_BASES"]
                + row["CON_REF_BASES"]
            )
            return round(max(1e-10, numerator / denominator), 6)

        outfile["ERROR_RATE"] = outfile.apply(compute_pre_error_rate, axis=1)
        outfile["QSCORE"] = outfile.apply(
            lambda x: int(-10 * math.log10(x["ERROR_RATE"])) if x["ERROR_RATE"] > 0 else 100,
            axis=1,
        )
    
    # Save pre adapter metrics output
    outfile.to_csv(
        "pre_adapter_metrics.tsv",
        float_format="%.6f",
        sep="\t",
        index=False,
    )

    # 2) Merge BAIT-BIAS ARTIFACT MERGING
    bait_bias_files = glob(join(artifacts, "*bait_bias_detail_metrics"))
    if not bait_bias_files:
        raise FileNotFoundError(
            f"No *bait_bias_detail_metrics files found in {artifacts}"
        )
    
    if not aggregate:
        # Do not compute metrics if these are provided
        outfile = pd.read_csv(bait_bias_files[0], sep="\t", skiprows=6)
    else:
        # Read the first file as a baseline
        df0 = pd.read_csv(bait_bias_files[0], sep="\t", skiprows=6)
        base_counts = df0[PICARD_BAIT_BIAS_COLS].copy()
        outfile = df0[PICARD_CONTEXT_COLS].copy()

        for fpath in bait_bias_files[1:]:
            tmp_df = pd.read_csv(fpath, sep="\t", skiprows=6)
            base_counts = base_counts.add(tmp_df[PICARD_BAIT_BIAS_COLS], fill_value=0)

        outfile = outfile.join(base_counts)

        # Compute columns
        def safe_rate(alt, ref):
            return max(1e-10, float(alt) / float(alt + ref))

        outfile["FWD_ERROR_RATE"] = outfile.apply(
            lambda x: round(safe_rate(x["FWD_CXT_ALT_BASES"], x["FWD_CXT_REF_BASES"]), 6),
            axis=1,
        )
        outfile["REV_ERROR_RATE"] = outfile.apply(
            lambda x: round(safe_rate(x["REV_CXT_ALT_BASES"], x["REV_CXT_REF_BASES"]), 6),
            axis=1,
        )
        outfile["ERROR_RATE"] = outfile.apply(
            lambda x: round(max(1e-10, x["FWD_ERROR_RATE"] - x["REV_ERROR_RATE"]), 6),
            axis=1
        )
        outfile["QSCORE"] = outfile.apply(
            lambda x: int(-10 * math.log10(x["ERROR_RATE"])) if x["ERROR_RATE"] > 0 else 100,
            axis=1,
        )
    
    # Save bait bias metrics output
    outfile.to_csv(
        join("bait_bias_metrics.tsv"),
        float_format="%.6f",
        sep="\t",
        index=False,
    )

    # If tmp picard files were used, clean tmp dir.
    if aggregate:
        shutil.rmtree(join(artifacts, "tmpPicard"), ignore_errors=True)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Consolidate Picard artifact metrics. If a tmpPicard subfolder "
            "exists with multiple files it will consolidate metrics, if not, "
            "it will look for a consolidated one to remove headers."
        )
    )
    parser.add_argument(
        "--dir", required=True, help="Path to search for Picard metrics files.",
    )
    args = parser.parse_args()
    
    print(f"[INFO] Getting Picard metrics...")
    get_picard_metrics(args.dir)
    print("[INFO] Done!")


if __name__ == "__main__":
    main()
