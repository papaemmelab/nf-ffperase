#!/usr/bin/env python3
from os.path import join
import math
import shutil
import argparse

import pandas as pd
from pysam import FastaFile
from scipy import stats

from microrep_python3 import mhcaller, repcaller, finalcaller


def parse_args():
    parser = argparse.ArgumentParser(description="Combine annotations to create classifier input.")
    parser.add_argument("--pileup", required=True, help="Variants pileups file")
    parser.add_argument("--picard_preadapter", required=True, help="Picard's pre-adapter metrics file")
    parser.add_argument("--picard_baitbias", required=True, help="Picard's bait bias metrics file")
    parser.add_argument("--outdir", required=True, help="Output directory for results")
    parser.add_argument("--reference", required=True, help="Path to reference FASTA")
    parser.add_argument("--coverage", type=int, required=True, help="Global coverage used for LOG_DEPTH_RATIO")
    parser.add_argument("--median_insert", type=int, default=300, help="Global median insert size used for insert ratio calculations")
    parser.add_argument("--indels", action="store_true", help="If set, run the Indel pipeline branch, else run SNV logic")
    return parser.parse_args()

def calculate_strand_bias_score(mutation):
    """
    Calculate fisher score p-value based on strand information and converts to phred.

    Arguments:
        mutation (object): pileup structure with strand information.

    Returns:
        float: fisher exact score in phred scale
    """
    _, pvalue = stats.fisher_exact(
        [[mutation["FR"], mutation["RR"]], [mutation["FA"], mutation["RA"]]]
    )
    fs = -10 * math.log10(pvalue)
    return abs(fs)

def get_indel_change_and_contexts(row, fasta):
    """Get Indel Change and 5' and 3' contexts."""
    change = max(row["REF"], row["ALT"])[1:]
    bases_offset = 25 + len(change)

    context_5_start = row["START"] - bases_offset
    context_5_end = row["START"]
    context_3_start = row["START"] + len(row["REF"])
    context_3_end = row["START"] + len(row["REF"]) + bases_offset

    # Validate context indexes are in the valid range
    chrom_lengths = dict(zip(fasta.references, fasta.lengths))
    max_chrom_length = chrom_lengths[str(row["CHR"])]

    context_5_start = max(context_5_start, 0)
    context_5_end = max(context_5_end, 0)
    context_3_start = min(context_3_start, max_chrom_length)
    context_3_end = min(context_3_end, max_chrom_length)

    # Fetch the contexts
    context_5 = fasta.fetch(
        reference=str(row["CHR"]), start=context_5_start, end=context_5_end
    )
    context_3 = fasta.fetch(
        reference=str(row["CHR"]), start=context_3_start, end=context_3_end
    )

    annotation = pd.Series([change, context_5, context_3])
    return annotation


def get_indel_length_type(row):
    """Get Indel Type."""
    ref_length = len(str(row["REF"]))
    alt_length = len(str(row["ALT"]))

    indel_length = abs(alt_length - ref_length)
    indel_type = ""
    if ref_length > 1 and alt_length == 1:
        indel_type = "D"
    elif ref_length == 1 and alt_length > 1:
        indel_type = "I"
    else:
        indel_type = "DI"

    return pd.Series([indel_length, indel_type])


SNV_CLASSIFIER_COLUMNS = [
    "CHR", "START", "END", "REF", "ALT", "5_BASE", "3_BASE", "VAF", "STRAND_BIAS", 
    "AVG_BQ", "AVG_ALT_BQ", "AVG_MQ", "AVG_ALT_MQ", "AVG_ALT_MATE_MQ", "LOG_IS_RATIO",
    "LOG_ALT_IS_RATIO", "AVG_EDIT_DIST", "AVG_READ_BAL", "DEPTH", "LOG_DEPTH_RATIO", 
    "VARIANT_READS", "VARIANT_ALLELES", "PA_BASE_CHANGE_ERROR", "PA_TRINUCLEO_ERROR", 
    "BB_BASE_CHANGE_ERROR", "BB_TRINUCLEO_ERROR",
]

INDEL_CLASSIFIER_COLUMNS = [
    "CHR", "START", "END", "REF", "ALT", "CHANGE", "5_BASE", "3_BASE", "VAF", 
    "STRAND_BIAS", "AVG_BQ", "AVG_ALT_BQ", "AVG_MQ", "AVG_ALT_MQ", "AVG_ALT_MATE_MQ", 
    "LOG_IS_RATIO", "LOG_ALT_IS_RATIO", "AVG_EDIT_DIST", "AVG_READ_BAL", "DEPTH", 
    "LOG_DEPTH_RATIO", "VARIANT_READS", "INDEL_LENGTH", "INDEL_TYPE", "MHCOUNT", 
    "REPCOUNT", "CLASSIFICATION", "INDEL_COUNT",
]

def main():
    args = parse_args()

    # Read in the pileup data
    df = pd.read_csv(args.pileup, sep="\t")
    df = df.drop_duplicates(subset=["CHR", "START", "REF", "ALT"])

    # Calculate strand bias
    df["STRAND_BIAS"] = df.apply(calculate_strand_bias_score, axis=1)

    # Calculate Log Depth ration with coverage
    coverage = int(args.coverage)
    df["LOG_DEPTH_RATIO"] = df.apply(
        lambda x: math.log(float(x["DEPTH"]) / coverage, coverage)
        if x["DEPTH"] != 0
        else 0,
        axis=1,
    )

    # Add insert size metrics to mutations using median insert size
    global_insert_size = int(args.median_insert)
    df["LOG_IS_RATIO"] = df.apply(
        lambda x: math.log(x["AVG_IS"] / global_insert_size, global_insert_size)
        if round(x["AVG_IS"], 300) > 0
        else 0,
        axis=1,
    )
    df["LOG_ALT_IS_RATIO"] = df.apply(
        lambda x: math.log(x["AVG_ALT_IS"] / global_insert_size, global_insert_size)
        if round(x["AVG_ALT_IS"], 300) > 0
        else 0,
        axis=1,
    )

    # Open reference FASTA
    ref_fasta = FastaFile(args.reference)

    if args.indels:
        # Indel-specific columns
        df[["INDEL_LENGTH", "INDEL_TYPE"]] = df.apply(
            get_indel_length_type, 
            axis=1
        )

        df[["CHANGE", "CONTEXT_5", "CONTEXT_3"]] = df.apply(
            lambda row: get_indel_change_and_contexts(row, ref_fasta),
            axis=1,
        )

        df[["MHCOUNT", "MH"]] = df.apply(
            lambda row: list(mhcaller(row["CHANGE"], row["CONTEXT_3"])),
            axis=1,
            result_type="expand"
        )

        df[["REPCOUNT", "REPEAT"]] = df.apply(
            lambda row: list(
                repcaller(
                    row["CHANGE"],
                    row["CONTEXT_3"],
                    row["CONTEXT_5"],
                    row["INDEL_LENGTH"]
                )
            ),
            axis=1,
            result_type="expand"
        )

        df["CLASSIFICATION"] = df.apply(
            lambda row: finalcaller(
                row["MHCOUNT"], row["REPCOUNT"] * len(row["REPEAT"]), row["REPEAT"]
            ),
            axis=1
        )

        # Keep only Indel-based columns
        df = df[INDEL_CLASSIFIER_COLUMNS]

    else:
        # SNV-specific columns
        df["5_BASE"] = df.apply(
            lambda x: ref_fasta.fetch(str(x["CHR"]), start=x["START"] - 2, end=x["START"] - 1),
            axis=1
        )
        df["3_BASE"] = df.apply(
            lambda x: ref_fasta.fetch(str(x["CHR"]), start=x["END"], end=x["END"] + 1),
            axis=1
        )

        # Read Picard pre-adapter metrics
        pa_df = pd.read_csv(args.picard_preadapter, sep="\t")

        # base-change error
        pa_cxt = pa_df.set_index(["REF_BASE", "ALT_BASE"])
        df["PA_BASE_CHANGE_ERROR"] = df.apply(
            lambda x: pa_cxt.loc[x["REF"], x["ALT"]].ERROR_RATE.sum(), axis=1
        )

        # tri-nucleotide context error
        pa_tri = pa_df.set_index("CONTEXT")
        df["PA_TRINUCLEO_ERROR"] = df.apply(
            lambda x: pa_tri.loc[x["5_BASE"] + x["REF"] + x["3_BASE"]].ERROR_RATE.sum(),
            axis=1
        )

        # Read Picard bait-bias metrics
        bb_df = pd.read_csv(args.picard_baitbias, sep="\t")

        # base-change error
        bb_cxt = bb_df.set_index(["REF_BASE", "ALT_BASE"])
        df["BB_BASE_CHANGE_ERROR"] = df.apply(
            lambda x: bb_cxt.loc[x["REF"], x["ALT"]].ERROR_RATE.sum(), axis=1
        )

        # tri-nucleotide context error
        bb_cxt = bb_df.set_index("CONTEXT")
        df["BB_TRINUCLEO_ERROR"] = df.apply(
            lambda x: bb_cxt.loc[x["5_BASE"] + x["REF"] + x["3_BASE"]].ERROR_RATE.sum(),
            axis=1
        )
        
        # Keep only SNV-based columns        
        df = df[SNV_CLASSIFIER_COLUMNS]

    # Write out annotated file
    output_path = join(args.outdir, "input_df.tsv")
    df.to_csv(output_path, sep="\t", index=False)

    print(f"[INFO] Done! Annotated results written to {output_path}")


if __name__ == "__main__":
    main()
