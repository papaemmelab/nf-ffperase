#!/usr/bin/env python3

import pickle
from pathlib import Path
import pandas as pd

pd.options.display.float_format = "{:.2f}".format

def classify_with_random_forest(
    model_path, model_name, mutation_type, annotated_tsv_path, output_dir
):
    """
    Classifies data using a Random Forest model.

    Args:
        model_path (str): Path to the trained model (pickle file).
        model_name (str): Name of the model for labeling outputs.
        mutation_type (str): Type of mutation ("snvs" or "indels").
        annotated_tsv_path (str): Path to the annotated TSV file.
        output_dir (str, optional): Directory to save the output files.

    Returns:
        None
    """
    # Ensure output directory exists
    classify_dir = Path(output_dir) / "classify"
    classify_dir.mkdir(parents=True, exist_ok=True)

    # Paths for output files
    out_classified_tsv = classify_dir / f"classified_df_{mutation_type}.tsv"

    # Load and validate model
    with open(model_path, "rb") as model_file:
        model = pickle.load(model_file)

    if not hasattr(model, "predict") or not hasattr(model, "predict_proba"):
        raise Exception("Invalid pickle file model: Missing necessary methods.")

    # Locate input dataframe
    input_df_path = None
    potential_paths = [
        out_classified_tsv,
        Path(output_dir) / "labeled" / "input_df.labeled.tsv",
        Path(output_dir) / "input_df.tsv",
    ]
    for path in potential_paths:
        print("Path to check:", path)
        if path.exists():
            input_df_path = path
            break

    if not input_df_path:
        raise FileNotFoundError("No input dataframe available!")

    # Load input dataframe
    input_df = pd.read_csv(input_df_path, sep="\t", low_memory=False)

    # Prepare features
    cols_to_drop = ["CHR", "START", "END"]
    for col in input_df.columns:
        if col == "ARTIFACT" or "predicts" in col:
            cols_to_drop.append(col)

    features = input_df.drop(cols_to_drop, axis=1)
    if mutation_type == "indels":
        features = features.drop(["REF", "ALT", "CHANGE"], axis=1)

    # Perform classification
    predicts = model.predict(features)
    raw_scores = model.predict_proba(features)[:, 1]

    # Add predictions to dataframe
    input_df[f"{model_name}_raw_predicts"] = raw_scores
    input_df[f"{model_name}_predicts"] = predicts.astype(bool)
    input_df.to_csv(out_classified_tsv, sep="\t", index=False)

    print(f"Classified TSV saved at: {out_classified_tsv}")

    # Annotate if annotated_tsv_path is provided
    if annotated_tsv_path:
        out_annotated_tsv = classify_dir / "annotated.tsv"
        annotated_tsv = pd.read_csv(
            annotated_tsv_path, sep="\t", comment="#", low_memory=False
        )
        cols = [
            "CHR",
            "START",
            "REF",
            "ALT",
            f"{model_name}_raw_predicts",
            f"{model_name}_predicts",
        ]
        annotated_tsv = annotated_tsv.merge(
            input_df[cols], how="inner", on=["CHR", "START", "REF", "ALT"]
        )
        annotated_tsv.to_csv(out_annotated_tsv, sep="\t", index=False)
        print(f"Annotated TSV saved at: {out_annotated_tsv}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Classify the FFPE mutations using the Random Forest model."
    )
    parser.add_argument(
        "--model", required=True, help="Path to the trained model (pickle file)."
    )
    parser.add_argument(
        "--model-name", required=True, help="Name of the model for labeling outputs."
    )
    parser.add_argument(
        "--mutation-type", required=True, help='Type of mutation ("snvs" or "indels").'
    )
    parser.add_argument(
        "--annotated-tsv",
        help="Path to the annotated TSV file (optional).",
        default=None,
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory to save the output files.",
    )

    args = parser.parse_args()

    classify_with_random_forest(
        model_path=args.model,
        model_name=args.model_name,
        annotated_tsv_path=args.annotated_tsv,
        mutation_type=args.mutation_type,
        output_dir=args.output_dir,
    )
