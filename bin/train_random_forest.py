#!/usr/bin/env python3

import pickle
from pathlib import Path

import pandas as pd
# from sklearn.model_selection import train_test_split
from imblearn.ensemble import BalancedRandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

pd.options.display.float_format = "{:.2f}".format

NUMERICAL_COL = [
    "VAF",
    "DEPTH",
    "LOG_DEPTH_RATIO",
    "AVG_MQ",
    "AVG_ALT_MQ",
    "AVG_ALT_MATE_MQ",
    "LOG_IS_RATIO",
    "LOG_ALT_IS_RATIO",
    "AVG_EDIT_DIST",
    "AVG_READ_BAL",
    "VARIANT_READS",
    "STRAND_BIAS",
]
SNV_NUMERICAL_COL = [
    "VARIANT_ALLELES",
    "PA_BASE_CHANGE_ERROR",
    "PA_TRINUCLEO_ERROR",
    "BB_BASE_CHANGE_ERROR",
    "BB_TRINUCLEO_ERROR",
]
INDEL_NUMERICAL_COL = [
    "INDEL_LENGTH",
    "MHCOUNT",
    "REPCOUNT",
    "INDEL_COUNT",
]
CATEGORICAL_COL = [
    "5_BASE",
    "3_BASE",
]
SNV_CATEGORICAL_COL = [
    "REF",
    "ALT",
]
INDEL_CATEGORICAL_COL = [
    "INDEL_TYPE",
    "CLASSIFICATION",
]

def get_brfc(categorical_columns, numerical_columns):
    categorical_encoder = OneHotEncoder(handle_unknown="ignore")
    numerical_pipe = Pipeline([
        ("imputer", SimpleImputer(strategy="mean"))
    ])
    
    preprocessing = ColumnTransformer(
        [("cat", categorical_encoder, categorical_columns),
         ("num", numerical_pipe, numerical_columns)])
    
    brfc = Pipeline([
        ("preprocess", preprocessing),
        ("classifier", BalancedRandomForestClassifier(random_state=42, n_estimators=100))
    ])
    
    return brfc

def train_random_forest(
    features_path, label_col, model_name, pretrained_model, mutation_type, outdir
):
    """
    Trains a Random Forest model.

    Args:
        features_path (str): Path to tsv with preprocessed features.
        label_col (str): Name of column with artifact labels.
        model_name (str): Name of the model for labeling outputs.
        model_path (str): Path to the trained model (pickle file).
        outdir (str): Directory to save the output files.

    Returns:
        None
    """
    # Ensure output directory exists
    train_dir = Path(outdir) / "train"
    train_dir.mkdir(parents=True, exist_ok=True)

    features = pd.read_csv(features_path, sep="\t", low_memory=False)
    targets = features[label_col].astype(int)

    numerical_columns = NUMERICAL_COL
    categorical_columns = CATEGORICAL_COL
    if mutation_type == "snvs":
        numerical_columns += SNV_NUMERICAL_COL
        categorical_columns += SNV_CATEGORICAL_COL
    else:
        numerical_columns += INDEL_NUMERICAL_COL
        categorical_columns += INDEL_CATEGORICAL_COL
    
    brfc = None
    if pretrained_model:
        # load pretrained model and train with double the estimators
        brfc = pickle.load(open(pretrained_model, 'rb'))
        n_estimators = brfc.named_steps['classifier'].n_estimators
        brfc.named_steps['classifier'].set_params(warm_start=True, n_estimators=n_estimators*2)
        brfc.fit(features[numerical_columns+categorical_columns], targets)
    else:
        brfc = get_brfc(categorical_columns, numerical_columns)
        brfc.fit(features[numerical_columns+categorical_columns], targets)
    
    # print(
    #     "Train accuracy: %0.3f" % brfc.score(
    #         features[numerical_columns+categorical_columns],
    #         targets
    #     )
    # )
    # print("n_estimators: %i" % brfc.named_steps["classifier"].n_estimators)

    outpath = Path(train_dir) / f"model_{model_name}.pkl"
    pickle.dump(brfc, open(outpath, "wb"))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Train the Random Forest model with FFPE mutations."
    )
    parser.add_argument(
        "--features",
        required=True,
        help="Path to tsv with preprocessed features.",
    )
    parser.add_argument(
        "--label-col",
        required=True,
        help="Column in feature tsv with artifact labels (1 = artifact, 0 = real mutation)."
    )
    parser.add_argument(
        "--model-name", required=True, help="Name of model for output."
    )
    parser.add_argument(
        "--mutation-type", required=True, help="Type of mutation ('snvs' or 'indels')."
    )
    parser.add_argument(
        "--pretrained-model",
        help="Path to pretrained model for combined training (optional).",
        default=None,
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Directory to save the output files.",
    )

    args = parser.parse_args()

    train_random_forest(
        features_path=args.features,
        label_col = args.label_col,
        model_name = args.model_name,
        mutation_type = args.mutation_type,
        pretrained_model=args.pretrained_model,
        outdir=args.outdir,
    )
