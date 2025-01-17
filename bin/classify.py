from os.path import exists, join, isdir
from os import makedirs
import pickle

import pandas as pd

# Input Params
outdir = ""  #
mutation_type = ""  # "snvs" | "indels"
model = ""  # path to the model
model_name = ""  # name of the model
annotated_tsv = ""  # path to the annotated tsv

classify_dir = join(outdir, "classify")
if not isdir(classify_dir):
    makedirs(classify_dir)

# first check for existing classified df
# input_df_path = output_path
# if not exists(input_df_path):
#     input_df_path = join(outdir, "labeled", "input_df.labeled.tsv")
# if not exists(input_df_path):
#     input_df_path = join(outdir, "input_df.tsv")
# if not exists(input_df_path):
#     raise Exception("No input dataframe available!")

# Load and Validate Model
with open(model, "rb") as model_obj:
    model = pickle.load(model_obj)

if not hasattr(model, "score"):  # need a better validation here
    raise Exception("Not a proper BRFC object.")

# Load Input Dataframe
output_path = join(classify_dir, f"classified_df_{mutation_type}.tsv")
for i in [
    output_path,
    join(outdir, "labeled", "input_df.labeled.tsv"),
    join(outdir, "input_df.tsv"),
]:
    if exists(i):
        input_df_path = i
        break

if not exists(input_df_path):
    raise Exception("No input dataframe available!")

input_df = pd.read_csv(input_df_path, sep="\t", low_memory=False)
cols_to_drop = ["CHR", "START", "END"]
for col in input_df.columns:
    if col == "ARTIFACT" or "predicts" in col:
        cols_to_drop.append(col)

features = input_df.drop(cols_to_drop, axis=1)
if mutation_type == "indels":
    features = features.drop(["REF", "ALT", "CHANGE"], axis=1)

predicts = model.predict(features)
raw_scores = model.predict_proba(features)[:, 1]

input_df[f"{model_name}_raw_predicts"] = raw_scores
input_df[f"{model_name}_predicts"] = predicts
input_df[f"{model_name}_predicts"] = input_df[f"{model_name}_predicts"].astype(bool)
input_df.to_csv(output_path, sep="\t", index=False)

if annotated_tsv:
    annotated_tsv = pd.read_csv(annotated_tsv, sep="\t", comment="#", low_memory=False)
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
    annotated_tsv.to_csv(join(classify_dir, "annotated.tsv"), sep="\t", index=False)
