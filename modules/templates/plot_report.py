#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import StrMethodFormatter
from pycirclize import Circos

# Inputs
TSV = "${classifiedTsv}"
MUTATION_TYPE = "${mutationType}"
CHROM_BED = "${workflow.projectDir}/assets/hg19.chrom.bed"
CYTOBAND = "${workflow.projectDir}/assets/cytoBand.txt"

# Constants
COLUMNS = [
    "CHR", 
    "START", 
    "END", 
    "REF", 
    "ALT", 
    "VAF",
    "ARTIFACT_raw_predicts",
    "ARTIFACT_predicts",
]
LABEL_MAP = {False: "Real", True: "Artifact"}

def plot_bars(df, mutation_type, outfile="distributions.png"):
    """Plot numbers of real vs artifacts."""
    main_color = "#000099" # blue
    colors = [main_color, f"{main_color}22"]
    threshold = 0.5
    bins      = 20
    bar_width = 0.9

    df["label"] = df["ARTIFACT_predicts"].map(LABEL_MAP)

    # Figure 1
    fig = plt.figure(figsize=(12, 8))
    gs  = GridSpec(2, 2, figure=fig,
                   width_ratios=[1, 10],
                   height_ratios=[1, 1],
                   hspace=0.3)
    axes = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, :]),
    ]

    # 1A: Stacked bar 
    counts = df['label'].value_counts()
    axes[0].bar("Variants", counts["Real"], 
                color=colors[0], label="Real", 
                width=bar_width)
    axes[0].bar("Variants", counts["Artifact"], bottom=counts["Real"], 
                color=colors[1], label="Artifact",
                width=bar_width)
    axes[0].set_ylabel("Count")
    axes[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    real_count = counts["Real"]
    art_count = counts["Artifact"]
    real_pct = f"{real_count*100/(real_count+art_count):.1f}%"
    art_pct = f"{art_count*100/(real_count+art_count):.1f}%"
    axes[0].text(0, real_count/2, f"{real_count:,}", ha='center', va='center', color='white')
    axes[0].text(0, real_count + art_count/2, f"{art_count:,}", ha='center', va='center', color='black')
    axes[0].set_title(f"Real ({real_pct}) vs Artifact ({art_pct})")
    
    # 1B: Probability Histogram Split 
    prob_col = "ARTIFACT_raw_predicts"
    below = df[df[prob_col] < threshold][prob_col]
    above = df[df[prob_col] >= threshold][prob_col]
    axes[1].hist(below, bins=bins,
                 label=f"Real:      Prob < {threshold}",
                 color=colors[0], alpha=1)
    axes[1].hist(above, bins=bins,
                 label=f"Artifact: Prob ≥ {threshold}",
                 color=colors[0], alpha=0.2)
    axes[1].axvline(threshold, color='gray', linestyle='--', linewidth=1)
    axes[1].set_xlabel("Predicted Probability of Artifact")
    axes[1].set_ylabel("Count")
    axes[1].set_title(f"Artifact Probability Distribution by Threshold ({mutation_type}) ")
    axes[1].legend(loc='upper center')
    
    # 1C: Per‐Chromosome Stacked Bar
    chrom_order = [str(i) for i in range(1, 23)] + ['X', 'Y']
    df['CHR'] = pd.Categorical(df['CHR'], categories=chrom_order, ordered=True)
    
    pct_chr = (
        df.groupby('CHR', observed=False)['label']
          .value_counts()
          .unstack(fill_value=0)
          .reindex(chrom_order)
          .loc[:, ["Real", "Artifact"]]
    )
    
    axes[2].bar(pct_chr.index, pct_chr["Real"],      color=colors[0], label="Real", width=bar_width)
    axes[2].bar(pct_chr.index, pct_chr["Artifact"],  bottom=pct_chr["Real"],
           color=colors[1], label="Artifact", width=bar_width)
    
    for i, chrom in enumerate(pct_chr.index):
        real_pct = pct_chr.loc[chrom, "Real"]
        art_pct  = pct_chr.loc[chrom, "Artifact"]
        axes[2].text(i, real_pct/2,     f"{real_pct:,}", ha='center', va='center', color='white')
        axes[2].text(i, real_pct + art_pct/2, f"{art_pct:,}", ha='center', va='center', color='black')
    
    axes[2].set_ylabel("Variant Counts")
    axes[2].set_title(f"Real vs Artifact by Chromosome ({mutation_type})")
    axes[2].legend()
    
    for ax in axes:
        ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    
    fig.subplots_adjust(
        left=0.05, right=0.95,
        top=0.95, bottom=0.05,
        wspace=0.4, hspace=0.4
    )
    fig.savefig(outfile, dpi=100, format="png", bbox_inches="tight")


def plot_circos_snvs(df, outfile="circos.png"):
    """Create circos of pre and post filtering with snvs track."""
    colors = {
        "mutation": {
            "C>T": "#ff0000",
            "T>C": "#07ee00",
            "C>A": "#4169e1",
            "C>G": "#000000",
            "T>A": "#bebebe",
            "T>G": "#ff68b4",
            # reverse strand
            "G>A": "#ff0000",
            "A>G": "#07ee00",
            "G>T": "#4169e1",
            "G>C": "#000000",
            "A>T": "#bebebe",
            "A>C": "#ff68b4",
        },
        "label": {
            "Real": "#000099", 
            "Artifact": "orange",
        }
    }
    df["mutation"] = df.apply(lambda x: x["REF"] + ">" + x["ALT"], axis=1)
    
    # Datasets for each circos
    subsets = [
        (df, "FFPE Raw"),
        (df[df["label"]=="Artifact"], "FFPE Artifacts"),
        (df[df["label"]=="Real"], "FFPE Filtered"),
    ]
    
    # Create one Matplotlib figure with 3 polar subplots
    fig = plt.figure(figsize=(6*len(subsets), 12), dpi=150)
    
    for i, label in enumerate(colors.keys()):
        for j, (sub_df, title) in enumerate(subsets):
            ax = fig.add_subplot(2, 3, 3*i+j+1, polar=True)
            
            # 1) Initialize Circos
            circos = Circos.initialize_from_bed(CHROM_BED, space=2)
            circos.add_cytoband_tracks((95, 100), CYTOBAND)
            circos.text(title, size=15)
            
            # 2) Scatter track (40–90 radius)
            for sector in circos.sectors:
                sector.text(sector.name.replace("chr", ""), size=10)
                track = sector.add_track((40, 90), r_pad_ratio=0.1)
                track.axis(alpha=0.4)
                chrom = sector.name.replace("chr", "")
                df_chr = sub_df[sub_df["CHR"] == chrom]
                if df_chr.empty:
                    continue
                track.scatter(
                    df_chr["START"].tolist(),
                    df_chr["VAF"].tolist(),
                    s=4, vmin=0.0, vmax=1.0,
                    color=df_chr[label].map(colors[label]).tolist(),
                )
            # 3) Render onto our existing polar axis
            circos.plotfig(ax=ax)
    
    # Create handles for mutation legend
    mutation_handles = [
        mpatches.Patch(color=color, label=mut)
        for mut, color in colors["mutation"].items()
    ]
    
    # Create handles for label legend
    label_handles = [
        mpatches.Patch(color=color, label=lbl)
        for lbl, color in colors["label"].items()
    ]
    
    # Top legend: mutations
    fig.legend(
        handles=mutation_handles[:6],
        title="Mutation Types",
        loc="upper center",
        ncol=6,
        bbox_to_anchor=(0.5, 1.05)
    )
    
    # Bottom legend: Real vs Artifact
    fig.legend(
        handles=label_handles,
        title="Classification Type",
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.5, -0.05)
    )
    
    plt.tight_layout()
    fig.savefig(outfile, dpi=200, format="png", bbox_inches="tight")

def plot_circos_indels(df, outfile="circos.png"):
    """Create circos of pre and post filtering with indels track."""
    colors = {
        "indel": {
            "del": "#cd5b45",
            "ins": "#caf374",
        },
        "label": {
            "Real": "#000099", 
            "Artifact": "orange",
        }
    }
    def classify_indel(ref, alt):
        if len(ref) > len(alt):
            return "del"
        elif len(ref) < len(alt):
            return "ins"
        else:
            return np.nan
    
    df['indel'] = df.apply(lambda x: classify_indel(x["REF"], x["ALT"]), axis=1)
    
    # Datasets for each circos
    subsets = [
        (df, "FFPE Raw"),
        (df[df["label"]=="Artifact"], "FFPE Artifacts"),
        (df[df["label"]=="Real"], "FFPE Filtered"),
    ]
    
    # Create one Matplotlib figure with 3 polar subplots
    fig = plt.figure(figsize=(6*len(subsets), 12), dpi=150)
    
    for i, label in enumerate(colors.keys()):
        for j, (sub_df, title) in enumerate(subsets):
            ax = fig.add_subplot(2, 3, 3*i+j+1, polar=True)
            
            # 1) Initialize Circos
            circos = Circos.initialize_from_bed(CHROM_BED, space=2)
            circos.add_cytoband_tracks((95, 100), CYTOBAND)
            circos.text(title, size=15)
            
            # 2) Scatter track (40–90 radius)
            for sector in circos.sectors:
                sector.text(sector.name.replace("chr", ""), size=10)
                track = sector.add_track((40, 90), r_pad_ratio=0.1)
                track.axis(alpha=0.4)
                chrom = sector.name.replace("chr", "")
                df_chr = sub_df[(sub_df['CHR'] == chrom) & (sub_df['indel'].notnull())]
                if df_chr.empty:
                    continue
                for _, row in df_chr.iterrows():
                    ymin, ymax, offset = 40, 90, 1
                    y =  (ymax - ymin) * row['VAF']
                    y1, y2 = max(ymin, (y-offset) + ymin), min(ymax, (y+offset) + ymin)   
                    x1, x2 = row["START"], row["END"]
                    color = colors[label][row[label]]
                    track.rect(x1, x2, r_lim=(y1, y2), 
                               fc=color, ec=color, lw=1)
    
            # 3) Render onto our existing polar axis
            circos.plotfig(ax=ax)
    
    
    # Create handles for mutation legend
    mutation_handles = [
        mpatches.Patch(color=color, label=mut)
        for mut, color in colors["indel"].items()
    ]
    
    # Create handles for label legend
    label_handles = [
        mpatches.Patch(color=color, label=lbl)
        for lbl, color in colors["label"].items()
    ]
    
    # Top legend: mutations
    fig.legend(
        handles=mutation_handles[:6],
        title="Mutation Types",
        loc="upper center",
        ncol=6,
        bbox_to_anchor=(0.5, 1.05)
    )
    
    # Bottom legend: Real vs Artifact
    fig.legend(
        handles=label_handles,
        title="Classification Type",
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.5, -0.05)
    )
    
    plt.tight_layout()
    fig.savefig(outfile, dpi=200, format="png", bbox_inches="tight")

def plot_circos(df, mutation_type, outfile="circos.png"):
    """Select proper tracks to plot in circos."""
    if mutation_type == "snvs":
        plot_circos_snvs(df, outfile)
    elif mutation_type == "indels":
        plot_circos_indels(df, outfile)


df = pd.read_csv(TSV, sep="\t")[COLUMNS]
df['CHR'] = df['CHR'].astype(str)
plot_bars(df, MUTATION_TYPE)
plot_circos(df, MUTATION_TYPE)
