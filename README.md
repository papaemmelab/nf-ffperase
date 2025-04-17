# 🧽 nf-ffperase

[![code formatting][black_badge]][black_base]
[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://github.com/askimed/nf-test)
[![nf-ffperase CI](https://github.com/papaemmelab/nf-ffperase/actions/workflows/ci.yaml/badge.svg)](https://github.com/papaemmelab/nf-ffperase/actions/workflows/ci.yaml)

> [!important]
> **You may use `FFPErase`, the underlying content, and any output therefrom for personal use, academic research and noncommercial purposes only. See [LICENSE](LICENSE) for more details.**

Tool for pre-processing and classifying FFPE artifact mutations, using nextflow.

## Contents

- [Run Pipeline](#-run-pipeline)
  - [Preprocess Variants](#-preprocessing-variants)
  - [Classify Variants](#-classifying-artifacts)
- [Contributing](#contributing)

## 🚀 Run Pipeline

You need [Nextflow](https://www.nextflow.io/docs/latest/install.html) installed.

```bash
nextflow run papaemmelab/nf-ffperase --help
```

### 0. ⚡️ Full pipeline

Default value: `--step full`. It runs both `Preprocessing` and `Classify` steps.

See this example:

```bash
nextflow run papaemmelab/nf-ffperase \
    -r main \
    --step full \
    --vcf {snvs.vcf} \
    --bam {tumor.bam} \
    --reference {grch37.fasta} \
    --bed {grch37.genome.bed} \
    --outdir {results} \
    --coverage {100} \
    --medianInsert {250} \
    --model {trained_models/snvs.pkl} \
    --modelName {name}
```

`nf-ffperase` has 2 steps to classify variants, `preprocess` and `classify`:

1. ✏️ `preprocess` takes an input of a VCF, BAM, median coverage and reference fasta and annotates mutations for classification. This step uses [hileup][hileup] and GATK's [Picard][picard] to calculate necessary metrics.

2. 🔮 `classify` takes an input of preprocessed mutations and a model and generates a classification as real or artifact for each mutation.


### 1. ✏️ Preprocessing Variants

`--step preprocess` runs the following processes:

- **Pileup**: mutations to calculate Variant Allele Frequency (VAF).
- **Picard**: its `CollectSequencingArtifactMetrics` command to calculate estimated error rates at the base change and trinucleotide levels. More details on this calculation can be found [here][csam]. The user has the option to run this during `preprocess` or to optionally pass in a directory with the following output files:
  - *.bait_bias_detail_metrics
  - *.pre_adapter_detail_metricsto compute the necessary features:
- **Annotation**: using pileup and picard's output estimates the following features:
  - `Variant Allele Frequency (VAF)`
  - `Average Base Quality (AVG_BQ)`
  - `Average Mapping Quality (AVG_MQ)`
  - `Number of Variant Reads`
  - `Number of distinct Variant Alleles (>= 2% VAF)`
  - `Strand Bias Fisher Score`

#### Example

```bash
nextflow run papaemmelab/nf-ffperase \
    -r main \
    --step preprocess \
    --vcf {snvs.vcf} \
    --bam {tumor.bam} \
    --reference {grch37.fasta} \
    --bed {grch37.genome.bed} \
    --outdir {results} \
    --coverage {100} \
    --medianInsert {250} \
```

Output is the features, located at: `{outdir}/preprocess/features.tsv`.

#### ⚡️ Optional Speed Improvements
Option `--splitPileup` corresponds to number of mutations to include in each pileup split and is set as default to 1000. `--splitReads` corresponds to number of reads to include within each picard split with a default of 7,500,000. If desired and resources are available, decreasing these will increase the number of split jobs optimizing the pileup and picard processes. Changes to these will impact how much memory is required per job so may require updates in nextflow config.

### 2. 🔮 Classifying Artifacts

`--step classify` takes an input of a model type, corresponding model and classifies preprocessed mutations based on their likelihood of being artifactual. Output should be directly from preprocess step, located in the output directory: `{outdir}/preprocess/features.tsv`.

See this example:

```bash
nextflow run papaemmelab/nf-ffperase \
    -r main \
    --step full \
    --features {results/preprocess/features.tsv} \
    --outdir {results} \
    --model {trained_models/snvs.pkl} \
    --modelName {name}
```

### 3. 🧠 Training/Retraining

`--step train` takes an input of preprocessed mutations and a boolean label column (0: real, 1: artifact), a model name, mutation type, and an optional pretrained model to train a new classifier.

See this example:

```bash
nextflow run papaemmelab/nf-ffperase \
    -r main \
    --step train \
    --features {results/preprocess/features.tsv} \
    --labelCol {column name} \
    --modelName {name} \
    --outdir {results} \
    --mutationType {snvs or indels} \
    --modelPath {trained_models/snvs.pkl} (optional)
```

## Contributing

Contributions are welcome, and they are greatly appreciated, check our [contributing guidelines](.github/CONTRIBUTING.md)!


<!-- References -->
[hileup]: https://github.com/brentp/hileup
[picard]: https://broadinstitute.github.io/picard/
[csam]: https://gatk.broadinstitute.org/hc/en-us/articles/360037429491-CollectSequencingArtifactMetrics-Picard-
[black_badge]: https://img.shields.io/badge/code%20style-black-000000.svg
[black_base]: https://github.com/ambv/black
