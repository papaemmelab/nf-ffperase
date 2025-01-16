# 🧽 nf-ffperase

[![code formatting][black_badge]][black_base]
[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://github.com/askimed/nf-test)

> [!important]
> **You may use `FFPErase`, the underlying content, and any output therefrom for personal for academic research and noncommercial purposes only. See [LICENSE](LICENSE) for more details.**

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

`nf-ffperase` has 2 steps, `preprocess` and `classify`:

1. ✏️ `preprocess` takes an input of a VCF, BAM, median coverage and reference fasta and annotates mutations for classification. This step uses [hileup][hileup] and GATK's [Picard][picard] to calculate necessary metrics/

2. 🔮 `classify` takes an input of preprocessed mutations and a model and generates a classification as real or artifact for each mutation.

### 1. ✏️ Preprocessing Variants

`nf-ffperase preprocess` runs the following processes:

- Pileup mutations to compute the necessary features:
  - `Variant Allele Frequency (VAF)`
  - `Average Base Quality (AVG_BQ)`
  - `Average Mapping Quality (AVG_MQ)`
  - `Number of Variant Reads`
  - `Number of distinct Variant Alleles (>= 2% VAF)`
  - `Strand Bias Fisher Score`

- Picard's `CollectSequencingArtifactMetrics` is used to calculate estimated error rates at the base change and trinucleotide levels. More details on this calculation can be found [here][csam]. The user has the option to run this during `preprocess` or to optionally pass in a directory with the following output files:
  - *.bait_bias_detail_metrics
  - *.pre_adapter_detail_metrics

#### Example

```bash
nextflow run papaemmelab/nf-ffperase \
    --step preprocess \
    --outdir {outdir} \
    --ffpe-bam {tumor.bam} \
    --reference {genome.gr37.fasta} \
    --input-vcf {snvs.any2vum.vcf} \
    --coverage 15 \
    --median-insert 254 \
    --min-mapq 0 \
    --min-baseq 0 \
    --min-depth 0
```

### 2. 🔮 Classifying Artifacts

`nf-ffperase classify` takes an input of a model type, corresponding model and classifies preprocessed mutations based on their likelihood of being artifactual. Output should be directly from preprocess step, located in the output directory: `{outdir}/input_df.tsv` or `{outdir}/labeled/input_df.labeled.tsv`

See this example:

```bash
nextflow run papaemmelab/nf-ffperase \
    --step classify \
    --outdir {outdir} \
    --ffpe-bam {tumor.bam} \
    --reference {genome.gr37.fasta} \
    --mode rf \
    --model {model_path} \
    --model_name {model_name}
```

## Contributing

Contributions are welcome, and they are greatly appreciated, check our [contributing guidelines](.github/CONTRIBUTING.md)!


<!-- References -->
[hileup]: https://github.com/brentp/hileup
[picard]: https://broadinstitute.github.io/picard/
[csam]: https://gatk.broadinstitute.org/hc/en-us/articles/360037429491-CollectSequencingArtifactMetrics-Picard-
[black_badge]: https://img.shields.io/badge/code%20style-black-000000.svg
[black_base]: https://github.com/ambv/black
