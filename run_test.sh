#!/bin/bash

nextflow run main.nf \
    -profile hpc_slurm \
    --step preprocess \
    --vcf tests/data/vcf/test_snv.vcf \
    --bam tests/data/tumor/tumor.bam \
    --reference tests/data/reference/reference.fasta