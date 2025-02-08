#!/bin/bash

nextflow run main.nf \
    -profile hpc_slurm \
    --step preprocess \
    --vcf tests/data/vcf/test_snv.vcf \
    --bam tests/data/tumor/tumor.bam \
    --reference tests/data/reference/reference.fasta



nextflow run main.nf \
    --step preprocess \
    --vcf /data1/papaemme/isabl/data/analyses/32/64/513264/any2lcc.vcf \
    --bam /data1/papaemme/isabl/data/analyses/56/81/425681/IID_H203539_T02_01_WG01.bam \
    --reference /data1/papaemme/isabl/ref/homo_sapiens/GRCh37d5/genome/gr37.fasta \
    --model /data1/papaemme/isabl/ref/homo_sapiens/GRCh37d5/formalinfixer/indel_model.pkl \
    --modelName ARTIFACT \
    --outdir results2 \
    --coverage 81 \
    --medianInsert 349


nextflow run main.nf \
    --step full \
    --vcf /data1/papaemme/isabl/data/analyses/32/64/513264/any2lcc.vcf \
    --bam /data1/papaemme/isabl/data/analyses/56/81/425681/IID_H203539_T02_01_WG01.bam \
    --reference /data1/papaemme/isabl/ref/homo_sapiens/GRCh37d5/genome/gr37.fasta \
    --model /data1/papaemme/isabl/ref/homo_sapiens/GRCh37d5/formalinfixer/indel_model.pkl \
    --modelName ARTIFACT \
    --outdir results2 \
    --coverage 81 \
    --medianInsert 349 \
    --mutationType indels
