// nextflow run main.nf \
//     --step preprocess \
//     --vcf tests/data/vcf/test_snv.vcf \
//     --bam tests/data/tumor/tumor.bam \
//     --reference tests/data/reference/reference.fasta \
//     --bed tests/data/reference/reference.bed \
//     --outdir results \
//     --coverage 76 \
//     --medianInsert 254 \
//     --picardMetrics tests/data/picard \
//     --model tests/data/models/snv_model.pkl \
//     --modelName ARTIFACT \
//     --resume