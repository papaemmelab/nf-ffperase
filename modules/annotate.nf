process annotate_variants {
    publishDir "${params.outdirPreprocess}", mode: "copy"

    input:
    path pileupOutput
    path picardPreAdapter
    path picardBaitBias
    path reference

    output:
    path "features.tsv", emit: featuresTsv


    script:
    """
    annotate_variants.py \\
        --pileup ${pileupOutput} \\
        --picard_preadapter ${picardPreAdapter} \\
        --picard_baitbias ${picardBaitBias} \\
        --reference ${reference} \\
        --coverage ${params.coverage} \\
        --median_insert ${params.medianInsert} \\
        --mutation_type ${params.mutationType} \\
        --outdir \$PWD

    rm -rf \\
        ${params.outdirPreprocess}/splits \\
        ${params.outdirPreprocess}/pileup \\
        ${params.outdirPreprocess}/picard
    """.stripIndent()
}
