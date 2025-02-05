process annotate_variants {
    publishDir "${params.outdirPreprocess}", mode: "copy"

    input:
    path pileupOutput
    path picardPreAdapter
    path picardBaitBias
    path reference
    val mutationType

    output:
    path "input_df.tsv", emit: featuresTsv


    script:
    def indelsOption = mutationType == "indels" ? "--indels" : ""
    """
    annotate_variants.py ${indelsOption} \\
        --pileup ${pileupOutput} \\
        --picard_preadapter ${picardPreAdapter} \\
        --picard_baitbias ${picardBaitBias} \\
        --reference ${reference} \\
        --coverage ${params.coverage} \\
        --median_insert ${params.medianInsert} \\
        --outdir \$PWD
    """.stripIndent()
}