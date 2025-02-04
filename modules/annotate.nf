process annotate_variants {
    publishDir "${params.outdirPreprocess}/annotate", mode: "copy"

    input:
    path pileupOutput
    path picardPreAdapter
    path picardBaitBias
    path reference
    val indels

    output:
    path "input_df.tsv", emit: annotatedVariants


    script:
    def indelsOption = indels ? "--indels" : ""
    """
    annotate_variants.py ${indelsOption} \
        --pileup ${pileupOutput} \
        --picard_preadapter ${picardPreAdapter} \
        --picard_baitbias ${picardBaitBias} \
        --reference ${reference} \
        --coverage ${params.coverage} \
        --median_insert ${params.medianInsert} \
        --outdir \$PWD
    """.stripIndent()
}