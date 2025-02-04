process classify_random_forest {
    publishDir "${params.outdir}", mode: "copy"
    
    input:
    path tsv
    path model
    val modelName
    val mutationType
    
    output:
    path "classify/classified_df_${mutationType}.tsv", emit: classifiedTsv
    
    script:
    """
    classify_w_random_forest.py \
        --model ${model} \
        --model-name ${modelName} \
        --mutation-type ${mutationType}
    """
}