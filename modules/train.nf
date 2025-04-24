process TRAIN_RANDOM_FOREST {
    publishDir "${params.outdir}", mode: "copy"
    
    input:
    path features
    val labelCol
    val modelName
    val mutationType
    val modelPath
    
    output:
    path "train/model_${modelName}.pkl", emit: modelOut
    
    script:
    def modelOption = modelPath != '' ? "--pretrained-model ${modelPath}" : ""
    
    """
    train_random_forest.py ${modelOption} \\
        --features ${features} \\
        --label-col ${labelCol} \\
        --model-name ${modelName} \\
        --mutation-type ${mutationType}
    """.stripIndent()
}