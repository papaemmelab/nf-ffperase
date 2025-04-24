process DOWNLOAD_MODEL {
    input:
    val mutationType

    output:
    path "model.${mutationType}.joblib", emit: model

    script:
    """
    #!/usr/bin/env python3

    from huggingface_hub import hf_hub_download
    import shutil

    print(
        f"Downloading papaemmelab/ffperase model for "
        f"${mutationType} from Hugging Face Hub..."
    )
    downloaded_file = hf_hub_download(
        repo_id="papaemmelab/ffperase",
        filename=f"model.${mutationType}.joblib",
    )
    shutil.copy(downloaded_file, f"model.${mutationType}.joblib")
    """.stripIndent()
}


process CLASSIFY_RANDOM_FOREST {
    publishDir "${params.outdir}", mode: "copy"
    
    input:
    path features
    path model
    val modelName
    val mutationType
    path tsv
    
    output:
    path "classify/classified_df_${mutationType}.tsv", emit: classifiedTsv
    path "classify/annotated.tsv", optional: true, emit: annotatedTsv
    
    script:
    def tsvOption = tsv.name != 'NO_FILE' ? "--annotated-tsv ${tsv}" : ""
    """
    classify_w_random_forest.py ${tsvOption} \\
        --features ${features} \\
        --model ${model} \\
        --model-name ${modelName} \\
        --mutation-type ${mutationType}
    """.stripIndent()
}
