process SPLIT_INTERVALS {
    input:
    path bam
    path bai
    path bed

    output:
    path "split_reads.bed", emit: splitBed

    script:
    """
    split_bed_by_index \\
        ${bed} \\
        ${bam} \\
        split_reads.bed \\
        -c ${params.splitReads}
    """.stripIndent()
}

process PICARD {
    publishDir "${params.outdirPreprocess}/picard/tmpPicard", mode: "copy"

    input:
    tuple val(bedSplitLine), path(bam), path(bai), path(reference), path(picard)

    output:
    path "picard_*", emit: picardFiles

    script:
    """
    LINE="${bedSplitLine}"
    CHR=\$(echo "\$LINE" | cut -f1)
    START=\$(echo "\$LINE" | cut -f2)
    END=\$(echo "\$LINE" | cut -f3)
    REGION="\${CHR}:\${START}-\${END}"
    OUTFILE="picard_\${CHR}_\${START}_\${END}"

    samtools view -h -M ${bam} \${REGION} \\
        | java -jar ${picard} CollectSequencingArtifactMetrics \\
            I=/dev/stdin \\
            O=\${OUTFILE} \\
            R=${reference} \\
            MINIMUM_MAPPING_QUALITY=${params.minMapq} \\
            MINIMUM_QUALITY_SCORE=${params.minBaseq} \\
            VALIDATION_STRINGENCY=LENIENT
    """.stripIndent()
}


process MERGE_PICARD {
    publishDir "${params.outdirPreprocess}/picard", mode: "copy"

    input:
    path picardFiles

    output:
    path "pre_adapter_metrics.tsv", emit: picardPreAdapter
    path "bait_bias_metrics.tsv", emit: picardBaitBias

    script:
    """
    collect_picard.py --dir ${params.outdirPreprocess}/picard
    """.stripIndent()
}


process COPY_PICARD {
    publishDir "${params.outdirPreprocess}/picard", mode: "copy"

    input:
    path picardInput, stageAs: "picardInput"

    output:
    path "pre_adapter_metrics.tsv", emit: preAdapterMetrics
    path "bait_bias_metrics.tsv", emit: baitBiasMetrics

    script:
    """
    collect_picard.py --dir ${picardInput}
    """.stripIndent()
}
