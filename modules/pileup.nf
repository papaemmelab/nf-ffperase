
process split_pileup {
    publishDir "${params.outdirPreprocess}/splits", mode: "copy"

    input:
    path vcf

    output:
    path "split_*.vcf", emit: splitVcfs

    script:
    """
    #!/usr/bin/env python
    import os
    
    ANNOT_LIMIT = 1  # Adjust this limit as needed

    with open('${vcf}', 'r', encoding='utf-8') as f_in:
        lines = list(f_in.readlines())
    header = [l for l in lines if l.startswith('#')]
    records = [l for l in lines if not l.startswith('#')]

    for ix, region in enumerate(range(0, len(records), ANNOT_LIMIT)):
        end = min(len(records), region + ANNOT_LIMIT)
        output_vcf_path = f'split_{ix}.vcf'
        with open(output_vcf_path, 'w', encoding='utf-8') as f_out:
            f_out.write(''.join(header))
            f_out.write(''.join(records[region:end]))
    """.stripIndent()
}

process pileup {
    publishDir "${params.outdirPreprocess}/pileup", mode: "copy"
    
    input:
    tuple path(splitVcf), path(bam), path(bai), path(reference)

    output:
    path "pileup_*.txt", emit: pileupVcfs

    script:
    """
    BASENAME_VCF=\$(basename ${splitVcf} .vcf)
    
    annotate_w_pileup \
        ${bam} \
        ${reference} \
        ${splitVcf} \
        pileup_\${BASENAME_VCF}.txt \
        --mapq ${params.minMapq} \
        --baseq ${params.minBaseq} \
        --depth ${params.minDepth}
    """
}

process merge_pileup {
    publishDir "${params.outdirPreprocess}/pileup", mode: "copy"

    input:
    path pileupVcfs

    output:
    path "pileup.txt", emit: pileupOutput

    script:
    """
    pileup_files=( ${pileupVcfs} )
    
    # Get the header from the first file
    head -n 1 "\${pileup_files[0]}" > pileup.txt

    # Concatenate all files, skipping headers from subsequent files
    for f in "\${pileup_files[@]}"; do
        tail -n +2 "\$f" >> pileup.txt
    done
    """
}